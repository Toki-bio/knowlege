# Input and output formats

All line numbers at commit `3c28777`.

## Inputs accepted

| Format | Reader | Notes / quirks |
|---|---|---|
| FASTA / FASTQ, plain or gzip, reference and query | `kseq` via `bseq.c` (`mm_bseq_open` `bseq.c:38`, `mm_bseq_read3` `bseq.c:80-119`) | Query read in mini-batches of `-K` bases (`bseq.c:99-100`); in fragment mode the batch is extended so paired records stay together (`bseq.c:101-109`). Non-ACGT bases map to code 4 (`seq_nt4_table`, `sketch.c:9-26`) and break k-mers (`sketch.c:115`). `-` = stdin, never treated as an index (`index.c:577`). |
| `.mmi` prebuilt index | `mm_idx_is_idx` magic check (`index.c:571-593`), `mm_idx_load` (`index.c:516-569`) | `k, w, HPC` from the index override the command line with a warning (`index.c:626-627`). Multi-part indexes are read one part per `mm_idx_reader_read` (`index.c:621-635`). |
| Two query files (`-x sr`) or `--frag=yes` | `mm_bseq_read_frag2` (`map.c:551`) | At most two query files in `sr` mode (`main.c:421-424`); fragments grouped by name (`mm_qname_same`, `map.c:567`). |
| `--junc-bed`, `-j`, `--pass1`, `--spsc` BED/score files | `mm_idx_bed_read` `index.c:796`, `mm_idx_jjump_read` `index.c:903`, `mm_idx_spsc_read2` `index.c:970` | Splice modes only. |
| `--alt FILE` | `mm_idx_alt_read` `index.c:648` | list of ALT contig names; affects primary selection (`mm_mark_alt`, `hit.c:90`). |

## Outputs produced

### PAF (default) — `mm_write_paf4` (`format.c:425-458`)

| Col | Content | Source |
|---|---|---|
| 1 | query name (`/1`, `/2` suffix in frag mode) | `format.c:428-430` |
| 2-4 | query length, start, end (0-based, end exclusive) | `format.c:436` |
| 5 | strand `+`/`-` | `format.c:436` |
| 6-9 | target name (or rid if names dropped), length, start, end | `format.c:437-443` |
| 10 | matches: `mlen` = sum of anchor spans before alignment (`hit.c:8-21`), exact matches from CIGAR after (`align.c:254`) | `format.c:444` |
| 11 | block length `blen` | `format.c:444` |
| 12 | MAPQ 0-60 (`hit.c:432-485`) | `format.c:445` |

Tags (`write_tags` `format.c:397-423`):

| Tag | Meaning | Present |
|---|---|---|
| `tp:A:P/S/I/i` | primary / secondary / primary inversion / secondary inversion | always |
| `cm:i` | number of anchors (minimizers) in the chain | always |
| `s1:i` | chaining score | always |
| `s2:i` | best secondary's chaining score (`subsc`) | primaries only |
| `NM:i`, `ms:i`, `AS:i`, `nn:i` | mismatches+gaps; peak segment DP score (`dp_max0`, what `-s` tests); total DP score; ambiguous bases | with `-c`/`-a` |
| `ts:A:+/-` | transcript strand | splice |
| `de:f` | gap-compressed divergence `1 - mlen / (blen + n_ambi - n_gap + n_gapo)` (`align.c:997-1003`) | with `-c`/`-a` |
| `dv:f` | minimizer-based divergence `1 - (n_match/n_tot)^(1/k)` (`esterr.c:62`) | without CIGAR, non-sr |
| `zd:i` | chain was split by Z-drop (1 = left piece, 2 = right piece, `hit.c:122`) | when split |
| `rl:i` | query bp covered by seeds dropped as too frequent (`seed.c:117-122`) | always (`format.c:447`) |
| `cg:Z` | CIGAR (`-c`) | `format.c:448-453` |
| `cs:Z`, `ds:Z`, `MD:Z` | difference strings (`--cs[=long]`, `--ds`, `--MD`) | `format.c:454-455` |

Unmapped PAF line only with `--paf-no-hit` (`map.c:620-624`, `format.c:431-434`).
Secondaries are printed unless `--secondary=no` (`map.c:612-613`). With `-P`
every hit is `tp:A:P`.

### SAM (`-a`) — `mm_write_sam3` (`format.c:522-679`)

Standard SAM; secondary hits get flag `0x100` (`format.c:556`) and by default
no SEQ/QUAL (`format.c:587,628`; `--secondary-seq` changes that); supplementary
`0x800` with hard clips unless `-Y`. Header from `mm_write_sam_hdr`
(`format.c:118`), omitted `@SQ` for multi-part indexes without `--split-prefix`
(`main.c:451-455`). Tags as in PAF plus `s1/s2/zd/rl` and optional `--MD`,
`--cs`, `--eqx` (`=/X` operators, `align.c:183`), `-L` long CIGAR in `CG` tag.

### `--write-junc` BED-like junction output (`map.c:601-607`, `format.c:263`)

Splice modes only; primaries with MAPQ >= 10.

## On-disk index (`.mmi`) — `mm_idx_dump` / `mm_idx_load` (`index.c:475-569`)

```
"MMI\2"                                  magic (index.c:481)
uint32 w, k, b, n_seq, flag              (index.c:480-482)
for each seq: uint8 name_len, name bytes, uint32 len    (index.c:483-494)
for each of 1<<b buckets:
    int32 n_p; uint64 p[n_p]             positions of multi-occurrence minimizers
    uint32 hash_size; (uint64 key, uint64 val) x hash_size
if !(flag & MM_I_NO_SEQ): uint32 S[(sum_len+7)/8]   4-bit packed bases
```

A multi-part index is simply several of these concatenated (`-I`), one per
`mm_idx_reader_read`. The file is architecture-dependent (host endianness,
raw structs). Minimizer key stored = `hash >> b << 1 | is_singleton`; a
singleton's value is its packed position, otherwise `start_in_p << 32 | count`
(`index.c:256-266`).

## `--split-prefix` temporary files (`splitidx.c`, `map.c:589-600`)

`PREFIX.NNNN.tmp` per index part: `k`, `n_seq`, names/lengths, then for every
query `n_reg, rep_len, frag_gap` and raw `mm_reg1_t` (+ `mm_extra_t` with CIGAR)
records (`map.c:590-600`), merged by `merge_hits` (`map.c:476-539`).

## Debug output to stderr

| Flag | Lines | Code |
|---|---|---|
| `--print-seeds` (forces `-t 1`) | `RS rep_len`; `SD ref pos strand qpos span gap`; `CN chain_id ...` | `main.c:220`, `map.c:255-260, 326-330` |
| `--print-chains` | `CN` lines only | `main.c:255`, `map.c:326` |
| `--dbg-seed-occ` | `SF qpos occ filtered` per query minimizer found in the index | `main.c:264`, `seed.c:115-116` |
| `--print-qname` | `QR name tid len`, `QT name tid seconds`, `QM arena stats` | `map.c:433-436, 472-473, 368-369` |
| `--print-aln-seq` (forces `-t 1`) | DP parameters, target and query bases, score and CIGAR for every ksw2 call | `main.c:224`, `align.c:339-346, 361-367` |
| `-v3` (default) | `mid_occ = N`, index stats (distinct minimizers, singleton %, mean occurrence, spacing), batch timings, peak RSS | `options.c:81-82`, `index.c:117,132-133`, `main.c:523` |
