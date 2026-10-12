# Input and output formats

## Inputs accepted

| Format | Reader | Notes / quirks |
|---|---|---|
| FASTA/FASTQ, optionally gzipped, for `mem`, `fastmap`, `bwasw` | `kseq.h` via `bseq_read` (`bwa.c:79-112`) | Read names ending `/1` `/2` are trimmed (`bwa.c:54-58`). Comment kept only with `-C`. Non-ACGT query bases become code 4 and break seeds (`bwt.c:296,316-319`) and score -1 (`bwa.c:142`). |
| Interleaved pairs | `-p` + `bseq_classify` (`bwa.c:114`) | adjacent identical names form a pair |
| FASTQ/BAM for `aln` | `bwaseqio.c`, `bamlite.c` | `-b -0 -1 -2` |
| Reference FASTA (gz ok) for `index` | `bns_fasta2bntseq` (`bntseq.c:280`) | `>name comment`; comment stored in `.ann` |
| `<prefix>.alt` (optional) | `bns_restore` (`bntseq.c:176-209`) | marks ALT contigs; `-j` ignores |

## On-disk index files (`<prefix>.*`; `.64.*` if `index -6`, probed first `bwa.c:253`)

Let G = forward genome length (`l_pac`), N = 2G = indexed text length.

| File | Writer | Layout | Size |
|---|---|---|---|
| `.pac` | `bns_fasta2bntseq` (`bntseq.c:314-324`) | 2 bits/base forward genome, then one pad byte if `G%4==0`, then a byte holding `G%4` | `G/4 + 1..2` bytes |
| `.ann` | `bns_dump` (`bntseq.c:70-83`) | text: `l_pac n_seqs seed` then per sequence `gi name [comment]` / `offset len n_ambs` | KBs |
| `.amb` | `bns_dump` (`84-94`) | text: `l_pac n_seqs n_holes` then `offset len char` per N-run | KBs-MBs |
| `.bwt` | `bwt_dump_bwt` (`bwt.c:385-394`) | `primary` (8 B), `L2[1..4]` (32 B), then `bwt_size` uint32 words: per 128 bases 4x uint64 Occ + 8 words of 2-bit text (`bwt.h:74-75`) | `40 + N/4 + 32*(N/128+1)` ~ `N/2 = G` bytes |
| `.sa` | `bwt_dump_sa` (`396-407`) | `primary`, `L2[1..4]`, `sa_intv`, `seq_len` (7x8 B), then `n_sa-1` uint64 | `56 + 8*(N/32)` ~ `N/4 = G/2` bytes |

Loader: `bwt_restore_bwt` derives `bwt_size` from file size (`bwt.c:451`), `bwt_restore_sa` checks `primary` and `seq_len` match (`429,433`). The `.sa` loaded array has `n_sa` entries (`sa[0]` = -1 sentinel, `437`).

Intermediate states: after `pac2bwt` the `.bwt` lacks Occ and "CANNOT be used with BWA" until `bwtupdate` (`bwtindex.c:128`).

`.sai` (`aln` output): magic `SAI\1` (`bwtaln.h:36`), the `gap_opt_t` struct, then per read `int n_aln` + `bwt_aln1_t[n_aln]` (`bwtaln.c:178-179,214-218`). Binary, struct-layout dependent.

## Outputs produced

### `bwa mem` SAM (`mem_aln2sam`, `bwamem.c:851-976`)

| Column/tag | Source | Flags that change it |
|---|---|---|
| FLAG | 0x1/0x8/0x20 from mate, 0x4 unmapped, 0x10 reverse (`858-866`), 0x100 secondary (`mem_reg2aln:1137`), 0x800 supplementary for extra primaries (`1057-1058`), 0x100 instead of 0x800 with `-M` | `-M`, `-5` |
| MAPQ | `mem_approx_mapq_se` (`982-1006`) for primaries; 0 for secondaries (`1136`); supplementary capped to primary MAPQ unless `-q/-5` (`1059-1060`) | `-q -5 -Q` |
| CIGAR | banded global DP (`bwa_gen_cigar2`), soft clip `S` on primary, hard clip `H` on supplementary unless `-Y` (`add_cigar` `838-849`) | `-Y` |
| SEQ/QUAL | `*` for secondary (`897-898`); hard-clipped bases removed for supplementary (`901-904`) | |
| `NM:i`, `MD:Z` | `bwa_gen_cigar2` (`bwa.c:196-225`) | |
| `MC:Z`, `MQ:i` | mate CIGAR / MAPQ (`934-935`) | PE |
| `AS:i` | best local score (`936`) | |
| `XS:i` | suboptimal score `max(sub, csub)` (`937`, `1185`); suppressed on secondaries (`1056`) | |
| `RG:Z` | from `-R` (`938`) | |
| `SA:Z` | `chr,pos,strand,CIGAR,mapq,NM;` for the other primary/supplementary pieces (`939-958`) | |
| `pa:f` | score / ALT score (`959-960`) | ALT only |
| `XA:Z` (or `XB:Z` with `-u`) | `chr,±pos,CIGAR,NM;` per alternative hit (`bwamem_extra.c:146-163`); XB adds `,score,mapq` | `-h -z -u -a` |
| `XR:Z` | reference FASTA comment (`967-974`) | `-V` |
| Header | `@HD` + `@SQ` (`AH:*` for ALT) + `-H`/`-R` lines + `@PG` (`bwa_print_sam_hdr`, `bwa.c:407-439`) | |

Unmapped reads get a record with `*\t0\t0\t*` (`878`, `1063-1067`).

### `bwa fastmap` (`fastmap.c:448-475`)

```
SQ  <name>  <len>  [<seq> with -p]
EM  <qbeg>  <qend>  <n_occ>  <chr:±pos>...   (positions only if n_occ <= -w, else "*")
//
```
`qbeg/qend` are 0-based half-open on the query; positions are 1-based leftmost on the forward strand of the hit, strand `+`/`-` (`466-469`).

### `bwa bwasw` SAM (`print_hits`, `bwtsw2_aux.c:491-540`)

MAPQ up to 250 (`427`); repetitive ("random") hits print MAPQ 0 and CIGAR `*` (`503`); all hits for a read are printed as separate lines, secondary flag only with `-M` (`497`).

### `bwa samse/sampe` SAM (`bwase.c:386-500`)

`XA:Z` multi-hits only when total occurrences `<= -n+1` (`bwase.c:57-58`); `X0/X1` counts of best/second-best hits.
