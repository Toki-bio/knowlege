---
title: "ssearch36-fasta36: Library handling: huge sequences, scores per entry, multiple HSPs"
type: reference
generated: true
---

> Part of the [ssearch36-fasta36](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/ssearch36-fasta36/references/library-chunking-and-hsps.md)

{% raw %}
# Library handling: huge sequences, scores per entry, multiple HSPs

## 1. A 200 Mbp chromosome is split into 150,000-residue chunks

| Step | Code |
|---|---|
| max chunk = `MAXLIB` = 150,000 residues (non-`SMALLMEM` build) | `src/defs.h:91-93`; `m_msp->max_tot = MAXLIB` `src/initfa.c:603` |
| `-N n` lowers the chunk size (`ldb_info.maxn`), never raises it above 150,000 | `src/doinit.c:413-415`; `reset_maxn()` `src/compacc2.c:155-158` |
| overlap fixed at 150 residues; continuation chunks read `maxt3 = maxn-150` new residues | `src/comp_lib9.c:893-894`; `compacc2.c:172-173`; `comp_lib9.c:2505-2508` |
| reader sets `lcont` when its buffer is full before the next `>` | `agetlib()` `src/nmgetlib.c:574-705` (flag at 699-701) |
| last 150 residues copied into `aa1save`, prepended to the next chunk | `comp_lib9.c:2644-2646`, `2505-2508` |
| chunk bookkeeping: `cont` = chunk number (1-based), `l_offset` = residues before this chunk, `libstr` kept | `comp_lib9.c:2546-2550`, `2655` |
| continuation chunks do not increment `db.entries`; `n1` includes the overlap | `comp_lib9.c:2616-2620` |
| `n1tot_p` -> full original length, printed as `(N nt)` in hit lines | `comp_lib9.c:2596-2611`; `mshowbest.c:467` |
| `-M lo-hi` filters on chunk length `n1`, not sequence length | `compacc2.c:2990-2996` |

So one chromosome of 200 Mbp becomes 1,334 entries of 150,000 (the last
shorter), each scored independently on both strands. The man page's `-N`
description (`doc/fasta36.1:296-302`) and the guide FAQ about duplicated hits
across overlapping pieces (`doc/fasta_guide.tex:2126-2136`) describe this.
The reverse strand is implemented by reverse-complementing the *query*
(`comp_lib9.c:1133-1141`; `qframe = 2`, `initfa.c:1583`), so each chunk is
scanned twice; `-3` = forward only, `-i` = reverse only (`initfa.c:1043-1052`,
`doinit.c:397`).

Consequences:
- A SINE copy that straddles a chunk boundary is intact in one of the two
  chunks only if it lies within the 150-residue overlap; a 300 bp element
  centred on the boundary is cut in both (each half aligns separately).
  `notes`: `-N` cannot enlarge the overlap (`SEQDUP` 150 is unused, `defs.h:100`).
- Hit lines carry the chunk index as a `:k` suffix on the name
  (`mshowbest.c:493-496`, `mshowalign2.c:418-421`) and `_r` for the reverse
  strand (`mshowbest.c:488-492`).
- Coordinates are still whole-chromosome coordinates:
  `calc_coord(n0, n1, q_offset..., loffset + (l_off-1) + sq1off-1, aln)`
  (`mshowbest.c:633-636`, `mshowalign2.c:484-486`; `calc_coord()` in
  `src/c_dispn.c:529`).

## 2. The scan keeps exactly one score per chunk per strand

- `do_work()` returns the single maximal SW score of the chunk
  (`src/dropgsw2.c:829`). There is no list of local maxima in the kernel.
- `save_best()` makes one `beststr` per (chunk, strand) result whose z-score is
  above `zbestcut` (`src/compacc2.c:2558-2611`); initially everything is kept
  (`zbestcut = -FLT_MAX`, `comp_lib9.c:917`), and `MAX_BEST` 60,000 entries are
  retained (`defs.h:104`, pruning `compacc2.c:2578-2585`).
- Therefore, at the *scan* level, ssearch36 reports at most 2 hits (one per
  strand) per 150 kb chunk, i.e. the question "how many SINE copies in a 200 Mbp
  chromosome" is answered by ~1,334 x 2 best scores, nothing more.

## 3. Multiple HSPs per chunk appear in the alignment phase

FASTA36's headline feature (`doc/README.versions:7-19`) is implemented in
`nsw_malign()` (`src/wm_align.c:133-297`):

1. Align the chunk once with `sw_walign()` (`185-187`); rescore with the
   `-S` matrix (`196-200`).
2. If `do_rep` is off or the score `<= score_thresh`, stop (`202`).
3. Otherwise recurse on the library prefix `[0, min1)` (`204-233`) and suffix
   `[max1, n1)` (`235-276`) if each is at least `min(MIN_LOCAL_LEN=33, n0)` long
   (`172`; `defs.h:89`). Sub-alignments that beat `score_thresh` are merged into a
   score-sorted chain (`293-294`, `merge_ares_chains`).
4. Recursion continues inside each sub-piece, so the result is *every*
   non-overlapping (in library coordinates) local alignment with score
   `> score_thresh`, found greedily best-first. There is **no cap on the count**
   for ssearch36/fasta36 (`-K` is parsed only `#ifdef LALIGN`,
   `src/initfa.c:1201-1204`).
5. Query overlap is allowed: two copies of the SINE in one chunk both align to
   the whole query. Only library positions are partitioned.

`score_thresh` = `bbp->repeat_thresh`, computed just before the alignment pass
(`src/comp_lib9.c:1534-1550`; same logic `mshowbest.c:375-383`):
- if the chunk's own E() `> e_cut_r`: `repeat_thresh = score*10` (never
  satisfied -> no extra HSPs);
- else `repeat_thresh = min(E1_to_s(e_cut_r, n0, n1, zdb_size, stats), score)`,
  i.e. the raw score that corresponds to E() = `e_cut_r` for a chunk of this length.

`e_cut_r` defaults: protein `e_cut/10` (`initfa.c:952`); **DNA 0.001**
(`initfa.c:1586`, applied when the query is DNA and `-E` was not given).
Setting it (`initfa.c:1113-1141`):

| Option | `e_cut` | `e_cut_r` | `do_rep` |
|---|---|---|---|
| none (DNA) | 2.0 | 0.001 | 1 |
| `-E e1` (single value) | e1 | `e1/e1 = 1.0` if e1 >= 1, else e1 (`1126,1137-1140`) | 1 |
| `-E "e1 e2"`, e2 >= 1 | e1 | `e1/e2` (`1138`) | 1 |
| `-E "e1 e2"`, 0 < e2 < 1 | e1 | `e2` (`1139`) | 1 |
| `-E "e1 0"` or `-E "e1 -1"` | e1 | 1e-100 | 0 (no extra HSPs, `1116-1119`) |

The man page says the default is `e_cut/10` (`fasta36.1:200-205`); for DNA the
code uses 0.001 and for a single `-E` value `>= 1` it uses 1.0.

For a SINE-vs-genome scan where every copy should be reported use e.g.
`-E "0.01 1"` (both thresholds 0.01) or `-E "1 1"`; the recursion then yields
every copy in the chunk with E() <= `e_cut_r`.

## 4. Which hits get aligned at all

- Only the top `nshow` chunks (E() < `e_cut`, `comp_lib9.c:1364-1367`) are
  aligned, and `-d N` (`ashow`, `initfa.c:1109-1112`) limits how many alignment
  *blocks* are printed; `-b N` limits hit lines. `-m 8`/`-m 9c`/`-m B` force
  the alignment pre-pass so every listed hit has its HSP chain
  (`comp_lib9.c:1524-1556`; `mshowbest.c:414-465`).
- A chunk whose best copy is below `e_cut` is never aligned, so weaker copies in
  that chunk are invisible even if they would individually pass `e_cut_r`
  (`repeat_thresh` logic above, `comp_lib9.c:1537-1539`).
- Rows: `showbest()` prints one line per `a_res` in the chain
  (`mshowbest.c:504-741`); `showalign()` prints each as a block, secondary ones
  headed `>--` (`mshowalign2.c:485-600`).

## 5. `-m` and other knobs that touch this

- `-m 8CBl`, `-m 9c`: per-HSP coordinates both sides (see `io-formats.md`).
- `-d 0`: no alignment blocks, but `-m 9c` rows still include all HSPs.
- `-a`: show full sequences around the alignment (`initfa.c:1054`).
- `-K n`, lalign36 only: maximum number of non-intersecting alignments passed as
  `max_count` to `SIM()` (`dropgsw2.c:1029`; `lsim4.c:109`); default 50
  (`doinit.c:288`).
- `-N n`: smaller chunks -> more entries, more statistics points, more boundary
  cuts; not a way to get more HSPs any more (`fasta36.1:299-302`).
{% endraw %}
