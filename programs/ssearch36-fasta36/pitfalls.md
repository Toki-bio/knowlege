---
title: "ssearch36-fasta36: Pitfalls and gotchas"
type: reference
generated: true
---

> Part of the [ssearch36-fasta36](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/ssearch36-fasta36/references/pitfalls.md)

{% raw %}
# Pitfalls and gotchas

- **One score per 150 kb chunk per strand in the scan.** A chromosome with 500
  SINE copies yields at most 2 x 1,334 hits from the scan; extra copies appear
  only through the recursive sub-alignment of chunks that passed `-E`
  (`src/compacc2.c:2558-2611`; `src/wm_align.c:202-296`). Copies in chunks whose
  best copy is above `e_cut` are never reported (`comp_lib9.c:1537-1539`).
- **DNA secondary-alignment threshold is 0.001, not e_cut/10.** `initfa.c:1586`
  overrides the man page (`doc/fasta36.1:200-205`). A single `-E 2` makes
  `e_cut_r = 1.0` (`initfa.c:1126,1138`); `-E 0.5` makes it 0.5. Use
  `-E "e1 e2"` explicitly.
- **`-M lo-hi` filters chunk length**, not sequence length (`compacc2.c:2990-2996`),
  so `-M 100-500` on a genome drops almost every chunk.
- **Chunk boundaries cut elements.** Overlap is fixed at 150 residues
  (`comp_lib9.c:894`); `-N` only shrinks chunks and cannot widen the overlap.
  Duplicates across the overlap show up as two hits with a shared coordinate
  (`doc/fasta_guide.tex:2126-2136`). Subject ids get `:k` suffixes
  (`mshowbest.c:493-496`).
- **`-K` is ignored by ssearch36/fasta36**; it is parsed only `#ifdef LALIGN`
  (`initfa.c:1201-1204`). There is no cap on HSPs per chunk in ssearch36.
- **ktup is positional** (`ssearch36`/`fasta36 [opts] query lib ktup`,
  `initfa.c:1469-1476`, `fasta36.1:235`); `-k` is the number of shuffles
  (`initfa.c:1206-1209`), not ktup.
- **Statistics silently change method.** `< 500` kept hits -> library shuffles
  (`comp_lib9.c:1079-1096`, doubles runtime, header shows `(shuffled [N])`);
  `< 10` scores -> Altschul-Gish (`scaleswn.c:289`); regression fit failure ->
  MLE (`scaleswn.c:468-470`). Read the `Expectation_n fit` header line.
- **E-value database size for DNA is residues / chunk length** unless `-Z`
  (`scaleswn.c:2694-2699`); `-Z` was ignored for DNA:DNA before July 2025
  (`doc/readme.v36:34-38`). Older binaries differ.
- **Query longer than 40,000 is truncated with a warning only**
  (`comp_lib9.c:710-715`, `defs.h:92`).
- **8-bit lane saturation is handled, but matrix values are not clamped:** any
  profile value > 255 after biasing aborts (`dropgsw2.c:596-601`); huge `-r`
  values (e.g. `+100/-50`) can hit this.
- **Alignment phase is scalar and per full chunk.** Thousands of hits means
  thousands of 150 kb x n0 Gotoh runs plus recursion (`wm_align.c:37,149`);
  this, not the SIMD scan, often dominates for repeats. `-N 20000` cuts the
  per-hit alignment cost ~7x at the price of more boundary cuts.
- **`-m 8` percent identity includes gaps** since Dec 2024 (`doc/readme.v36:57-62`);
  `-XG` restores ungapped identities. `-XB` was repurposed (`readme.v36:67-71`).
- **fasta36 DNA alignments are banded** (width 16, `initfa.c:1699,1711-1715`)
  unless `-A`; a banded alignment of a diverged copy with indels can be worse
  than the score suggests.
- **`-c` thresholds in the man page are stale** for DNA (0.1/0.02 vs code
  0.25/0.05, `fasta36.1:160-161`, `initfa.c:1708,2904-2906`).
- **Reverse-strand coordinates** are on the reverse-complemented query
  (`comp_lib9.c:1133-1141`); the hit is marked `_r` / `[r]`
  (`mshowbest.c:488-492,536-539`).
- **`-S` requires a lower-case-masked library**; `-S` on an upper-case genome
  does nothing.
- **Multi-query runs keep the whole library in RAM** up to 16 GB (`defs.h:70-74`);
  `-XM 2G` to limit (`initfa.c:1398-1407`).

## Defaults that surprise

| Parameter | Default | Why it matters |
|---|---|---|
| `-E` (DNA) | 2.0 (`initfa.c:1585`) | ~2 random hits per genome scan by construction |
| `e_cut_r` (DNA) | 0.001 (`initfa.c:1586`) | weak secondary copies hidden |
| `-b` | 50 hit lines when E-values are off (`doinit.c:296`); unlimited below `-E` otherwise (`fasta36.1:147-154`) | `-b =N` to force |
| `-d` | -1 = all shown hits (`doinit.c:299`); `-m 8` sets 0 (`doinit.c:925`) | |
| `-T` | all CPUs (`doinit.c:255`) | |
| `-N` | 150,000 (`defs.h:93`) | only lowers |
| `-z` | 1 (`initfa.c:939`); 0 for gg/glsearch (`936-937`) | |
| `-k` | 500 shuffles (`initfa.c:969`) | also the "few hits" trigger |
| `-f/-g` DNA | -12/-4 (`initfa.c:1680-1686`) | header prints open/ext; table in `initfa.c:82` uses the old -16 convention |
| `-r` | +5/-4 (`initfa.c:638-639`) | |
| ktup DNA | 6 (`initfa.c:1696`) | 3-4 for sensitivity (`fasta_guide.tex:127-128`) |

## Known bugs / open upstream issues (as of pinned commit)

- `doc/readme.v36:21-26`: tfastx36 reverse-frame annotation coordinate still
  off by one (translated searches only).
- `doc/readme.v36:90-94`: `last_thresh()`/`E1_to_s()` zsflag handling was fixed
  Aug 2024 / Feb 2025; older binaries may compute `repeat_thresh` with the
  wrong method.
- `SEQDUP` 150 is defined "future - overlap" (`defs.h:100`) but the overlap is
  hard-coded in `comp_lib9.c:894`.
{% endraw %}
