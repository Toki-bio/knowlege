# Smith-Waterman implementation (ssearch36)

## Files

| File | Role |
|---|---|
| `src/dropgsw2.c` | driver: `init_work()` 132-620 builds profiles; `do_work()` 696-870 scores one library chunk; `FLOCAL_ALIGN()` 875-990 scalar fallback (no-SIMD build); `do_walign()` 1001-1040 alignment entry; `#ifdef LALIGN` branch 1025-1029 |
| `src/smith_waterman_sse2.c` | Farrar 2006 striped kernel: `smith_waterman_sse2_word()` 34-218 (16-bit), `smith_waterman_sse2_byte()` 230-427 (8-bit) |
| `src/smith_waterman_altivec.c` | Lindahl 2004 Altivec/VMX equivalents (`_word` 23, `_byte` 1063), compiled with `-DSW_ALTIVEC` |
| `src/wm_align.c` | scalar Gotoh with traceback `sw_walign()` 37-131, recursive `nsw_malign()` 149-297 |
| `src/dropgsw2.h` | `struct f_struct`; SIMD fields guarded by `#if defined(SW_ALTIVEC) || defined(SW_SSE2)` (line 18) |
| `src/simde/x86/sse2.h` | SIMD-everywhere mapping, included with `SIMDE_ENABLE_NATIVE_ALIASES` (`smith_waterman_sse2.c:27-28`) so the SSE2 source compiles on ARM NEON (`doc/README.md`, Nov 2020 entry) |

There is **no AVX2 / AVX-512 kernel**: no `avx`/`__AVX` symbol exists outside
`src/simde/`, and the only SIMD defines are `SW_SSE2` and `SW_ALTIVEC`
(`make/Makefile.fcom:211-238`, `make/Makefile.sse_alt`). Register width is 128
bits: 16 lanes x 8-bit or 8 lanes x 16-bit.

## How the striped profile is laid out (Farrar)

- The *query* is striped, the *library chunk* is streamed. For 8-bit the query
  is divided into 16 segments of length `col_len = (n0+15)/16`; register lane i
  holds query positions i*col_len + j (`dropgsw2.c:458-475` comment,
  profile fill `551-612`). For 16-bit, 8 segments of `(n0+7)/8` (`551-552`).
- Profile per alphabet symbol: `word_score[(symbol)*iter*8 ...]` and
  `byte_score[...]`, so the inner loop does one aligned load per segment
  (`smith_waterman_sse2.c:108,310`). Segment padding beyond `n0` is score 0
  (`dropgsw2.c:558-562,578-583`).
- Bias: 8-bit lanes are unsigned, so `bias = min over pam2[i][j]`
  (`dropgsw2.c:535-547`) is subtracted from every profile value; the kernel adds
  profile then `_mm_subs_epu8(H, v_bias)` (`smith_waterman_sse2.c:331`).
  `f_str->bias = -bias` (`dropgsw2.c:605`). A profile value `> 255` after
  biasing is fatal (`596-601`). For DNA `+5/-4` the range is 9, far below the limit.
- Kernel per library residue (`smith_waterman_sse2.c:105-205`): load previous
  H column shifted by one lane (`116-122`), loop over `iter` segments computing
  `H = max(H+score, E, F)`, `E = max(E-ext, H-open)`, `F = max(F-ext, H-open)`
  (`124-160`), then the "lazy F" fix-up loop that re-propagates F across
  segment boundaries only while it can still improve H (`160-205`). Max score
  is accumulated in `v_maxscore` and reduced at the end (`207-214`).
- 16-bit kernel biases scores to -32768 (`81-85`) and returns `score + 32768`.

## 8-bit vs 16-bit and overflow fallback

1. `init_work()` sets `try_8bit = 1` (`dropgsw2.c:612`).
2. `do_work()` runs the byte kernel first (`785-795`). The byte kernel returns
   255 whenever `score + bias >= 255` (`smith_waterman_sse2.c:421-424`).
3. If `score >= 255`, `do_work()` reruns the word kernel (`801-813`) and counts
   `done_16bit++`; after `done_8bit > 500`, if more than one third of chunks
   needed a rerun, `try_8bit = 0` and all later chunks use 16-bit directly
   (`820-826`, comment says 8-bit is "roughly 50% faster").
4. Net effect for SINE vs genome: random 150 kb chunks score far below 255 and
   stay 8-bit; true hits (a 300 bp SINE at 80% identity scores ~900 under
   +5/-4) overflow and are redone, costing one extra 16-bit pass for that chunk
   only. `notes`: the 16-bit kernel saturates at 32767, which a 40 kb query
   could not exceed with +5 matches, so there is no third level.

## What the kernel does and does not produce

- Output is a single integer: the best local score in the chunk
  (`rst->score[0]`, `dropgsw2.c:829`). No coordinates, no second-best score.
- Coordinates and the alignment come from `do_walign()` -> `nsw_malign()` ->
  `sw_walign()` (`wm_align.c:37-131`), a scalar O(n0 x n1) Gotoh with an `n0`
  row buffer and full traceback via `nw_align()` (`wm_align.c:326`). This runs
  once per displayed hit and once more per recursive sub-alignment; for a 150 kb
  chunk x 300 bp query that is ~45 M cells per call.

## Scoring parameters for DNA and the flags that set them

| Flag | Code | Default for `-n` / DNA | Notes |
|---|---|---|---|
| `-n` | `initfa.c:1227-1232` | forces DNA alphabet; otherwise composition is guessed (`fasta36.1:106-110`) | required for stdin queries |
| `-r "+m/-n"` | `initfa.c:1248-1258`, matrix built by `mk_n_pam()` `1731` | `+5/-4` (`initfa.c:638-639`; name set `1733-1737`) | values between max/min rescaled, -1 stays -1 (`fasta36.1:323-327`); `-r` also forces DNA |
| `-s file` | `initfa.c:1265-1286` | BL50 for protein; for DNA an explicit BLAST-format matrix file can be given (`fasta36.1:343-345`) | `-s ?...` variable matrices are protein-only in practice |
| `-f N` | `initfa.c:1143-1146` (sign forced negative) | `-12` (`initfa.c:1680-1684`, `OLD_FASTA_GAP` is never defined in `make/`) | gap open |
| `-g N` | `initfa.c:1148-1151` | `-4` (`initfa.c:1686`) | gap extend per residue |
| displayed | `get_param()` `dropgsw2.c:655-691` | prints `+5/-4 matrix (5:-4), open/ext: -12/-4` | the header table in `initfa.c:82` lists ssearch DNA as "-16 -4" in the old first-residue convention |

Gap cost model: a gap of length k costs `-(gdelval + k*ggapval)` = 12+4k. The
Farrar kernels take `gap_open` = cost of the first gap residue, so `do_work()`
passes `-(gdelval+ggapval)` = 16 and `-ggapval` = 4 (`dropgsw2.c:795-800`);
`nsw_malign()` / `sw_walign()` take `-gdelval`, `-ggapval` and add them
internally (`dropgsw2.c:1021`, `wm_align.c:299`).

Program defaults table (`initfa.c:150-200`): ssearch36 `BL50`, protein gaps come
from the matrix (-10/-2), `e_cut` 10.0, no ktup; the DNA override happens only
in `resetp()` when the query is DNA.

## `-S` (low-complexity) interaction

`pam2[1]` treats lower-case residues as `X`; the scan uses `pam2p[ip]` with
`ip = ext_sq_set` (`dropgsw2.c:154-160,245-250`), alignments are re-scored with
the plain matrix (`wm_align.c:189-200`). Requires a lower-case-masked library.

## Memory per thread

`workspace` 3*16*(40,000+150,000+32)+256 = ~9.1 MB (`dropgsw2.c:456`), plus
profiles `(n0+32)*(nsq+1)*3` bytes, allocated once per strand (`f_str[0]`,
`f_str[1]`) per worker (`work_thr2.c:313-315`; `comp_lib9.c:1115-1141` for the
main thread).
