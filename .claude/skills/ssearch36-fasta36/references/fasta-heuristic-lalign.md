# fasta36 heuristic vs ssearch36; lalign36

## fasta36 for DNA (`src/dropnfa.c`)

Pipeline per library chunk, `do_fasta()` (`dropnfa.c:567-832`):

1. **ktup word lookup.** The query is hashed into all `ktup`-mers at
   `init_work()` (`dropnfa.c:202-300`); the library chunk is scanned, each
   ktup-mer looks up query positions (`671-745`) and runs are accumulated on
   diagonals (`diag[]`, scored with `pamh1/pamh2`). Low-complexity (`NMAP`)
   residues break words (`672-690`).
2. **init1 / spam.** The best diagonal regions (`MAXSAV` of them, `savemax`
   896) are rescored with the real matrix by `spam()` (`928-1024`); the
   x-drop extension of a region is **disabled for DNA** unless either sequence is
   < 50 residues (`do_extend = !ppst->nt_align || n0<50 || n1<50`, `607`).
   `init1` = best region (`766-769`).
3. **initn / joining.** If `init1 >= c_gap` the regions are joined with gap
   penalties by `sconn()` (`780-786`, `1026-1090`); `initn` = joined score.
4. **opt.** If `initn > opt_cut` and `optflag`, a banded SW of width `optwid`
   around the best diagonal (`dmatch()` 1131, called `812-819`) gives `opt`,
   the score used for statistics (`score_ix` points to opt).
5. Thresholds `c_gap`/`opt_cut` are derived from E-values via `ELK_to_s()` with
   the ungapped lambda/K (`631-634`) when `use_E_thresholds` (default), scaled
   by `ktup_sq` (= ktup for ktup > 3, `610-616`).

Defaults for DNA set in `resetp()` (`src/initfa.c:1693-1716`):

| Parameter | DNA default | Code | Flag |
|---|---|---|---|
| ktup | 6 (`mktup=6`) | `initfa.c:1696`; positional 3rd argument `argv[optind+3]` `initfa.c:1469-1476` | `fasta36 -n q.fa lib.fa 4` |
| optwid (band) | 16 | `initfa.c:1699` | `-Xy N` (`1431-1436`) |
| E_band_opt (E-opt) | 0.05 | `initfa.c:1708` | `-c "E-opt E-join"` (`1082-1108`) |
| E_join | 5 x E_band_opt = 0.25 | `initfa.c:2904-2906` | second value of `-c` |
| classic optcut thresholds | off | `-c O` or `-c N` with N > 1 (`1085-1092`) | |
| opt calculation | on | `-Xo` turns off (`1420-1425`) | |
| alignment display | banded SW (`sw_flag=0`, label "banded Smith-Waterman") | `initfa.c:1711-1715`; `bd_walign()` `dropnfa.c:1559-1600` | `-A` forces full SW for the displayed alignment (`1055-1058`) |
| gaps / matrix | -12/-4, +5/-4 | same `resetp()` block as ssearch | `-f -g -r` |
| e_cut | 2.0 | `initfa.c:1584-1586` | `-E` |

The man page says DNA thresholds are E-join 0.1 / E-opt 0.02 (`doc/fasta36.1:160-161`)
and the header table in `initfa.c:80` says 0.01; the code uses 0.05 / 0.25.

Multiple HSPs: `bd_malign()` (`dropnfa.c:1325-1500`) recurses left/right like
`nsw_malign()`, re-running `do_fasta()` on each sub-piece, so a sub-alignment is
only found if the heuristic re-finds it in the sub-piece.

## Does fasta36 with small ktup approach ssearch36 for diverged SINEs?

Verified facts:
- Sensitivity is bounded by step 1: a copy must contain an exact `ktup`-mer on
  a diagonal that accumulates enough score to enter `MAXSAV` (`savemax`), and
  for DNA regions are not extended (`607`), so the region score is the sum of
  matches/mismatches within the exact-word runs.
- The guide states ktup 3-4 "provides higher sensitivity" for DNA
  (`doc/fasta_guide.tex:125-128`) and that ktup=1 may find what ktup=2 misses for
  proteins (`2071-2082`). Smallest allowed ktup is 1 (`query_parm`,
  `initfa.c:1791`), hash size grows as `nsq^ktup` (`dropnfa.c:219-223`).
- Cost: the diagonal pass is O(n1 x hits); with ktup 1-2 on DNA nearly every
  position hits (4^2 = 16 words) so the heuristic degenerates toward O(n0 x n1)
  with poorer constants than the SIMD kernel.

`notes` (interpretation, untested here): for 10-30% diverged 100-500 bp SINEs,
ktup 3 or 4 is the realistic sweet spot; ktup 1-2 buys little over 3 and costs a
lot because of hash saturation. Expect fasta36 ktup 3-4 to recover most copies
that ssearch36 finds with E() < 1e-3 but to lose the most diverged ones, and its
E-values are less reliable because only chunks passing `opt_cut` get a real
`opt` score (`valid_stat`, `dropnfa.c:639-641,820-826`). If ssearch36 time is the
problem, the cheaper wins are `-T`, 8-bit lanes (default), restricting
alignment work with `-E`/`-b`, and running strands separately with `-3`/`-i`
only if you need to split jobs.

## lalign36 (`-DLALIGN`)

- Same `dropgsw2.c` driver and SIMD scan (`make/Makefile.fcom:217-223`), but
  `do_walign()` calls `SIM()` (`dropgsw2.c:1025-1029`), the Huang-Miller
  implementation of Waterman-Eggert non-intersecting alignments
  (`src/lsim4.c:28-31`): `big_pass()` fills the matrix once (`lsim4.c:101`), then
  up to `max_count` = `-K` (default 50, `doinit.c:288`) alignments are extracted
  best-first, each followed by a local recomputation (`locate()`/`small_pass()`,
  `lsim4.c:109-180`). Alignments must not intersect in *either* sequence
  (unlike `nsw_malign`, which partitions only the library).
- `mini_score` = `repeat_thresh` from `e_cut_r` (`mshowalign2.c:260-267`);
  for lalign `e_cut_r = e_cut` (`initfa.c:954,1600-1602`), default `e_cut` 1.0
  (`initfa.c:197`), 0.2 for DNA (`1585`).
- Self-comparison: `same_seq()` sets `nseq=1` so the identity diagonal is
  excluded unless `-J` (`dropgsw2.c:707-712,1026`; `initfa.c:1190`).
- Library = usually one sequence; statistics come from shuffles (`fasta36.1:49-51`).
  Non-threaded (`Makefile.pcom_t:50`). `-m 11` lav output and `lav2ps` plots
  (`fasta36.1:273-276`).
- `notes`: lalign36 is the tool for "all local alignments between two
  sequences" (e.g. SINE consensus vs one BAC), but it is a pairwise program and
  pays `big_pass()` over the whole pair; for genome scale, ssearch36 with a low
  `e_cut_r` gives the same library-side partitioning cheaper.

## ggsearch36 / glsearch36

Global (`dropnnw2.c -DGLOBAL_GLOBAL`) and global-query/local-library
(`dropnnw2.c`), SIMD kernels `global_sse2.c`/`glocal_sse2.c`, `-z 0` statistics
by default (`initfa.c:936-938`). glsearch36 requires the *whole query* to align,
so it rejects truncated SINE copies; not suitable for the SINE-vs-genome task
except to find full-length copies.
