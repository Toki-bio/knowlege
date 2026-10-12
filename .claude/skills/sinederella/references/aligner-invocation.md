# Every aligner invocation (Q1)

Line numbers at `90d8ac4`. ssearch36 is FASTA 36.3.8g on the author's servers
(`docs/SATELLITES.md:327`). Flag semantics in the notes column are my reading
of the FASTA36 manual, not run here (ssearch36 is not installed in this
container): `-g #` gap-extension penalty (so `-g -3` sets it to -3),
`-f #` gap-open (never set -> default), `-n` DNA query, `-Q` quiet,
`-z 11` shuffle-based statistics, `-E` E-value cut, `-Z N` pretend library
size N, `-d N` max alignments shown, `-w`/`-W` display widths (`-w` also sets the subject-ID width in
`-m 8` output, 60 by default; see `projects/sear-audit-2026-10.md`), `-m 8`/`8C`/`8CB` BLAST tabular (plain / with comments / with BTOP),
`-T` threads, trailing positional = ktup (ignored by ssearch).
`MANUAL.md:1218-1230` glosses `-g` as "gapped" and `-3` as "both strands";
both strands are searched in every call here either way (minus-strand hits
appear in all outputs).

## A. Genome search (step 1)

| # | Call | Where | Query / library | Output parsing |
|---|---|---|---|---|
| A1 | `ssearch36 -g -3 -T ${THREADS:-nproc} -Q -n -z 11 -E 2 -w 95 -W 70 -m 8C "$query" "$bank" 3` | `sear:404` (called from `step1_search_extract.sh:111` and `SINEderella:1004`) | one consensus vs one of 20 part files of 2500 bp windows (step 2000), `seqkit sliding -g -s 2000 -W 2500` + `split2 -p 20` (`sear:379-380`) | awk `/^[^#]/`: `$2` window, `$9-1,$10` subject coords, `$10-$9+1` length, strand from sign of `$8-$7`, `$3` %id, `$12` bits (`sear:405-409`); `sort -V`; `bedtools merge -s -c 4,5,6,7,8 -o max,max,distinct,max,max` (`:413`); identity cut `$7>=65` -> `rawhits.tsv` (`:417`); + `$4>=int(0.8*len)` -> `plus_to_add` (`:418`) |
| A2 | `ssearch36 -g -3 -T T -Q -n -z 11 -E 2 -w 95 -W 70 -m 8C "$query_in" "$bank" 3` | `sear_multi:267-268` (standalone) | whole bank FASTA as multi-query vs each part | awk keeps `$1` qname and splits per query (`:269-293`); same merge/filter per query (`:305-308`) |
| A3 | `minimap2 -x asm10|asm20 -t T --secondary=yes -p 0.3 -N 10000 genome query` | `sear:298-301`, only query >= 1000 bp | whole genome | PAF: identity = `matches/min(qspan,tspan)`, `AS:i` tag (`:302-322`); merge then `$4>=65 && span>=0.8*len` (`:337-344`) |

Search genome: `genome.clean.fa` or `genome.masked.fa` (`SEARCH_GENOME`,
`step1:30-38`). No scoring matrix, no gap-open, no `-Z`, no `-d/-b`, no `-S`.

## B. Assignment vote (step 2) and its reproductions

| # | Call | Where | Query / library |
|---|---|---|---|
| B1 | `ssearch36 -g -3 -T "$THREADS" -Q -n -z 11 -E 2 -w 95 -W 70 -m 8 consensus_clean.fa part >> raw_<cycle>.m8` | `step2_asSINEment.sh:186-188` (full), `:634-636` (incremental Tier 2) | all consensuses (query) vs a part of <= 20,000 extracted copies (library), 10 cycles |
| B2 | same | `asSINEment:162` (old standalone), `tools/vote_eval/vote_full.sh:36` (keeps top-6 per copy with coordinates for re-tallying) | |

Parsing (`step2:193-213`): per `(copy=$2, cons=$1)` keep max `$12`; per copy
the best and the best distinct runner-up -> `copy cons bits runner runner_bits`.
`-E 2` applies to a 20,000-sequence library. Audit D18 measured that `-z 11`
bitscores do not depend on library size (340 vs 346 bits for the same
self-alignment, `docs/AUDIT_2026-10-05.md:234`).

## C. Similarity, plots, diagnostics

| # | Call | Where | Purpose / parsing |
|---|---|---|---|
| C1 | `ssearch36 -Q -n -z 11 -E 100 -T T -m 8 cons cons` | `step3_postprocess.sh:80-81` | self bits: max `$12` where `$1==$2` (`:83-88`) |
| C2 | `ssearch36 -Q -n -z 11 -E 100 -T T -m 8 cons_single.fa copies.fa` | `step3:128-130` | best bits per copy (`:138-143`) -> `sim_ratio = bits/self` (`:146-157`) |
| C3 | `ssearch36 -g -3 -Q -n -z 11 -E 100 -T T -m 8 cons copies` | `step4_plots.sh:139-141` | % identity of best hit per copy (`$3` at max `$12`) for divergence histograms |
| C4 | `ssearch36 -Q -n -z 11 -E 100 -T T -m 8 ...` and `ssearch36 -m 8 -T T all.fa all.fa` | `step4_diagnostic.py:152`, `:1061` (historical) | diagnostic positions / SINEplot matrix |

## D. Satellite stage (inside step 1)

| # | Call | Where | Purpose |
|---|---|---|---|
| D1 | `ssearch36 -m 8 -E 0.01 -z 11 -T T units2.fa cons.fa` (cwd = temp dir, bare names) | `tools/satellite_trf_verify.py:225` (`MAX_E`) | kind A: TRF unit (doubled) vs consensus, >= 45 bp aligned |
| D2 | `ssearch36 -m 8 -E 10 -z -1 -T 1 query_unit partners` per unit, `--threads` processes in parallel | `tools/satellite_kindB_verify.py:137` (`SSEARCH_OPTS`), `:181-219` | kind B: identity of lag-1/lag-2 unit pairs; `-z -1` (no statistics) because with statistics ssearch36 keeps searching alternative alignments (56 s vs 5.5 s per 26 x 20 kb pair, `:131-136`) |
| D3 | `ssearch36 -m 8 -E 10 -z 11 -T T units library` all-vs-all per run | `:104` (`verify_runs_serial`, `--serial`) | reference implementation (13 h 24 min on rsi vs 9.5 min, `docs/SATELLITES.md:282-306`) |
| D4 | `trf windows.fa 2 5 7 80 10 40 300` (approx.; `TRF_ARGS`) | `tools/satellite_trf_verify.py:177` | kind A tandem detection on +-1 kb windows |

## E. flankscan, length variants, composites (post-assignment)

| # | Call | Where | Note |
|---|---|---|---|
| E1 | `ssearch36 -m 8 -E $EVAL(1e-5) -Z 20000 -z 11 -T T cons.masked.fa chunk` over `seqkit split2 -s 20000` chunks of windows | `flankscan/fs3_partners.sh:85-96` | **the one place `-Z` is used**: one search over 622k windows lost 91 % of VES core hits; 20k chunks with fixed `-Z` recovered them (`:86-89`, `HANDOFF.md:150-153`) |
| E2 | `ssearch36 -m 8 -E 1e-5 -z 11 -Z 1000 all.fa all.fa`; `-E 1e-3 -Z 1000` fold/ends checks | `flankscan/fs_lib.sh:124,140,178` | candidate self-comparison, same-element fold |
| E3 | `ssearch36 -m 8 -E 1 -Z 1 -z 11 -T T cons.masked.fa flanks.fa` | `flankscan/fs8_singletons.sh:83` | |
| E4 | `minimap2 -x asm20 -k12 -w4 -N 10 -p 0.2 --secondary=yes -t T genome q.fa` | `flankscan/fs9_dupregion.sh:30` | segmental duplications around copies (the only other minimap2 use) |
| E5 | `ssearch36 -m 8CB -E 1e-3 -d 5 -z 11 -T 1 cons.fa copies` | `tools/length_variants.py:391` | BTOP parsed per query column (`:77`) for end modes |
| E6 | `ssearch36 -m 8 -E 1e-5 -z 11 -T T cons windows` | `tools/composite_scan.py:106` (prototype) | |
| E7 | `ssearch36 -m 8 -E 1e-5 -z 11 -Z 1000 ...` | `flankscan/tests/checks_5.sh:15`, `test_chain.sh:18` | toy checks |

`-z 11` is required on all-homolog libraries: "default statistics return
nothing" (`flankscan/HANDOFF.md:48`, `docs/COMPOSITES.md:24`).

## F. MAFFT (for completeness)

- step 1 chunk consensus alignment: `mafft --thread T --threadit 0 --localpair --maxiterate 1000 --ep 0.123 --nuc --reorder --preservecase --quiet` (`step1_search_extract.sh:300-308`).
- plates: same L-INS-i with `--adjustdirection --thread 1`, subfamilies in parallel (`step8a_extract_alignments.sh:50-59`); `--threadit 0` everywhere iterative since 2026-10-02 (`PUBLISH_WORKFLOW.md:121-123`).
- `tests/bench/run_bench.sh` compares L-INS-i 1 thread / 16 threads / `--maxiterate 2` / FFT-NS-i / FFT-NS-2 by time and Q score.

## G. blastn

Not called anywhere in the code. Measured once against ssearch36 on a planted
copy (`docs/MASKING.md:26-32`): blastn word 11 recovers 20/20 at 100 %, 4/20 at
85 %, 0/20 at <= 75 %; blastn-short 20/20 at 85 %, 17-20/20 down to 45 %;
ssearch36 20/20 throughout.
