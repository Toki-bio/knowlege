---
title: "sinederella: Pitfalls and gotchas"
type: reference
generated: true
---

> Part of the [sinederella](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/sinederella/references/pitfalls.md)

{% raw %}
# Pitfalls and gotchas

- **`CHUNK_BP` / `FLANK` are misnamed.** The orchestrator documents them as
  "Step1 chunk size (30000)" and "Step1 flank size (50)" (`SINEderella:51-52,112-113`)
  and passes them as positionals 3 and 4 to step 1 (`:837`), where they are
  `SAMPLE_SIZE` (copies sampled for SubFam) and `BINSIZE` (SubFam chunk size)
  (`step1_search_extract.sh:7-8,253,261`). The real sear flank is hard-coded
  `50` at `step1:111`. Setting `FLANK=100` changes the SubFam chunk size.
- **Locus merge ignores strand.** `bedtools merge -c 4,5,6 -o max,max,distinct`
  without `-s` (`step1:218-219`): overlapping + and - hits make one locus with
  strand `+,-`, extracted as `+` (`step8a:57`); the per-consensus BEDs are
  strand-aware (`sear:537`).
- **Bitscores in `assignment_full.tsv` are sums over 10 cycles**; the 0.45
  threshold applies to the sum. The report once divided by 100 instead of 10
  (audit D17, `docs/AUDIT_2026-10-05.md:233`).
- **A missing subfamily is absorbed by its nearest relative.** With AluY
  absent, 1,305 AluY copies were firmly assigned to AluSx; with AluJb excluded,
  ~2,000 AluJ copies went to AluSx unanimously and above threshold
  (`docs/RELEASE_CHECK_2026-10-05.md:33-37`). A firm label means "closest bank
  member", not "member".
- **Fragments are not found.** 0.8 x consensus length is required of the best
  merged hit (`sear:418`); 97-99 % of RepeatMasker Alu copies missed on chr21
  are < 240 bp (`RELEASE_CHECK:30-32`). Lower the `0.8` in `step1:111` if
  fragments matter.
- **`-z 11` noise:** counts move ~1 % between identical runs; a locus on both
  the `MIN_ALN` and `MAX_E` cut-offs of the kind-A check is in or out by
  chance (`docs/SATELLITES.md:310-313`).
- **ssearch36 file-name limit (~120 characters) and id truncation (~60).**
  Every satellite call silently failed under a long run path on Sicista
  (`docs/SATELLITES.md:327-332`); fixed by running in the temp dir with bare
  names (`satellite_trf_verify.py:222-225`, `satellite_kindB_verify.py:160`).
  Other callers "have the same weakness".
- **`-z 11` with statistics keeps searching alternative alignments** of long
  pairs (56 s vs 5.5 s with `-z -1`, `satellite_kindB_verify.py:131-136`); use
  `-z -1` when only the best alignment matters.
- **One huge library loses hits** when a family dominates (622k windows: 91 %
  of VES core hits lost; fixed with 20k chunks + `-Z 20000`,
  `flankscan/HANDOFF.md:150-153`). Keep library parts small.
- **Default statistics return nothing on all-homolog libraries**; `-z 11` is
  mandatory there (`HANDOFF.md:48`).
- **Lowercase masking does nothing**; ssearch36 uppercases DNA libraries and
  `-S` is unfinished for DNA (`docs/MASKING.md:36-42`). Mask with N via
  `tools/mask_genome.py`.
- **`--add`/`--exclude` refuse masks** (`SINEderella:1347-1349`) because the
  re-vote needs the overlap between old and new families.
- **`--resume` and `results/`:** older code `rm -rf results` and destroyed
  publish alignments (audit D25/D26, fixed: `build_results_links` `:434-464`).
- **Splits are directory-global**: `sear -k` reuses any `*.2k.part_*.bnk` in
  `searches/`; `step1:127` and `SINEderella:1013` delete them after the loop.
- **Step 1 MAFFT and SubFam ignore `THREADS`** (`step1:300`, audit D8).
- **The flag table in `MANUAL.md:1218-1230`** glosses `-g` as "gapped" and
  `-3` as "both strands"; in FASTA36 `-g` takes the next token as the
  gap-extension penalty. Unverified here; harmless for results.
- **`SINEderella_multi`, flankscan on real data, `--mask-*`** were not covered
  by the release check (`RELEASE_CHECK:137-139`).
- **Executable bits** were missing from every file in a fresh clone until
  audit D11 (`docs/AUDIT_2026-10-05.md:204`); check `ls -l` after cloning.

## Defaults that surprise

| Parameter | Default | Where | Why it matters |
|---|---|---|---|
| sear Slf / Homology / flank | 0.8 / 65 / 50 | `step1:111` | not configurable by env; edit the line |
| `SPLIT_SIZE` (vote library part) | 20,000 copies | `step2:10` | library size shifts unanimity by ~1-2 % |
| `MIN_REL_BITSCORE` | 0.45 of the 10th-best summed score | `step2:11,308` | per consensus; per family it would reject 41,014 of 66,567 rsi copies (`FAMILY_SUBFAMILY:56-58`) |
| LEAK | runner >= 0.90 | `step3:300` | |
| `SATELLITE_EXCLUDE_B` | `verified` | `step1:134,142` | docs once said `flagged` (audit D2) |
| `--max-verify-runs` | 2000 | `satellite_stage.py:89` | hit-dense families leave runs `untested` |
| step 7 | base 50, step 50, cap 1000, cutoff 45 %, +5 points over background | `step7:8,120-124` | `undetermined` is not used for plates |
| step 1c | P 25, bin 20, sample 30000, min length 0.8 | `step1c:7-8` | |
| canonicalize | min-id 80, same-orientation merge 98 %, length ratio 0.9 | `SINEderella:353`, `README.md:52` | |
| `CANON_MIN_ID`, `SKIP_*` env | see `PUBLISH_WORKFLOW.md` | audit D14 | not in the usage text |

## Known bugs / open items (as of pinned commit)

- Satellite tools stream the genome once per consensus (audit D10).
- flankscan fs5-fs9 open candidates; fs9 wired only as the twin step
  (`HANDOFF.md`, audit perspectives).
- Two-layer assignment, composites' place in the hierarchy, subfamily
  criterion: evaluated, undecided (`docs/FAMILY_SUBFAMILY_ASSIGNMENT.md:120-135`).
- `benchmark_sear.sh` never reported; `sear_multi` unused.

## Added 2026-10-12 (benchmark in knowlege `projects/sine-search-strategy.md`)

- The hit loss "when the library is huge and one family dominates"
  (`flankscan/HANDOFF.md:150-153`) is FASTA36's `MAX_BEST` = 60,000 cap on
  retained scan hits (`defs.h:104`; weakest quarter dropped when full,
  `compacc2.c:2592-2602`). One ssearch36 call over 370k windows of a 1.67 Gbp
  genome kept 43 % of dar_squam1 copies, almost all above 90 % identity.
  `sear`'s 20 parts (~42k windows each for a 2 Gbp genome, ~5k hits per part
  for a 10^5-copy family) stay under the cap; never collapse them into one run.
- `fasta36 ... 4` (ktup 4) in place of `ssearch36 ... 3` at `sear:396` gave
  99.8 % of the ssearch36 copies on D. valentini in 1/8 of the time, same
  `-m 8C` columns.

## Added 2026-10-13 (KIT install + Timema rerun at `38d811a`)

- **The satellite stage is gated on `trf` being on PATH** (`step1_search_extract.sh:138`); without it the log says
  "Satellite screen skipped" and the run goes on. `manifest.txt` can say `trf not found` even when the stage ran
  (KIT: `trf` lives in `/data/V/toki/bin`). Two runs on the same inputs are not comparable unless both did or both
  skipped the screen. Pin it explicitly: `SKIP_SATELLITES=1`, or put `trf` on PATH.
- **Kind-B unit check cost on a hit-rich genome:** Timema (1.23 Gb, 14 consensuses, THREADS 24): satellite stage
  1 h 35 min, of which ~85 min is kind B (2,000 runs per big family, units median 5 kb, max 24 kb, full SW) to find
  2–12 arrays per family; whole screen removed 1,310 hits (~0.1 %). Searches took 28 min. `SATELLITE_EXCLUDE_B=none`
  keeps kind A (~8 min) and skips the unit check. Owner's verdict: long alignments are against the design.
  Replacement spec: skill `tandem-array-unit-check`.
- **Do not pull or edit a checkout a run executes from.** bash reads scripts incrementally; a run started from
  `repos/SINEderella` breaks if another session updates that checkout mid-run. Develop in a separate clone.
{% endraw %}
