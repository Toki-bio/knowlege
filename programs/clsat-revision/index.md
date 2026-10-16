---
title: "clsat-revision — digest"
type: reference
topic: Program digests
generated: true
---

> Type: reference · Claude Code skill `clsat-revision` · [source on GitHub](https://github.com/Toki-bio/knowlege/tree/main/.claude/skills/clsat-revision)

**When to use:** Resume the CLsat satellite-DNA manuscript major revision (Mobile DNA, Reviewer 1; started 2026-10-15). Author's rules, the reviewer's objections, what was found about the Darevskia assemblies (sample 215 "nairensis" is D. unisexualis; unm/unp well phased with swapped labels), CLsat measured from raw reads, the read-first locus approach, data and server locations, running jobs. Use at the start of any session on the CLsat revision, the response to reviewers, Darevskia hybrid haplotypes, satellite validation with reads, or the Amazon molly control.

---

{% raw %}
# CLsat manuscript revision

> Type: notes (project state; unpublished results, verify before citing). Last rewritten 2026-10-15.

## Read first, every session
1. `C:\work\CLsat\review\REVISION_HEADER.md`: task, status log (newest at bottom), open questions.
2. Notebook page for the author (keep updated after every step, plain language, collapsible details):
   https://claude.ai/artifact/RRVZNfHfqaAGkm2oJPYMHM, source `C:\work\CLsat\review\revision_log\index.html`.
3. Reviewer: `C:\work\CLsat\review\Reviewer1.txt`. Manuscript: `CLsat_manuscript_03_08_26.docx` (same folder).
   Not to be edited until the data are rechecked (author).

## Author's rules
- Short, atomic replies; one question per message at most; no jargon (explain every term).
- No big task or run without approval; toy-test every new pipeline before the full run.
- Never claim something is missing ("I could not find it"); never guess.
- Document every step in the notebook page and the header the same day.

## The reviewer, in one paragraph
unm (GCA_032766585.1 "maternal") and unp (GCA_032766555.1 "paternal") were treated as parental subgenomes without
a test (hifiasm HiFi-only, ~17x, no Hi-C/trio). Reviewer saw >90% of unp nairensis-like, 1/4 of unm valentini-like,
16-35 switches in large unm contigs. Asks: (i) parental affinity per contig/window, (ii) map HiFi back (collapses,
unsupported regions, switches at array flanks), (iii) re-phase with parental k-mers, (iv) hybrid claims only where
reads prove phase. Also: narrow Darevskia focus (suggests Amazon molly, Aspidoscelis), only one assembly peer
reviewed, 1.4 vs 1.9 Gb sizes, AI-heavy text.

## What we now know (2026-10-15)
- **Sample 215 (GCA_034642065.1, "D. raddei nairensis", Amberd 2021) is a D. unisexualis genome.** 51% of its
  10-kb pieces identical to unm, 39% to unp; 5 real nairensis animals (17_78, 10_78, 1_78 Pyunik; 216, 200 Amberd
  10x) fit unp at error floor (0.60-0.86%) and 215 worse (1.33-1.94%); 215's own HiFi reads (2021 zip
  X201SC21070615) give the same profile as uni 224 reads. Lab mix-up not yet confirmed (author asking colleagues).
  -> 215 out; every "nai" result in the manuscript is void.
- **unm/unp are well phased, labels swapped.** 30M short reads of real nai 17_78 and val 5_38 on 100-kb blocks:
  unp nairensis-like 99.0%, unm valentini-like 99.3%; 2 and 6 call changes along contigs. Reviewer's switches likely
  came from using 215 as "nairensis". Parent of origin can return (unp = mother, unm = father). ~140 blocks to
  check with reads. Mixta is NOT a uni parent (reads 3.1-3.4% from both halves). Arm long reads split ~50/50
  mixta/valentini.
- **CLsat from raw reads** (150-bp pieces vs 4 family + 9 sub-family consensuses x10; >=80 bp at >=80%):
  real nai 1.23-1.28%; unp 1.26%; val reads 0.76-0.88%; unm 0.70%; uni reads 0.91-1.10%; val assembly 0.78%;
  arm assembly lost all CLsat2 (reads 0.27%); 10x dvl assembly ~0. A families-only reference missed CLsat1.2/1.4
  (22/5491 pieces): never use it.
- **Read-first prototype, val218** (one HiFi file, 967k reads ~12x): 9,248 CLsat reads, 93% inside arrays,
  608 array-end reads, 50 whole-array reads (arrays 0.15-17 kb). 406 locus ends (351 single-read), 20 loci joined
  by one read. 39/47 catalogue val loci have a read-defined end within 5 kb. Limit = coverage.
  Spreadsheet/paper say 34.8 Gb, 1.93M reads for val; KIT file has exactly half -> second run or duplicated fastq?
- Lab chromosome-level mixta/valentini paper draft (`D_mixta_genome_Animals_14_10_2026.docx`): contigs <300 kb
  removed -> loses 72% of mixta CLsat (7% val); use unfiltered/our own re-assemblies for satellites.
  dvl (JAIXNF, isolate 245) is the 10x Genes 2022 assembly, NOT chromosome-level.

## Plan (read-backed satellites)
1. Panel: val, mix (Monsoon/KIT re-assemblies, hifiasm 0.25), uni (one animal; unp = nai half, unm = val half),
   arm (one animal); nairensis from reads only.
2. Assembly-free CLsat amounts (done, above).
3. Read depth and end-crossing reads per long array (KIT readmap job).
4. Read-proven array catalogue (prototype); stage 2: measure within-array unit divergence, then chain inside reads.
5. Orthology from flanks; empty site proven only by one read spanning both flanks.
6. Hybrids locus by locus; parent-split re-assemblies as a cross-check.
7. Amazon molly control (skill `amazon-molly-paper`).
Code-level method profile of satpipe at `f77448a` (stages, thresholds, 20 reviewer-style weaknesses with file:line,
plug-in points for read-based checks): `references/satpipe-method-profile.md`. HitSV-based ideas and reviewer
objections C1–C12: `projects/hitsv-vs-own-tools.md` §3.4.
Research: `C:\work\CLsat\review\RESEARCH_LOCAL_REASSEMBLY.md` (Khost 2017, Vondrak 2020, NucFreq/Flagger/VerityMap).

## Where things are
- Raw reads: `KIT:/data/V/toki/Darevskia_reads/` (1.2 TB). Origin: IBG3 `/stor1/Lizards/RAW/` (needs IBG VPN;
  `s.kosushkin@10.219.173.28`, key `C:\Users\T\.ssh\ibg_ed25519`).
- KIT analyses: `/data/V/toki/clsat_map/` (naicheck/, discrim/, clsat_reads/, readmap/, readfirst/, kitasm/).
- Monsoon: `/scratch/sk4386/darasm/` (re-assemblies, gsize; watcher darasm_submit.sh), `/scratch/sk4386/molly_asm/`.
- Molly families: github Toki-bio/CLsat-viewer `molly/alignments/` (ViewAlign links in README).
- Viewer: https://toki-bio.github.io/CLsat-viewer/ (repo Toki-bio/CLsat-viewer). VIEWER_MANUAL numbers (707 loci,
  L191/L425) are stale; manuscript numbers 627/541 loci.

## Traps met
- Long-read divergence: use minimap2 `de:f` (each indel once); NM/blocklen inflates with large indels.
- `python -I` ignores PYTHONIOENCODING; write UTF-8 files instead of printing.
- PowerShell ssh with inline quotes breaks; pipe a script with `bash -s` from the Bash tool.
- aria2 preallocates: file size is not progress. Conda installs on Monsoon take ~10 min; don't chain jobs on a guess.
- Copying a private key between servers is blocked by the permission classifier; ask the author.
{% endraw %}

## Reference pages

- [satpipe-method-profile](satpipe-method-profile/) — CLsat_workflow (satpipe): method profile for comparison with HitSV
