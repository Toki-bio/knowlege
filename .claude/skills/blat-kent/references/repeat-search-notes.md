# Notes: using BLAT as a pre-filter for diverged, high-copy repeats (SINE consensus vs genome)

> Type: notes. Interpretations and recommendations built on the verified facts in the other reference
> files; probability numbers are from a simple independent-substitution model (script in the digest session,
> not from the code). Nothing here was benchmarked.

## The problem restated against the code

Query: SINE consensus 100-500 bp; copies 10-30 % diverged (70-90 % identity, with indels); copy number
10^4-10^6 per genome. Four BLAT mechanisms, each verified in `index-and-seeding.md` / `scoring-and-filters.md`
/ `pitfalls.md`, work against this:

| Stage | Default behavior | Effect on SINE search |
|---|---|---|
| Index | tiles occurring >= 1024 times dropped (`genoFind.c:857`, `:2751-2752`) | consensus 11-mers of an abundant family are all overused: near-zero seeds |
| Seeding | 2 exact 11-mers on the 11-bp target grid, diagonal within 2 (`genoFind.c:1924-1938`) | ~0.3 chance per 300 bp copy at 75 % identity, ~0.7 at 80 % (table in `index-and-seeding.md`) |
| Alignment | fuzzyFind `ffCdna`, aborts if the exact-tile skeleton is improbable (`fuzzyFind.c:1222-1229`) | tuned for near-exact cDNA |
| Output | `-minIdentity=90` (`gfOut.c:134`), 16 alignments per target sequence per strand (`gfBlatLib.c:22`, `supStitch.c:837-841`), `-maxIntron=750000` chaining (`supStitch.c:815-820`) | 75 % hits discarded; at most 16 copies per chromosome; fragments 100s of kb apart merged |

## Tuning levers, with what the code allows

1. **repMatch**: set `-repMatch` to a value above the family copy number times the grid-sampling fraction
   (`-repMatch=5000000` effectively disables dropping; memory = kept entries x 4 B, see table). Do not use
   `-ooc`. Alternatively tileSize > 12 disables filtering outright (`genoFind.c:673`) but then seeds are
   13-mers (worse sensitivity) and entries cost 6 B.
2. **tileSize / stepSize**: `-tileSize=8 -stepSize=4` gives P(>= 2 hits) ~ 0.93 at 70 % and ~1.0 at 75 % for
   300 bp (model); index = 65536 tile space, 750 M entries for 3 Gbp (3 GB). `-tileSize=9 -stepSize=9`
   (~0.72 at 75 %) costs 1.3 GB. Query cost scales with total list length touched per query position:
   with K=8 and no dropping, average list length for a random 8-mer is 3e9/4/65536 ~ 11k entries per
   position (750 M / 65536), i.e. ~4-5 M hits for a 400 bp query, sorted per 64 kb bucket. Expect queries to
   be slow (seconds+), but for a handful of consensus sequences that is fine.
3. **oneOff**: `-oneOff` raises per-tile sensitivity ~4-6x at 70-80 % (p^K + K p^(K-1)(1-p)); P(>= 2 hits) for
   300 bp with K=11 s=11 goes from 0.32 to 0.98 at 75 %. Costs 34 lookups per query position. Works only in
   standalone blat (`genoFind.c:1050-1051`). Combining `-tileSize=11 -oneOff` keeps the index at 1.1 GB.
4. **minMatch**: `-minMatch=1` makes every surviving hit a clump (`genoFind.c:1938`) -> huge clump counts
   for an abundant repeat; keep 2 and rely on smaller K / oneOff instead, unless the index already
   suppresses chance hits.
5. **minIdentity / minScore**: `-minIdentity=60` (or 0) and `-minScore=30-50`. Identity is on block bases
   with target gaps free (`gfOut.c:131-135`); score is +1/-1 minus log gap penalties, so a 300 bp copy at
   75 % with a few indels scores ~130-150.
6. **maxIntron**: `-maxIntron=50` (or any non-default small value) activates `cutAtBigIntrons` so each copy
   is its own psl line (`supStitch.c:815-820`). Also `gfClumpsToBundles` groups clumps within 2000 bp
   (`gfBlatLib.c:649`): tandem copies closer than that are aligned together and then re-split by the
   intron cut.
7. **The 16-per-target cap**: split the genome into windows (e.g. 100 kb-1 Mb FASTA records, with unique
   names, `genoFind.c:1429-1442`) before indexing; the cap then applies per window. Or patch `ssAliCount`
   (`gfBlatLib.c:22`). gfServer+gfClient instead caps at `-maxDnaHits` clumps per query (raise it;
   `gfServer.c:1506`) and 16 per `ffIntronMax` cluster.
8. **-fine instead of fuzzyFind**: `-fine` extends the index hits directly (gapless +1/-3 with drop 9, then
   banded affine SW scoring 100/-200/300/300 limited to twice the smaller gap, `ffSeedExtend.c:500-530`,
   `:424-496`). Break-even identities: gapless extension 75 %, banded extension 66.7 %. It also re-queries the
   index in gaps (`:347-420`), useful for indel-rich copies. Downsides: splice-site refinement and polyA
   logic are mRNA-oriented, and it was written for high-quality mRNAs (`blat.c:152-153`). Worth trying;
   not verified for this use.
9. **Translated modes** (`-t=dnax -q=dnax`): seeds are 5-aa exact tiles (15 nt) with minMatch 1 and protein
   scoring (+2/-1, `ffLoose`), identity default 25 %. For non-coding repeats stop codons ('Z') break ~1 in 21
   codons per frame and every 5-aa window spanning one (`index-and-seeding.md`), and nucleotide
   substitutions change the amino acid ~2/3 of the time, so aa identity of a 75 % copy is ~45-50 % and the
   chance a given 5-aa window is exact is ~0.03-0.05. With ~100 codons per 300 bp and minMatch 1 the
   in-frame pair still seeds often, but each query/target frame pair is searched (3x3 per strand,
   `gfBlatLib.c:1294-1307`), alignments are reported in codon-multiple blocks, and `ffLoose` chaining is
   permissive. Treat as a last resort, not a sensitivity trick, for non-coding elements.

## Suggested starting configurations (untested)

| Goal | Command sketch | Expected cost (3 Gbp, model) |
|---|---|---|
| Fast, moderate sensitivity (>= 80 %) | `blat -tileSize=11 -oneOff=1 -repMatch=5000000 -minIdentity=70 -minScore=40 -maxIntron=50 -noHead chunks.fa sine.fa out.psl` | 1.1 GB index + 3 GB target; 34x lookups |
| Sensitive (>= 70 %) | `blat -tileSize=8 -stepSize=4 -repMatch=5000000 -minIdentity=60 -minScore=40 -maxIntron=50 -fine -noHead chunks.fa sine.fa out.psl` | ~3 GB entries + 3 GB target; millions of hits/query |
| Server for many queries | `gfServer start host port -tileSize=9 -stepSize=5 -repMatch=5000000 -maxDnaHits=100000 genome.2bit` then `gfClient -minIdentity=60 -maxIntron=50` | index once (~minutes per docs), but alignment still per client and 16/cluster cap |

`chunks.fa` = genome split into named windows to defeat the per-sequence cap. Verify sensitivity on a
RepeatMasker-annotated chromosome before trusting as a pre-filter; measure what fraction of annotated
copies (binned by divergence) produce a psl line.

## When BLAT is the wrong tool

BLAT has no substitution matrix, no affine gaps in its main path, no e-value model (blast8 e-values assume a
3 Gbp database and a blastz score rescale, `lib/blastOut.c:115-121`), and hard output caps. For
exhaustive repeat annotation the design intent is RepeatMasker-style tools; BLAT's design target is
near-identical mRNA/genome mapping (usage text, `blat.c:140-153`; fuzzyFind header `inc/fuzzyFind.h:8-11`).
As a pre-filter it is viable only with the index retuned (repMatch, K/step or oneOff) and the output caps
defeated (window splitting, small maxIntron, low minIdentity).
