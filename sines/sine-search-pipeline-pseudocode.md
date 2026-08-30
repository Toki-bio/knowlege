---
title: SINE search pipeline pseudocode
type: explanation
---

> Type: explanation

# Multi-genome SINE search pipeline — pseudocode

This is a high-level record of a **Toki / SINEderella-style QuickSearch** workflow:

search many SINE consensuses against many genomes → extract flanks → cluster / align → rebuild consensuses → report.

It is a **targeted library search** of known or candidate SINE models, not a de novo TE discoverer like HiTE, and not a mappability tool like Newmap. See also [Newmap / HiTE / NSIT roles](https://toki-bio.github.io/knowlege/methods/newmap-hite-nsit-for-sine-search/) and [SINE-de-novo-genome-scan](https://github.com/Toki-bio/SINE-de-novo-genome-scan) for related but different pipelines.

---

## Inputs / outputs

```text
INPUT
  query_consensuses.fa     # SINE family / subfamily consensuses
  genome_list.txt          # one assembly path per line
  flanking_length          # e.g. 50–100 bp
  min_match_length
  min_identity
  out_dir/

OUTPUT
  out_dir/<genome>/<consensus>/
      hits.bed
      search.log
      flanks.fa
      clusters/
      alignments/
      rebuilt_consensuses.fa
  out_dir/summary.tsv
  out_dir/warnings.txt
```

---

## Compact one-page pseudocode

```text
PIPELINE QuickSearch_SINE_discovery

LOAD queries  ← read_fasta(query_consensuses.fa)
LOAD genomes  ← read_lines(genome_list.txt)

FOR each genome IN genomes:
    FOR each consensus IN queries:
        hits, log ← sensitive_local_search(
            query   = consensus,
            target  = genome,
            min_len = min_match_length,
            min_id  = min_identity)

        RECORD hit_count(consensus, genome) FROM log   # "Found X hits"
        SAVE hits AS bed

        IF hits empty: CONTINUE

        flanks ← []
        FOR each interval IN hits:
            flanks ← flanks + extract(genome, interval ± flanking_length)
        SAVE flanks.fa

        clusters ← identity_cluster(flanks.fa, thresholds = 95% .. 75%)
        FOR each cluster:
            aln ← multiple_align(cluster)
            cons ← plurality_consensus(aln)
            STORE aln, cluster_size, cons

WRITE summary.tsv
WRITE sibling-similarity warnings
END
```

---

## Expanded modular pseudocode

```text
BEGIN PIPELINE  "QuickSearch SINE discovery"

############################################################
# CONFIG
############################################################
query_consensuses_fasta ← ...
genome_list_file        ← ...
flanking_length         ← 50..100
min_match_length        ← ...
min_identity            ← ...
cluster_thresholds      ← [0.95, 0.90, 0.85, 0.80, 0.75]
output_directory        ← ...

############################################################
# LOAD
############################################################
queries ← READ_FASTA(query_consensuses_fasta)
genomes ← READ_LINES(genome_list_file)

ASSERT each genome file exists and is indexed if required
ASSERT query IDs are unique and filesystem-safe

############################################################
# PER-GENOME SEARCH
############################################################
FOR each genome IN genomes:

    genome_id ← basename(genome)
    LOG("Searching genome:", genome_id)

    # ------------------------------------------------------
    # Step 1. Search every consensus against this genome
    # ------------------------------------------------------
    FOR each consensus IN queries:

        hit_bed ← TEMP(".bed")
        log_file ← TEMP(".log")

        RUN sensitive_local_search(
                query         = consensus.sequence,
                target_genome = genome,
                min_len       = min_match_length,
                min_id        = min_identity,
                output_bed    = hit_bed,
                log           = log_file)
        # Typical engine in this ecosystem: ssearch36 / similar
        # sensitive local alignment, not exact k-mer uniqueness.

        hit_count ← PARSE_HIT_COUNT(log_file)   # prefer log "Found X hits"
        RECORD metrics[genome_id][consensus.id].n_hits ← hit_count

        STORE search_results[genome_id][consensus.id].bed ← hit_bed
        STORE search_results[genome_id][consensus.id].log ← log_file

    END FOR

    # ------------------------------------------------------
    # Step 2. Extract flanking sequences around each hit
    # ------------------------------------------------------
    FOR each consensus IN queries:

        hits ← search_results[genome_id][consensus.id].bed
        IF EMPTY(hits): CONTINUE

        extracted ← EMPTY_LIST

        FOR each interval IN hits:
            # Keep strand; pad safely at contig ends
            window ← EXPAND(interval, flanking_length, clipped_to_contig=TRUE)
            seq ← EXTRACT_FASTA(genome, window)
            seq.header ← MAKE_HEADER(genome_id, consensus.id, interval, strand)
            APPEND seq TO extracted
        END FOR

        SAVE extracted TO
            output_directory / genome_id / consensus.id / "flanks.fa"

        RECORD metrics[...].n_extracted ← COUNT(extracted)

    END FOR

    # ------------------------------------------------------
    # Step 3. Cluster → align → rebuild consensus
    # ------------------------------------------------------
    FOR each consensus IN queries:

        flank_file ← ... / "flanks.fa"
        IF COUNT_SEQUENCES(flank_file) == 0: CONTINUE

        # Sweep identity thresholds; finer first
        FOR each thr IN cluster_thresholds:

            clusters ← VSEARCH_CLUSTER(
                flank_file,
                id = thr,
                strand = both)   # or plus-only if already stranded

            FOR each cluster IN clusters:
                IF cluster.size < min_cluster_size: CONTINUE

                aln ← MAFFT(cluster.sequences)   # or preferred MSA tool
                cons ← BUILD_CONSENSUS(aln, plurality_pct)

                STORE:
                    alignment
                    cluster_size
                    threshold
                    rebuilt_consensus

                # Optional quality flags
                FLAG if cons is much shorter than query
                FLAG if cons is nearly identical to a sibling query
                FLAG if cluster is tiny / unstable across thresholds

            END FOR
        END FOR

    END FOR

END FOR   # genome loop

############################################################
# Step 4. Cross-genome report
############################################################
GENERATE REPORT:
    FOR each genome:
        FOR each consensus:
            PRINT:
                total hits
                total extracted sequences
                cluster distribution by threshold
                consensus stability metrics
                warnings if rebuilt consensus collides with sibling queries

WRITE:
    summary.tsv
    per-genome logs
    per-consensus hit tables
    extracted FASTA sets
    rebuilt consensus library

END PIPELINE
```

---

## Why the structure looks like this

### Multi-genome outer loop

One assembly at a time keeps I/O and temp space manageable and makes per-genome abundance comparable.

### Multi-consensus inner search

Each consensus is searched independently so hit tables stay attributable to a named model / subfamily query. That matches SINEderella-style “search this subfamily consensus across genomes” work.

### Search → flanks → cluster → align → consensus

Hits alone are coordinates. Flanks recover insertion context and more complete copies. Clustering + MSA turns raw hits into inspectable sequence groups and rebuilt models.

### Log-driven hit counts

In the real runs, counts often come from search logs (`Found X hits`) rather than recounting BED lines. That is fine if logs and BEDs are guaranteed to match; otherwise validate once that they do.

### Threshold sweep (95% → 75%)

High identity finds near-identical young copies; lower identity reveals broader families / older clouds. Do not treat one cutoff as the biological truth — see [HDBSCAN vs vsearch](https://toki-bio.github.io/knowlege/methods/hdbscan-vs-vsearch-for-sine-grouping/).

---

## Suggested elaborations / improvements

These are optional upgrades that fit the same architecture.

### 1. Make `sensitive_local_search` an explicit interface

Document the actual engine and filters:

```text
sensitive_local_search =
    ssearch36 (or equivalent)
    + min length / identity / coverage filters
    + BED conversion
```

That avoids later confusion with RepeatMasker, HiTE, or exact-k-mer tools.

### 2. Deduplicate overlapping hits before extraction

Overlapping / nested intervals from related consensuses can inflate flank sets. A merge or “best-by-score per locus” step before extraction often cleans abundance stats.

### 3. Separate locus abundance from copy abundance

Report both:

* number of genomic intervals / loci,
* number of extracted monomers / clustered copies.

Huge arrays and scattered singletons are different biology.

### 4. Sibling-consensus collision checks are important

If two query consensuses are very close, their hit sets and rebuilt models will bleed into each other. The report should warn when:

```text
rebuilt_cons(A) ≈ query(B)
```

or when hit intervals for A and B overlap heavily.

### 5. Keep original queries as the archive

Rebuilt consensuses are derived products. Do not overwrite the input library. That matches the [MSA-as-catalog](https://toki-bio.github.io/knowlege/sines/msa-as-sine-consensus-catalog/) principle: sequences + alignments remain source of truth; labels and rebuilt models are layers.

### 6. Optional fourth stage: profile refinement

After accepting a rebuilt consensus:

```text
build HMM → nhmmer / hmmsearch genome
```

Useful for ultra-degenerate copies the primary local search missed. Keep it optional so the core pipeline stays fast.

### 7. Optional demultiplexing note

If one query pulls mixed structural variants, clustering at multiple thresholds is only a first split. Hard cases may need MSA inspection or density-based audits rather than a single vsearch cut.

---

## ASCII block diagram

```text
query_consensuses.fa          genome_list.txt
         │                            │
         └────────────┬───────────────┘
                      ▼
              FOR each genome
                      │
          ┌───────────┴───────────┐
          ▼                       ▼
   FOR each consensus      (parallelizable)
          │
          ▼
  sensitive local search ──► hits.bed + log
          │
          ▼
   extract ± flanks ──► flanks.fa
          │
          ▼
   vsearch 95%..75%
          │
          ▼
   MAFFT + consensus
          │
          ▼
   per-genome / per-consensus products
                      │
                      ▼
                 summary.tsv
            + sibling-collision warnings
```

---

## What this pipeline is / is not

| Is | Is not |
| --- | --- |
| Targeted multi-genome search with known/candidate SINE models | De novo TE library building (HiTE / RepeatModeler) |
| Hit extraction + identity clustering + consensus rebuild | Species-tree phylogenetics |
| Practical discovery / abundance / model-update loop | Exact unique k-mer mappability (Newmap) |
| Compatible with later HMM refinement | Automatic final taxonomy |

---

## Key takeaways

* Outer loop = genomes; inner loop = consensuses; then flanks → cluster → align → rebuild.
* Prefer log-validated hit counts, but keep BED as the coordinate truth.
* Sweep identity thresholds; do not worship one cutoff.
* Warn when rebuilt models collide with sibling queries.
* Treat rebuilt consensuses as derived layers; keep the input library intact.
* Optional HMM rescue and hit-dedup are the highest-value extensions without changing the architecture.
