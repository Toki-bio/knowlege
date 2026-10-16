---
title: "bwa: Control flow"
type: reference
generated: true
---

> Part of the [bwa](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/bwa/references/control-flow.md)

{% raw %}
# Control flow

## `bwa index [-a is|bwtsw|rb2] [-p prefix] [-b block] ref.fa`

1. Parse (`bwtindex.c:213-231`); default `algo_type = BWTALGO_AUTO`, `block_size = 10,000,000`.
2. `bwa_idx_build` (`bwtindex.c:255-323`):
   1. Pack FASTA **with** reverse complement -> `<prefix>.pac` (+`.ann`, `.amb`) via `bns_fasta2bntseq(fp, prefix, 0)` (`271`; revcomp appended `bntseq.c:306-312`; N -> random base `bntseq.c:266`).
   2. Choose algorithm if auto: `l_pac > 50,000,000 ? bwtsw : is` (`275`). Note `l_pac` here is already the doubled length.
   3. BWT: `bwtsw` -> `bwt_bwtgen2(pac, bwt, block_size)` (`281`, `bwt_gen.c:1592`); `is`/`rb2` -> `bwt_pac2bwt(pac, use_is)` (`284`, `bwtindex.c:64-126`: SA-IS over a 1-byte-per-base copy, `is.c:208-223`, or rope insertion `95-119`) then `bwt_dump_bwt`.
   4. Reload `.bwt`, interleave Occ checkpoints (`bwt_bwtupdate_core`, `296`, `bwtindex.c:150-172`), rewrite `.bwt`.
   5. Re-pack FASTA **forward only** -> overwrite `.pac/.ann/.amb` (`305`). So the final `.pac` is forward-only while the BWT covers both strands.
   6. Reload `.bwt`, `bwt_cal_sa(bwt, 32)` (`316`; sequential walk of all `seq_len` LF steps, `bwt.c:77-81`), `bwt_dump_sa` -> `.sa` (`bwt.c:396-407`).
3. Reads/writes: data structures `bwt_t`, `bntseq_t`.

## `bwa mem [opts] idx reads.fq [mates.fq]`

1. `main_mem` (`fastmap.c:141`): `mem_opt_init` (`157`), getopt (`159-258`), `-x` presets or `update_a` (`330-359`), `bwa_fill_scmat` (`360`).
2. Load index: shm first (`362`), else `bwa_idx_load(.., BWA_IDX_ALL)` (`364`) = `.bwt` + `.sa` + `.ann/.amb` + `.pac` into RAM.
3. Open FASTQ(s) via `kopen` (`371-392`), set `MEM_F_PE` if two files; print header (`393`).
4. `aux.actual_chunk_size = -K or chunk_size * n_threads` (`394`); `kt_pipeline(2 or 1 with -1, process, &aux, 3)` (`395`).
5. `process` (`fastmap.c:64-123`): step 0 `bseq_read` (`bwa.c:79-112`, stops after `chunk_size` bp, even count); step 1 `mem_process_seqs` (`bwamem.c:1235-1264`); step 2 write `seqs[i].sam` to stdout (`113-121`).
6. `mem_process_seqs`: per-thread `smem_aux_t` (`1249-1251`), `kt_for(worker1)` over reads (`1252`), PE insert-size inference (`1256-1259`), `kt_for(worker2)` (`1260`).
7. `worker1` -> `mem_align1_core` (`1081-1117`) per read:
   1. Convert to 2-bit (`1087-1088`).
   2. `mem_chain` (`277-341`): `mem_collect_intv` (three seeding passes, `140-188`), `frac_rep` (`291-298`), for each seed interval iterate up to `max_occ` sampled occurrences (`304-305`), `bwt_sa` -> `rbeg` (`309`), drop seeds bridging sequences/strands (`312-313`), insert into B-tree keyed by `rbeg` and try to merge into the nearest lower chain (`314-325`, `test_and_merge` `216-237`).
   3. `mem_chain_flt` (`353-411`): weight, sort, pairwise shadowing, `max_chain_extend`.
   4. `mem_flt_chained_seeds` (`624-641`): SW-tests individual seeds only for long queries (skipped when `5.5*ln(l) > 0.05*l`, i.e. for anything shorter than ~1 kb).
   5. For each chain `mem_chain2aln` (`658-812`): fetch reference window `rmax` (`668-685`), extend each seed (longest first) left/right with `ksw_extend2` unless already covered (`697-732`), band doubling once if `max_off` is large (`749-760`, `777-788`), pick local vs to-end (`762-768`, `790-796`).
   6. `mem_sort_dedup_patch` (`463-515`): remove near-identical regions (`mask_level_redun`), merge co-linear pieces (`mem_patch_reg` `432-461`), sort by score.
8. `worker2` (`1217-1233`): `mem_mark_primary_se` (`547-584`), optional `-5` reorder (`1008-1030`), `mem_reg2sam` (`1033-1079`) -> `mem_reg2aln` (CIGAR by banded global DP `1138-1152`, clipping ops `1167-1181`, MAPQ `1136`) -> `mem_aln2sam` (`851-976`). PE: `mem_sam_pe` (`bwamem_pair.c:276`).

## `bwa fastmap [-l minLen] [-w maxIntv] [-i minIntv] [-I maxIntv] [-L maxLen] idx reads`

1. Defaults `min_len=17, min_iwidth=20, min_intv=1, max_intv=0` (`fastmap.c:410-411`).
2. Loads BWT+SA+BNS, no PAC (`444`).
3. For each read: `smem_set_query`, loop `smem_next` (`bwamem_extra.c:86-96` -> `bwt_smem1a` with `min_intv`, `max_intv`), print `EM <qb> <qe> <nocc>` and positions only if `nocc <= -w` (`fastmap.c:456-473`), then `//`.
4. No re-seeding or third pass: fastmap = pass 1 of `mem_collect_intv` only.

## `bwa bwasw [opts] idx reads.fa [mates.fa]`

1. `bwa_bwtsw2` (`bwtsw2_main.c:11-89`): defaults from `bsw2_init_opt` (`bwtsw2_aux.c:54-64`); `t` and `coef` scaled by `a` (`80-81`); loads BWT+BNS (`83`).
2. `bsw2_aln` -> batches of `chunk_size*n_threads` bp (`bwtsw2_aux.c:748`), pthreads split reads round-robin (`684-699`).
3. `bsw2_aln_core` (`561-644`) per read: `update_opt` sets threshold `t = max(T, coef*ln(l))` and clamps band to query length (`545-557`); N -> random base (`587`); forward strand `bsw2_aln1_core` (`598`); reverse strand only if some hit has `n_seeds < t_seeds` (`599-602`).
4. `bsw2_aln1_core` (`248-295`): build the read's own FM-index `bwtl_seq2bwtl` (`255`); `bsw2_core(bns, opt, read_bwt, genome_bwt)` (`256`) does affine-gap DP between the read DAWG and the genome prefix trie with Z-best pruning (`bwtsw2_core.c:449-619`, heap of size `z` `467-468,544-547`, `cut_tail` `583,598`, band `491-495`), collecting hits with `G >= t` (`save_hits` `567`) and separately "narrow" hits with SA interval `<= is` (`save_narrow_hits` `582,597`, `248-271`). `bsw2_resolve_duphits` expands intervals `<= is` into every genomic copy and collapses wider ones into one representative flagged "random" (`273-310`).
5. Only narrow hits are chained (`bsw2_chain_filter`, `281`, threshold `2*t_seeds` `bwtsw2_chain.c:51`) and left-extended (`283`); all hits are right-extended (`286`) with `ksw_extend` (`133,161`); overlaps resolved by `mask_level` (`294`).
6. `write_aux` (`gen_cigar` via `bwa_gen_cigar`, MAPQ `423-429`) and `print_hits` (`491-540`): a "random"/repetitive hit prints MAPQ 0 and CIGAR `*` (`503`); extra hits are flagged secondary only with `-M` (`497`).

## `bwa aln` / `samse`

1. `bwa_aln` parses (`bwtaln.c:230-`), `bwa_aln_core` (`159-228`) loads **only `.bwt`** (`171-175`), reads 0x40000 reads per batch (`180`), threads by `i % n_threads` (`102`).
2. Per read: `bwt_cal_width` lower bounds (`57-81`), `max_diff` from read length and `-n` (`111`), `bwt_match_gap` (`bwtgap.c:109`) bounded backtracking over the FM-index, stack capped at `max_entries` (`bwtgap.c:139`); writes `(n_aln, bwt_aln1_t[])` to `.sai` (`214-218`).
3. `samse -n N` (`bwase.c:581-589`): `bwa_aln2seq_core` emits XA multi-hits only if total occurrences `<= N+1`, else none (`bwase.c:57-58`).

## Threading / parallelism model

- `kt_for(n, func, data, N)` (`kthread.c:49-61`): `n` pthreads, each starts at index `tid` and strides by `n`; when done it steals from the least-advanced worker (`steal_work`, `25-33`). Work item = one read (or one pair). Per-thread scratch `smem_aux_t` via `tid` (`bwamem.c:1208`).
- `kt_pipeline(n_threads=2, func, shared, n_steps=3)` (`119-147`): workers advance batches through steps in order; a batch may enter step `s` only if no batch with a smaller index is at step `<= s` (`96-101`). With 2 workers, reading batch i+1 overlaps aligning batch i; `-1` forces 1 worker (`fastmap.c:395`).
- Output order equals input order (step 2 is serialized by the index rule).
- Memory per batch: `chunk_size * n_threads` bases of reads plus all `mem_alnreg_v` for the batch (README FAQ 8 warns coordinate-sorted input makes batches of repeats blow up time/memory).
- `bwasw` and `aln` use their own ad-hoc pthread fan-out (`bwtsw2_aux.c:673-706`, `bwtaln.c:186-205`).

## Error handling conventions

`xopen/err_*` wrappers in `utils.c` abort on I/O error; `xassert` for invariants (`bwt.c:68,429,433`); `bwa_verbose` levels 1-4 (`-v`), level 4 prints seeds/chains/extension traces to stdout and breaks SAM (`bwa.1:364-369`). Return codes: `main` prints version/CMD/time on success (`main.c:121-127`).
{% endraw %}
