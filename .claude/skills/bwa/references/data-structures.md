# Core data structures

## `bwt_t` (`bwt.h:48-60`) — the FM-index

| Field | Type | Meaning | Invariant |
|---|---|---|---|
| `primary` | `bwtint_t` (uint64) | SA position of the `$` row | removed from the BWT text; `k -= (k >= primary)` everywhere (`bwt.c:114`) |
| `L2[5]` | uint64 | cumulative counts C(c); `L2[4] == seq_len` (`bwt.c:457`) | |
| `seq_len` | uint64 | length of indexed text = `2 * l_pac` (fwd + revcomp) | |
| `bwt_size` | uint64 | number of uint32 words in `bwt` incl. interleaved Occ | `seq_len/16 + 8*(seq_len/128+1)` (`bwtindex.c:74,155-156`) |
| `bwt` | uint32* | per 128 bases: 4 x uint64 Occ counts (8 words) then 8 words of 2-bit text (`bwt.h:74-75`) | layout fixed to OCC_INTERVAL=128 (`bwt.h:36`) |
| `cnt_table[256]` | uint32 | byte -> per-base popcounts (`bwt_gen_cnt_table`, `bwt.c:42-51`) | computed on load |
| `sa_intv`, `n_sa`, `sa` | int, uint64, uint64* | sampled SA: `sa[k/intv]`, `n_sa = (seq_len+intv)/intv` (`bwt.c:73`) | `intv` power of 2 (`bwt.c:68`); `sa[0] = -1` sentinel (`bwt.c:83`) |

Created by: `bwt_pac2bwt` + `bwt_bwtupdate_core` + `bwt_cal_sa` (index) or `bwt_restore_bwt` + `bwt_restore_sa` (`bwt.c:443,421`). Freed by `bwt_destroy`. Read-only during alignment (shared across threads).

## `bwtintv_t` (`bwt.h:62-66`) — a bidirectional SA interval = one exact match

| Field | Meaning |
|---|---|
| `x[0]` | SA interval start of the match |
| `x[1]` | SA interval start of its reverse complement |
| `x[2]` | interval size = number of genomic occurrences (both strands) |
| `info` | `(query_start << 32) | query_end` (`bwt.c:302,334`) |

`bwtintv_v` is a kvec of these. `bwt_set_intv` (`bwt.h:82`) initialises a 1-base interval. `smem_aux_t` (`bwamem.c:119-121`) holds per-thread scratch vectors.

## `bntseq_t`, `bntann1_t`, `bntamb1_t` (`bntseq.h:41-64`) — reference annotation

- `l_pac`: forward genome length in bases. The FM-index covers `2*l_pac`; positions `>= l_pac` are reverse-strand, mapped back by `bns_depos` (`bntseq.h:87-90`).
- `anns[n_seqs]`: per sequence `offset`, `len`, `n_ambs`, `gi`, `is_alt` (from `<prefix>.alt`, `bntseq.c:176-209`), `name`, `anno` (FASTA comment).
- `ambs[n_holes]`: runs of non-ACGT bases (`offset`, `len`, `amb` char), built in `add1` (`bntseq.c:246-263`).
- `seed = 11`: srand48 seed for N-replacement (`bntseq.c:295`).

## `bwaidx_t` (`bwa.h:48-56`)

`bwt`, `bns`, `pac` (2-bit forward genome, `l_pac/4+1` bytes, `bwa.c:308`), plus `mem/l_mem/is_shm` when the whole index lives in one shm blob (`bwa_mem2idx`, `bwa.c:337-360`).

## `mem_opt_t` (`bwamem.h:52-84`) — all `bwa mem` parameters

Defaults from `mem_opt_init` (`bwamem.c:74-110`); CLI mapping in `main_mem` (`fastmap.c:160-258`).

| Field | Flag | Default | Used at |
|---|---|---|---|
| `a`, `b` | `-A`, `-B` | 1, 4 | `bwa_fill_scmat` (`bwa.c:136-145`; N scores -1) |
| `o_del,e_del,o_ins,e_ins` | `-O`, `-E` | 6, 1, 6, 1 | ksw calls |
| `pen_clip5/3` | `-L` | 5, 5 | end bonus in `ksw_extend2` (`bwamem.c:757,785`) |
| `w` | `-w` | 100 | band; chaining diagonal tolerance (`bwamem.c:228`) |
| `zdrop` | `-d` | 100 | `ksw.c:494-500` |
| `T` | `-T` | 30 | output filter `bwamem.c:1048` |
| `min_seed_len` | `-k` | 19 | `bwamem.c:153,166,176,286,384,984` |
| `split_factor` | `-r` | 1.5 | re-seed length `split_len = k*r` (`bwamem.c:144`) |
| `split_width` | `-s` (hidden) | 10 | re-seed only if occ <= this (`bwamem.c:163`) |
| `max_occ` | `-c` | 500 | occurrence subsampling (`bwamem.c:304-305`), `frac_rep` (`291-298`) |
| `max_mem_intv` | `-y` | 20 | third-pass seeding (`bwamem.c:170-185`, `bwt.c:358-379`) |
| `max_chain_gap` | `-G` | 10000 | chaining (`228`), filtering (`381`), dedup (`471-472`) |
| `drop_ratio` | `-D` | 0.5 | chain filter (`384`) and secondary output (`1050`) |
| `min_chain_weight` | `-W` | 0 | `bwamem.c:363`, seed-SW filter `626` |
| `max_chain_extend` | `-N` | 1<<30 | `bwamem.c:399-404` |
| `mask_level` | `-X` (hidden) | 0.5 | chain/alignment overlap test (`381`, `534`) |
| `mask_level_redun` | — | 0.95 | dedup (`481`) |
| `XA_drop_ratio` | `-z` | 0.8 | XA eligibility (`bwamem_extra.c:119`) |
| `max_XA_hits`, `_alt` | `-h` | 5, 200 | all-or-nothing cap (`bwamem_extra.c:144`) |
| `mapQ_coef_len`, `_fac` | `-Q` (hidden) | 50, ln 50 | MAPQ length factor (`bwamem.c:992-996`) |
| `n_threads`, `chunk_size` | `-t`, (`-K`) | 1, 10,000,000 bp | batch = chunk*threads (`fastmap.c:394`) |
| `pen_unpaired`, `max_ins`, `max_matesw` | `-U`, —, `-m` | 17, 10000, 50 | PE only |
| `flag` | `-a -P -p -M -S -Y -V -5 -q -u` | 0 | `MEM_F_*` (`bwamem.h:40-50`) |
| `mat[25]` | — | filled from a,b | 5x5 with N row/col = -1 |

`update_a` (`fastmap.c:125-139`): if `-A` is given, every unspecified penalty (`B T O E d L U`) is multiplied by A. `-x pacbio|ont2d|intractg` presets at `fastmap.c:330-359`.

## `mem_seed_t`, `mem_chain_t` (`bwamem.c:194-206`)

Seed: `rbeg` (position in the 2*l_pac concatenated space), `qbeg`, `len`, `score` (= len unless seed-SW ran). Chain: seed array, `rid`, `w` (weight = non-overlapping seed bases, min of query/ref projections, `mem_chain_weight` `239-258`), `kept` (0 drop / 1 shadowed-first / 2 overlapping / 3 clean), `is_alt`, `frac_rep`, `pos` (first seed `rbeg`, B-tree key `212-213`).

## `mem_alnreg_t` (`bwamem.h:86-104`) — one extended alignment region

`rb,re,qb,qe,rid`, `score` (best local), `truesc` (score of the reported region, may be lower when to-end extension chosen), `sub` (best overlapping competitor score), `csub` (tandem/sub chain), `sub_n`, `w` (band used), `seedcov`, `secondary` (index of the parent primary, <0 if primary), `secondary_all`, `seedlen0`, `n_comp`, `is_alt`, `frac_rep`, `hash` (tie-break). Produced by `mem_chain2aln` (`658-812`), culled by `mem_sort_dedup_patch` (`463-515`), classified by `mem_mark_primary_se` (`547-584`).

## `mem_aln_t` (`bwamem.h:114-124`) — SAM-ready record

`pos` (forward-strand, 0-based within `rid`), `rid`, `flag`, `is_rev`, `mapq`, `NM`, `cigar` (BAM encoding, MD string appended after the ops: `bwa.c:198,223`), `XA`, `score`, `sub`, `alt_sc`. Built by `mem_reg2aln` (`bwamem.c:1119-1189`).

## `bseq1_t` (`bwa.h:58-61`)

`l_seq, id, name, comment, seq, qual, sam`. `seq` is converted in place to 2-bit codes by `mem_align1_core` (`bwamem.c:1087-1088`); `sam` receives the output lines.

## BWA-SW and backtrack

- `bsw2opt_t` (`bwtsw2.h:14-20`), defaults `bwtsw2_aux.c:54-64`: a=1 b=3 q=5 r=2 t=30 bw=50 z=1 is=3 t_seeds=5 mask_level=0.5 coef=5.5.
- `bsw2hit_t` (`bwtsw2.h:22-27`): `k,l` SA interval (after `bsw2_resolve_duphits`, `l==0` and `k` is a genome position), `G` score, `G2` second best, `beg,end` on query, `n_seeds`, `is_rev`, `flag` (bit 1 = "random hit" from a wide interval, `bwtsw2_core.c:304`).
- `bwtl_t` (`bwt_lite.h:6-11`): 32-bit FM-index of a single read (`bwtl_seq2bwtl`, `bwt_lite.c:13`).
- `gap_opt_t` (`bwtaln.h:105-115`), defaults `bwtaln.c:24-40`: s_mm=3 s_gapo=11 s_gape=4, fnr=0.04 (-> max_diff via `bwa_cal_maxdiff`, `bwtaln.c:42-54`), max_gapo=1, max_gape=6, seed_len=32, max_seed_diff=2, max_entries=2,000,000, max_top2=30.
