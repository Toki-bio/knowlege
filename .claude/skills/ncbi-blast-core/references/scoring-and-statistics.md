# Scoring and statistics (`blast_stat.c`, `blast_setup.c`, `blast_parameters.c`)

## 1. Nucleotide scoring matrix

`BlastScoreBlkNuclMatrixCreate` (`blast_stat.c:1060-1130`): 16x16 blastna
matrix; A/C/G/T pairs get `reward` or `penalty`; ambiguity codes get the
degeneracy-weighted average `((degen-1)*penalty + reward)/degen` rounded
(:1105-1113); row/column 15 (gap/sentinel) = `INT4_MIN/2` (:1123-1127).
Default macros: `BLAST_REWARD 1`, `BLAST_PENALTY -3` (`blast_options.h:151-152`)
but the task presets override (megablast 1/-2, blastn 2/-3, see task-presets.md).

## 2. Karlin-Altschul tables for reward/penalty pairs

Each row: `{gap_open, gap_extend, Lambda, K, H, alpha, beta, <8th col>}`
(indices used at :3852-3856 and :3970-3973). Row `{0,0,...}` = ungapped and
linear (megablast 0/0) values (`s_SplitArrayOf8` :3153). Tables
(`blast_stat.c`): 1/-5 :611, 1/-4 :617, 2/-7 :629, 1/-3 :638, 2/-5 :651,
1/-2 :660, 2/-3 :674, 3/-4 :687, 4/-5 :697, 1/-1 :706, 3/-2 :717, 5/-4 :722.
Reward/penalty are first divided by their gcd (:3248-3254), so 2/-6 == 1/-3.

Selected rows (Lambda, K, H):

| reward/penalty | ungapped or linear (0/0) | gapped rows available (open,extend) |
|---|---|---|
| 1/-2 (megablast) | 1.28, 0.46, 0.85 (:661) | 2/2, 1/2, 0/2, 3/1, 2/1, 1/1 (:662-667); max 2/2 |
| 2/-3 (blastn, dc-megablast) | 0.55, 0.21, 0.46 (:675) | 4/4, 2/4, 0/4, 3/3, 6/2, **5/2 = 0.625, 0.41, 0.78** (:681), 4/2, 2/2 (:676-683); max 6/4 |
| 1/-3 (blastn-short) | 1.374, 0.711, 1.31 (:639) | 2/2, 1/2, 0/2, 2/1, 1/1 (:640-644); max 2/2 |
| 1/-1 | none | 3/2 ... 2/1 (:707-713); max 4/2 |

`round_down` is TRUE for 2/-7, 2/-5, 2/-3, 3/-4 (:3278, :3300, :3308, :3316):
odd gapped scores are rounded down to even for E-value purposes
(`Blast_HSPListGetEvalues` `blast_hits.c:1858-1861`).

Allowed gap costs: `Blast_KarlinBlkNuclGappedCalc` (:3847-3954) looks the
(open,extend) pair up; 0/0 uses the linear row (:3882-3888); costs
`>= gap_open_max && >= gap_extend_max` copy the ungapped block (:3908-3910);
anything else is an error listing the supported values (:3911-3947).
`BLAST_GetNucleotideGapExistenceExtendParams` (:3403-3453) *silently*
replaces unsupported costs that are below the maxima by the maxima
(:3443-3446). Unsupported reward/penalty -> "Substitution scores %d and %d are
not supported" (:3340-3348). Validation in `BlastScoringOptionsValidate`
(`blast_options.c:862-954`): penalty must be negative unless 0/0 (rmblastn);
gap_open > 0 with gap_extend 0 is rejected (:897-902); the reward/penalty check
itself is disabled with `if (0 && ...)` (:890).

## 3. Karlin blocks

- Ungapped `kbp_std` from score frequencies (`Blast_ScoreBlkKbpUngappedCalc`
  :2738; `Blast_KarlinBlkUngappedCalc` :2700).
- Gapped `kbp_gap_std` per context from the tables
  (`Blast_ScoreBlkKbpGappedCalc`, `blast_setup.c:41-134`, blastn branch
  :82-102); Gumbel block for blastn explicitly disabled (:54-62), so `sbp->gbp`
  is NULL and `BLAST_SpougeStoE`/`EtoS` are never used for nucleotide.
- rmblastn (reward=penalty=0) borrows the 1/-3 statistics (:91-96).

## 4. Effective search space and length adjustment

`BLAST_CalcEffLengths` (`blast_setup.c:699-850`): `db_length` = `-dbsize`
override or real total (:711-714), `db_num_seqs` = `-num_seqs`/real;
`Blast_GetNuclAlphaBeta` (`blast_stat.c:3966-4035`, alpha/beta from the
gapped row or `Lambda/H`, beta -2 for 1/-1 and 2/-3 ungapped :3956-3963);
`BLAST_ComputeLengthAdjustment` (:5042-5127) iterates
`ell = alpha/lambda*(logK + log((m-ell)(n-N*ell))) + beta`;
`eff_searchsp = (db_length - N*ell) * (qlen - ell)` (:826-839), floored at
`db_length_eff >= 1`. `-searchsp` overrides everything (:719-720, `blast_args.cpp:300-304`).
For `-subject` in legacy mode (`db_length == 0`) this is recomputed per
subject (`blast_engine.c:1430-1439`, `BLAST_OneSubjectUpdateParameters`
`blast_setup.c:1001-1027`); in default dbscan mode total subject length acts
as database length (`seqsrc_multiseq.cpp:176-182`).

## 5. Cutoff scores from E-value

`BlastHitSavingParametersUpdate` (`blast_parameters.c:831-996`): per context
`BLAST_Cutoffs(&S, &E, kbp_gap, eff_searchsp, FALSE, 0)` (:936) where
`S = ceil(log(K*searchsp/E)/Lambda)` (`BlastKarlinEtoS_simple` `blast_stat.c:4041`,
`BLAST_Cutoffs` :4091-4155); `-min_raw_gapped_score` replaces it (:876-896).
`cutoff_score_max` is the E-value cutoff; sum-statistics may lower
`cutoff_score` (:943-968, not for gapped blastn). The ungapped stage then uses
`min(gap_trigger, cutoff_score_max)` (see extension-and-hsp-saving.md).
E-value of an HSP: `E = searchsp * exp(-Lambda*S + logK)`
(`BLAST_KarlinStoE_simple` :4158-4173) with the context's `eff_searchsp`
(`blast_hits.c:1878-1883`); bit score `(Lambda*S - logK)/ln2`.

Notes (arithmetic on the cited constants, not executed): a 300 bp query vs a
3 Gb genome at E=10 needs raw ~39 with blastn 2/-3 5/2 (~36 bits); megablast
1/-2 needs raw ~19 (19 identities); blastn-short at E=1000 needs ~15.

## 6. Corrections applied for nucleotide

- Length adjustment only (above). No composition-based statistics:
  `compositionBasedStats = eNoCompositionBasedStats` unless PSSM
  (`blast_options.c:682-688`) and `CCompositionBasedStatsArgs` is not among
  blastn arguments (`blastn_args.cpp:46-121`). `blast_kappa.c` is protein-only.
- Ambiguity re-evaluation after traceback (`Blast_HSPReevaluateWithAmbiguitiesGapped`).
- rmblastn complexity adjustment (`blast_traceback.c:521-583`) and
  `-mask_level` (`Blast_HSPResultsApplyMasklevel`, `blast_hits.c:3467`).
