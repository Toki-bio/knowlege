# One tig through utgcns (commit 86d7d3b88)

## Top loop (processTigs, utgcns-processTigs.C:92-173)

1. `params.loadPartitionedReads()` - all partition reads into RAM
   (utgcns-parameters.C:30-79). RSS ~ partition size + slack.
2. For each tig (line 123): `tig->filterContains(params.maxCov, false)`
   (line 132) caps evidence per window BEFORE consensus.
3. `unitigConsensus::generate(tig, 'P', aligner, iter, reads)` (line 140).
   On success: saveToStream to `-O` output; on failure: stderr warning; either
   way the progress-row second half prints (line 135) and the next tig loads.
4. Footer + "Consensus finished successfully" / WARNING per numFailures.

## generate() (unitigConsensus.C:1732-1790)

1. Copy child placements: _utghpc (homopolymer-compressed), _utgpos (working),
   _cnspos/_adjpos (outputs); classify low-quality reads by ident <= errorRateMaxID.
2. `switchToUncompressedCoordinates()`.
3. Dispatch: 1 usable read -> generateSingleton; 'Q' -> generateQuick (template
   only); 'P'/'p' -> generatePBDAG (the normal path, line 1098).
4. After: findCoordinates, trimCircular.

## generatePBDAG (1098-~1200)

1. `generateTemplateStitch()` (line 1105) - build a linear template by
   walking children in layout order (see below).
2. Align EVERY read to the finished template (edlibAlign at line 1496,
   CIGAR built at 1543).
3. Build partial-order multialignment / DAG from those alignments
   (line ~1156, "not thread safe"), vote per column, consensus string.

## generateTemplateStitch detail (386-860), the retry/backtrack machinery

Per read 'nr' walking the tig:

- `alignAgain:` (line 586): edlib template-tail vs read-prefix, HW mode,
  k = olapLen * bandErrRate (610-612). On noResult / hitTheStart:
  templateSize -= 0.10; on hitTheEnd: extensionSize += 0.10.
- When templateSize < 0.01 (698): either accept failure (k exhausted) or
  bandErrRate += _errorRate/ERROR_RATE_FACTOR and reset sizes (702-714)
  -> `goto alignAgain` (717-720).
- On final failure (noResult at 725): backtrack blocks:
  - 738-764: mark nr bad (`badToAdd`); if not back at firstCandidate or if
    !allowLargerShift -> flip allowLargerShift, CLEAR badToAdd (752),
    `goto retryCandidate` (763) to re-place earlier reads with looser shift.
  - 765-792: else trim trimbp off the template end, reset error rate and
    badToAdd, `goto retryCandidate` (792). Repeatable while template > minOlap.
  - 795-809: give up on alignment, append read verbatim to template.
- Whole-tig restart `retryTig` (522-528): if reads are exhausted before the
  tig end, savedErrorRate *= 1.25 and the entire template is rebuilt from
  scratch. Bounded only by _errorRateMax (assert at 524).

Reads that fail to extend the template are still appended (823-825), so
stitch "success" can be a chimera template; per-column voting can still
mismatch writes - see pitfalls.md.

## Logging (where progress goes)

- stdout rows (consensus.NNNN.out): first half at tig load (skipTig),
  second half after generate returns. Format:
  `tigID length reads used coverage ignored coverage doveUsed doveCoverage`.
- showPlacement()/showAlgorithm() are compile-time-ish debug prints (off by
  default); the alignAgain inner retries are therefore SILENT in production
  logs - a stall shows nothing but CPU.
