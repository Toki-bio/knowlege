---
title: Orchestrating a limited-context model (GLM via glm.js / aider) on code bigger than its window
date: 2026-10-17
tags: [llm, orchestration, glm, aider, code-audit, context-window]
---

# Orchestrating a limited-context model on code bigger than its window

Setting: a strong orchestrator (Claude) with a scarce budget, and a cheap worker
(glm-5.3 on `api-llm.gpu.uz`, 262k nominal context, one request at a time per key)
driven through `glm.js` (read-only, windowed `read_file`, `grep`) or `aider` (edits).
The concrete trigger: ViewAlign's `script.js` is ~416k tokens, so aider cannot load it at
all, and a whole-file prompt would be useless even if it fit.

## What the evidence says

1. **The usable window is far smaller than the advertised one.** RULER defines effective
   length as the longest input where a model still scores 85 %; for most models it is about
   half the claimed length or less, and most models start degrading between 16k and 64k
   tokens. NoLiMa (no lexical overlap between question and answer) is harsher: at 32k, 11
   models fall below half their short-context score. Material in the middle of a long prompt
   is used worst. Consequence: budget **30-60k tokens per call** whatever the window says,
   and put the instruction and the code that matters at the start or end.
2. **Split by program unit, not by byte count, and navigate on demand.** RepoAudit (ICML
   2025) audits whole repositories one function at a time, follows data flow to the next
   function only when needed, carries findings between functions in a small memory, and runs
   a separate validator over every candidate bug: 78 % precision at about $2.5 per project.
   SWE-agent found that a windowed file viewer (about 100 lines, with line numbers, search
   that returns at most 50 hits) beats showing the whole file, and that whole-file views make
   the model lose focus.
3. **Map-reduce needs a protocol, or the joins are lost.** LLM×MapReduce splits long input
   into chunks, but has each chunk return a fixed structure (extracted information,
   rationale, answer, calibrated confidence) so the reduce step can join facts that span
   chunks and settle conflicts. Without it, cross-chunk dependencies and conflicts are where
   the errors come from.
4. **The worker's self-report is not evidence; verify every finding independently.** For
   bug reports from an LLM: majority voting over several samples, a cross-model verifier,
   and above all a mechanical oracle (a test, a reproducer, a grep for the quoted line)
   remove most false positives (Tencent industrial study: 94-98 % of false positives removed
   with hybrid LLM + static checks). Same-model resampling misses errors the model makes
   consistently, so the last check should not be the same model again. Long, demanding
   verifier prompts increase wrong rejections; keep them narrow.
5. **Orchestrator-worker:** the orchestrator decomposes, keeps the task list, briefs each
   worker narrowly, and merges; workers do one sub-task, return a compact result and forget.
   Decomposition errors propagate to every worker, so the decomposition is the orchestrator's
   job, not the cheap model's.

## Best practice for glm.js + aider (what to actually do)

**Reading tasks (audits, reviews, "does X hold"):**

- Slice by unit: a function, a nested function, or the diff hunks of one feature. Give the
  worker **line ranges and anchors**, not the file: "read lines 5910-5997 (draw) and
  5998-6026 (redrawCells); grep for `contentSig` to find its definition". `glm.js` already
  has `read_file` with `offset`/`limit` and `grep`, which is the SWE-agent viewer.
- Put the **contract** in the task: what the code must guarantee ("a copied frame equals a
  full redraw pixel for pixel", "a click on a name never selects a cell"). Workers find bugs
  against a stated invariant far better than "find bugs".
- Ask for a **fixed output schema** per finding: title, file:line, the exact quoted line,
  a concrete failure scenario (inputs, then wrong result), severity, confidence 1-5. This is
  the map-side protocol that makes merging possible.
- Keep each task to **1-3 claims or one area**, under ~40k tokens of reading. Reference
  facts the worker needs from other areas inline (two or three lines), instead of letting it
  wander.
- **Reduce in the orchestrator, mechanically first:** drop findings whose quoted line does
  not exist (grep it); dedupe by location; then verify the rest.
- **Verify with an oracle before believing:** a reproducer in the browser harness or a
  failing regression check. A finding without a reproducible failure scenario is a
  question, not a bug. Earlier runs in this harness quoted code that does not exist and
  gave wrong reasons for right quotes, so check both.
- The endpoint allows one request at a time per key: run tasks **sequentially with retry
  and back-off** (`run_*_tasks.sh` pattern), unattended in the background.

**Editing tasks with aider when the target file is too big:**

- **Extract, edit, splice.** Copy the target function(s) into a small scratch module file
  (with a header comment listing the globals they use), give aider only that file
  (`--file scratch.js`, `--read` for one or two small reference files), then splice the
  result back by exact start/end anchors and run the syntax check plus the regression gate.
  Never `--file script.js`.
- Always `--edit-format diff --map-tokens 0 --no-check-update`, run in its own git
  worktree, with `PYTHONUTF8=1` on Windows.
- Gate every run with a deterministic check (`aider-loop.sh` + `BROWSER_CHECK_CMD`); the
  wrapper, not the model, decides "done".
- Long-term fix for the root cause: split `script.js` into modules of a few thousand lines
  (aider's own advice for large projects); every tool, cheap or expensive, works better on
  small files.

**What to keep for the strong model:** decomposition, the contracts, the final verification
of findings, and any change to code with pixel-exact or concurrency invariants (here: the
Canvas frame copy, hit-testing arithmetic). Text edits, doc updates, single-function
changes behind a regression check, and first-pass audits go to GLM.

## Sources

- [RULER effective length, survey table (arXiv 2503.17407)](https://arxiv.org/pdf/2503.17407)
- [Context Windows Are Not Context Capacity (NeurIPS 2026)](https://neurips.cc/virtual/2026/173749)
- [Why 1M-token context is not reliable](https://memx.app/blog/why-1m-token-context-is-not-reliable/)
- [RepoAudit (arXiv 2501.18160, ICML 2025)](https://arxiv.org/html/2501.18160v3), code: https://github.com/PurCL/RepoAudit
- [SWE-agent ACI (arXiv 2405.15793)](https://hyper.ai/papers/2405.15793)
- [LLM×MapReduce (arXiv 2410.09342, ACL 2025)](https://arxiv.org/html/2410.09342v1)
- [Reducing false positives in static bug detection with LLMs (ICSE-SEIP 2026)](https://arxiv.org/html/2601.18844v1)
- [Utilizing precise code context for false-positive mitigation (arXiv 2411.03079)](https://arxiv.org/pdf/2411.03079)
- [Are LLMs reliable code reviewers? Overcorrection (arXiv 2603.00539)](https://arxiv.org/html/2603.00539v1)
- [Orchestrator-worker pattern (AgentPatterns.ai)](https://agentpatterns.ai/patterns/multi-agent/orchestrator-worker/)
- [Aider token limits](https://aider.chat/docs/troubleshooting/token-limits.html), [Aider repo map](https://aider.chat/docs/repomap.html), [HN: aider cannot split large files](https://news.ycombinator.com/item?id=38842643)
- This harness's own incident record: `C:\work\glm-harness\ORCHESTRATION-MANUAL.md`, `GLM_PLAYBOOK.md`, `AIDER-PLAYBOOK.md`
