# Whole-repo audit with a cheap limited-context model

The scripts of the 2026-10-17 ViewAlign v255 audit (glm-5.3 through the local `glm.js` harness,
which is not included: it holds a private endpoint). Method, numbers and false-positive patterns:
[../llm-orchestration-limited-context.md](../llm-orchestration-limited-context.md).

| order | script | does |
|---|---|---|
| 1 | `make_viewalign_global_audit.py <frozen-worktree>` | one read-only task per ~450-line chunk (cut at top-level boundaries) of every product file, plus cross-cutting tasks; writes `audit/index.txt` (symbol -> file:line) into the worktree |
| 2 | `run_va_global.sh` | runs the tasks one at a time (one request per key), resumable, retries with back-off |
| 3 | `reduce_va_global.py <wt> reduce` | keeps a finding only if its quoted line exists in the code; dedupes |
| 4 | `reduce_va_global.py <wt> verify` + `run_va_verify.sh` | two independent re-checks per finding: refute it / trace it blind from the user steps |
| 5 | `reduce_va_global.py <wt> merge` | CONFIRMED only if both re-checks confirm; REFUTED if one refutes and none confirms; else UNCERTAIN |
| - | `chain_va_global.sh` | runs 3-5 automatically after 2 finishes |

Paths are those of the original machine (`/c/work/...`); the file list and contracts in step 1 are
ViewAlign's. Audit a frozen `git worktree add --detach` snapshot, never the tree you are editing.
Every finding still needs a reproducer (a check that fails) before it is called a bug.
