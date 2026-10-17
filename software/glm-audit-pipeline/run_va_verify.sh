#!/bin/bash
# Global ViewAlign audit, stage 2 (skeptical verifier per finding): every task in tasks/va-global, sequential (one request per key),
# resumable (skips tasks with a good out/*.json), retried with back-off on 429/network.
cd /c/work/MSA-viewer-glm-global || exit 1
export PYTHONIOENCODING=utf-8
for f in /c/work/glm-harness/tasks/va-verify/*.md; do
  t=$(basename "$f" .md)
  for try in 1 2 3 4 5 6; do
    if [ -s "/c/work/glm-harness/out/$t.json" ] && ! grep -q "\"error\"" "/c/work/glm-harness/out/$t.json"; then break; fi
    echo "=== $t try $try $(date +%H:%M:%S)"
    GLM_READ_ROOTS="C:/work/MSA-viewer-glm-global" timeout 1500 node /c/work/glm-harness/glm.js "$f" 2>&1 | tail -1
    if [ -s "/c/work/glm-harness/out/$t.json" ] && ! grep -q "\"error\"" "/c/work/glm-harness/out/$t.json"; then break; fi
    sleep $((60 * try))
  done
done
echo "=== done $(date +%H:%M:%S)"
