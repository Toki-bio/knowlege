#!/bin/bash
# Waits for stage 1 to finish, then reduces, generates two independent checks per finding and runs them, then merges.
cd /c/work/glm-harness || exit 1
export PYTHONIOENCODING=utf-8
until grep -q "^=== done" logs/va-global.log 2>/dev/null; do sleep 120; done
python reduce_va_global.py C:/work/MSA-viewer-glm-global reduce
python reduce_va_global.py C:/work/MSA-viewer-glm-global verify
bash run_va_verify.sh > logs/va-verify.log 2>&1
python reduce_va_global.py C:/work/MSA-viewer-glm-global merge
echo "=== chain done $(date +%H:%M:%S)"
