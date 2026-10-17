"""Global ViewAlign audit, stage 1 (find): one glm.js task per ~450-line chunk of every product file,
split at top-level boundaries, plus cross-cutting tasks. Method: knowlege/software/llm-orchestration-limited-context.md.

    python make_viewalign_global_audit.py <frozen-worktree>
Writes tasks/va-global/*.md and, in the worktree, audit/index.txt (every top-level symbol -> file:line).
"""
import io, os, re, sys

WT = sys.argv[1]
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'tasks', 'va-global')
os.makedirs(OUT, exist_ok=True)
os.makedirs(os.path.join(WT, 'audit'), exist_ok=True)
for f in os.listdir(OUT):
    os.remove(os.path.join(OUT, f))

FILES = {
  'script.js': 'the main viewer: parsing, rendering (Full/Block DOM views, Canvas view), editing, selection, search, consensus, codon analysis, trimming, clustering UI, snapshots, exports',
  'macse-align.js': 'exact port of MACSE v2.07 codon-aware alignment (must give identical output to MACSE)',
  'block-bicluster.js': '2D bicluster block detection over the alignment',
  'index.html': 'the page: controls and their inline handlers, script loading',
  'server.js': 'optional local Node server (static files, BLAST/MAFFT endpoints); must not expose files outside its root or run injected commands',
  'cluster.js': 'sequence clustering',
  'tree-draw.js': 'tree drawing',
  'blast-worker.js': 'web worker: local BLAST-like search',
  'macse-dp-wasm.js': 'WebAssembly dynamic-programming kernel generator for the MACSE port',
  'block-mask.js': 'Show 2D block mask computation',
  'realign-region.js': 'realign a selected region',
  'bam-parser.js': 'BAM file parser (reads view)',
  'codon-align.js': 'fast codon alignment engine',
  'kmer-tree.js': 'k-mer distance tree / grouping',
  'chromatogram-viewer.js': 'AB1 chromatogram viewer',
  'a11y.js': 'accessibility helpers',
  'peel.js': 'peel-by-diagnostic-columns grouping',
  'macse-worker.js': 'web worker running the MACSE port',
  'mafft-wasm.js': 'MAFFT WebAssembly wrapper',
  'ab1-parser.js': 'AB1 trace file parser',
  'doter-worker.js': 'web worker: dot plot',
  'doter-word-worker.js': 'web worker: word dot plot',
  'mafft-worker.js': 'web worker running MAFFT',
  'codon-align-worker.js': 'web worker for codon alignment',
}

CHECKLIST = """Look for these failure classes (each needs a concrete trigger):
1. Crash: null/undefined access, empty alignment, one sequence, zero-length or gap-only rows/columns,
   missing DOM element, NaN from parseInt/parseFloat of user input.
2. Wrong result: off-by-one (0- vs 1-based positions, inclusive vs exclusive ends), wrong row after
   reorder/delete/undo, gaps ('-' and '.') or lowercase or IUPAC codes or protein letters mishandled,
   wrong strand/frame, integer vs float, sorting instability changing output.
3. Stale state: a cache, index, selection, or result not invalidated after an edit, undo/redo, new file
   load, mode switch or setting change; an async result (worker, idle callback, timeout, fetch) applied
   after the data it was computed for has changed.
4. Undo/redo asymmetry: an action that cannot be undone, or whose undo restores the wrong thing.
5. Leaks/duplication: event listeners or timers added on every render/open and never removed; workers
   never terminated; DOM nodes accumulating.
6. Security: text from a loaded file (names, headers, annotations, BAM/AB1 fields, BLAST results) or a
   URL put into innerHTML/insertAdjacentHTML/attributes without escaping (XSS); in server.js path
   traversal, command injection, unbounded request bodies.
7. Data loss: export/save writing different data than shown (wrong rows, missing columns, wrong line
   wrapping, lost names), a parse that silently drops records or characters.
8. A user-visible text (message, tooltip, button) that states something the code does not do.
Not wanted: style, naming, refactoring ideas, performance (unless it freezes the page for seconds on a
normal 600 x 600 alignment), anything without a concrete trigger."""

COMMON = """You audit ViewAlign, a browser multiple-alignment viewer. Working directory = the repository.
File {file} ({purpose}). NEVER read a file whole: use read_file with offset/limit (at most 150 lines per
call) and grep. audit/index.txt lists every top-level function/const of the product files with file:line;
grep it to find a callee or caller instead of guessing what it does.

## Your chunk: {file} lines {a}-{b}
Read the whole chunk (in windows of 150 lines). Follow calls into other code only as far as needed to
decide whether a suspected defect is real (read the callee's lines, not its whole file).

{checklist}

## Output
For each real defect call append_items with list "findings" and objects with exactly these fields:
  "title": one line;
  "file": "{file}" or the other file the defect is in;
  "line": line number;
  "quote": the exact text of that line, character for character (it will be grep-checked);
  "class": the checklist number 1-8;
  "scenario": the concrete trigger (input or user steps) and the wrong result;
  "severity": "high" (crash, wrong data, data loss, security), "medium" (wrong picture, stale display), "low";
  "confidence": 1-5 (5 = traced through the code, certain).
At most 8 findings, the most serious first. If nothing in the chunk is a real defect, report nothing.
Then call finish.
"""

BOUND = re.compile(r'^(?:async\s+)?function\s|^(?:const|let|var|class)\s|^\}\s*;?\s*$|^// ={3,}|^/\*\*|^<script|^<!--|^\s{0,4}<(?:div|section|details|dialog|fieldset)\b')

index_lines = []
tasks = []
for fname, purpose in FILES.items():
    path = os.path.join(WT, fname)
    if not os.path.exists(path):
        print('missing', fname); continue
    lines = io.open(path, encoding='utf-8', errors='replace').read().split('\n')
    n = len(lines)
    for i, l in enumerate(lines, 1):
        m = re.match(r'^(?:async\s+)?function\s+([\w$]+)|^(?:const|let|var|class)\s+([\w$]+)', l)
        if m: index_lines.append(f'{fname}:{i}: {(m.group(1) or m.group(2))}')
    target = 450
    starts = [1]
    while starts[-1] + target < n:
        want = starts[-1] + target
        # nearest boundary after want-120 .. want+120 (a line that starts a top-level unit)
        cand = [k for k in range(want - 120, min(n, want + 120)) if BOUND.match(lines[k - 1] if k >= 1 else '')]
        nxt = min(cand, key=lambda k: abs(k - want)) if cand else want
        if nxt <= starts[-1] + 50: nxt = want
        starts.append(nxt)
    for ci, a in enumerate(starts):
        b = (starts[ci + 1] - 1) if ci + 1 < len(starts) else n
        tasks.append((fname, purpose, a, b))

CROSS = [
  ('xss-script', """Cross-cutting: every place in script.js that writes HTML (innerHTML, outerHTML, insertAdjacentHTML,
document.write, template strings put into those) - grep for them. For each, decide whether any interpolated
value can come from a loaded file or the URL (sequence names/headers, annotation names/descriptions, BAM
or AB1 fields, BLAST hit titles, snapshot JSON, tree labels) and whether it is escaped (_escapeHtml or
equivalent; grep its definition and check it escapes & < > " '). Report each unescaped path with the line
and a concrete file content that would inject markup (class 6). Work through the grep hits in order; at
most 8 findings."""),
  ('load-reset', """Cross-cutting: loading a new file must reset every per-alignment state. grep for where a new
alignment is installed (parseAndRender and the functions it calls) and list every state.* field and
module-level cache (let _... = ...) that holds per-alignment data (selections, caches, search history,
trims, groups, cluster results, block masks, codon data, annotation layout, undo history, Canvas state,
var-sites ranges, snapshots in memory). Report every one that survives a new load and would then show or
act on the previous alignment's data (class 3), with the concrete scenario."""),
  ('workers', """Cross-cutting: every Worker created in script.js (grep 'new Worker') and its message protocol with
the worker files (blast-worker.js, macse-worker.js, mafft-worker.js, doter-worker.js, doter-word-worker.js,
codon-align-worker.js). Check: results from a superseded run applied after a new file or a new run started
(class 3); workers never terminated or created on every run (class 5); error messages from the worker
not handled so the UI stays busy forever (class 1); message fields mismatched between sender and worker."""),
  ('undo', """Cross-cutting: the undo/redo system (grep deletedHistory, redoHistory, pushUndo, undo, redo).
For every action that pushes an undo entry, check that undo restores exactly the previous sequences,
names, order, selections it claims to, and that redo re-applies it; list actions that change sequences
or order but push no undo entry (class 4)."""),
  ('server', """Cross-cutting: server.js as a whole (read it in 150-line windows). It may run on a lab machine:
check path traversal and symlink escapes for static files, command/argument injection into BLAST/MAFFT
child processes, request size limits, CORS, temporary files never deleted, errors that crash the
process (class 6, 1)."""),
]

made = 0
for fname, purpose, a, b in tasks:
    key = f"{re.sub(r'[^a-z0-9]+', '-', fname.lower())}-{a:05d}"
    body = COMMON.format(file=fname, purpose=purpose, a=a, b=b, checklist=CHECKLIST)
    io.open(os.path.join(OUT, f'g-{key}.md'), 'w', encoding='utf-8', newline='\n').write(body)
    made += 1
for key, text in CROSS:
    body = COMMON.split('## Your chunk')[0].format(file='(several)', purpose='cross-cutting task') + '## Your task\n' + text + '\n\n' + CHECKLIST + '\n\n## Output' + COMMON.split('## Output')[1].format(file='(the file)')
    io.open(os.path.join(OUT, f'x-{key}.md'), 'w', encoding='utf-8', newline='\n').write(body)
    made += 1
io.open(os.path.join(WT, 'audit', 'index.txt'), 'w', encoding='utf-8', newline='\n').write('\n'.join(index_lines) + '\n')
print(f'{made} tasks ({len(tasks)} chunks + {len(CROSS)} cross-cutting), index {len(index_lines)} symbols')
