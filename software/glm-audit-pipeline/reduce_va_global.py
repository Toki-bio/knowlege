"""Global ViewAlign audit, stage 1 reduce + stage 2 task generation.

    python reduce_va_global.py <frozen-worktree> reduce   -> out/va-global-reduced.json, summary
    python reduce_va_global.py <frozen-worktree> verify   -> tasks/va-verify/v-NNN.md (one per kept finding)
    python reduce_va_global.py <frozen-worktree> merge    -> out/va-global-final.json (findings + verdicts)

Reduce: every quote must exist in the named file at/near the stated line (else anywhere in that file,
else in any audited file); fabricated quotes are dropped; duplicates (same file, line within 3, similar
title) merged, keeping the highest severity/confidence and counting how many tasks reported it.
"""
import glob, io, json, os, re, sys

WT, MODE = sys.argv[1], sys.argv[2]
HERE = os.path.dirname(os.path.abspath(__file__))
OUTD = os.path.join(HERE, 'out')
norm = lambda s: ' '.join((s or '').split())
_files = {}
def file_lines(f):
    if f not in _files:
        p = os.path.join(WT, f)
        _files[f] = io.open(p, encoding='utf-8', errors='replace').read().split('\n') if os.path.isfile(p) else None
    return _files[f]

def locate(it):
    q = norm(it.get('quote'))
    if not q: return None
    try: ln = int(it.get('line') or 0)
    except Exception: ln = 0
    names = [it.get('file')] if it.get('file') and file_lines(it.get('file')) else []
    names += [f for f in ('script.js', 'index.html', 'server.js') if f not in names]
    for f in names:
        L = file_lines(f)
        if not L: continue
        near = [k for k in range(max(1, ln - 6), min(len(L), ln + 6) + 1) if norm(L[k - 1]) == q or (len(q) > 25 and q in norm(L[k - 1]))]
        if near: return f, min(near, key=lambda k: abs(k - ln))
        hits = [k for k, l in enumerate(L, 1) if norm(l) == q or (len(q) > 25 and q in norm(l))]
        if hits: return f, min(hits, key=lambda k: abs(k - ln))
    return None

sev = {'high': 0, 'medium': 1, 'low': 2}
if MODE == 'reduce':
    kept, dropped, errors, empty = [], [], [], 0
    for f in sorted(glob.glob(os.path.join(OUTD, '[gx]-*.json'))):
        if f.endswith('.transcript.json'): continue
        task = os.path.basename(f)[:-5]
        try: d = json.load(io.open(f, encoding='utf-8'))
        except Exception as e: errors.append((task, str(e))); continue
        if d.get('error'): errors.append((task, d.get('error'))); continue
        items = d.get('findings') or []
        if not items: empty += 1
        for it in items:
            it['task'] = task
            loc = locate(it)
            if not loc: it['why_dropped'] = 'quote not found'; dropped.append(it); continue
            it['file'], it['line_checked'] = loc
            kept.append(it)
    kept.sort(key=lambda x: (sev.get(x.get('severity'), 3), -int(x.get('confidence') or 0)))
    merged = []
    for it in kept:
        words = set(re.findall(r'[a-z]{4,}', norm(it.get('title')).lower()))
        dup = next((m for m in merged if m['file'] == it['file'] and abs(m['line_checked'] - it['line_checked']) <= 3
                    and len(words & m['_w']) >= max(2, len(words) // 3)), None)
        if dup: dup['reports'] += 1; dup['tasks'].append(it['task']); continue
        it['_w'] = words; it['reports'] = 1; it['tasks'] = [it['task']]; merged.append(it)
    for m in merged: m.pop('_w', None)
    for n, m in enumerate(merged, 1): m['id'] = n
    json.dump({'kept': merged, 'dropped': dropped, 'errors': errors}, io.open(os.path.join(OUTD, 'va-global-reduced.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    by = {}
    for m in merged: by[m.get('severity')] = by.get(m.get('severity'), 0) + 1
    print(f'tasks with output: {len(glob.glob(os.path.join(OUTD, "[gx]-*.json"))) - len(glob.glob(os.path.join(OUTD, "[gx]-*.transcript.json")))}, empty {empty}, errors {len(errors)}')
    print(f'findings kept {len(merged)} {by}, dropped (quote not found) {len(dropped)}')

elif MODE == 'verify':
    d = json.load(io.open(os.path.join(OUTD, 'va-global-reduced.json'), encoding='utf-8'))
    VD = os.path.join(HERE, 'tasks', 'va-verify'); os.makedirs(VD, exist_ok=True)
    for f in os.listdir(VD): os.remove(os.path.join(VD, f))
    for m in d['kept']:
        a = max(1, m['line_checked'] - 60)
        body = f"""You are a skeptical reviewer. Another reviewer claims a defect in ViewAlign (browser alignment viewer).
Working directory = the repository. Never read a file whole: read_file with offset/limit (max 150 lines), grep.
audit/index.txt lists every top-level function with file:line.

## The claim
Title: {m.get('title')}
File: {m['file']}, line {m['line_checked']}: {norm(m.get('quote'))}
Scenario given: {m.get('scenario')}

## Your job
Try to REFUTE it. Read {m['file']} lines {a}-{m['line_checked'] + 60}, then the callers (grep the function name) and
anything that guards the path (checks earlier in the call chain, UI that makes the state impossible, a reset
elsewhere). Decide:
- "CONFIRMED" only if you traced the scenario step by step and the wrong result really happens;
- "REFUTED" if a guard, a reset, a caller or the UI prevents it, or the claim misreads the code (quote the line);
- "UNCERTAIN" if it depends on something you cannot determine from the code.
Call append_items once, list "verdict", one object: {{"id": {m['id']}, "verdict": "...", "evidence": "the decisive lines (file:line + quote) and why", "repro": "for CONFIRMED: exact steps or a JavaScript snippet for the browser console that shows the wrong result; else empty"}}.
Then call finish.
"""
        io.open(os.path.join(VD, f"v-{m['id']:03d}-a.md"), 'w', encoding='utf-8', newline='\n').write(body)
        # Second, independent check: no claim title, no reasoning - only the place and the user steps.
        trace = f"""You check one suspected problem in ViewAlign (browser alignment viewer) by tracing it yourself.
Working directory = the repository. Never read a file whole: read_file with offset/limit (max 150 lines), grep.
audit/index.txt lists every top-level function with file:line.

## What to trace
Code location: {m['file']} line {m['line_checked']}: {norm(m.get('quote'))}
User steps or input: {m.get('scenario')}

## Your job
Start from the user's action (find the event handler or entry point with grep), follow the real call path
to the code location and beyond, and write down what the program actually does at each step with these
inputs. Then say whether the outcome is wrong (crash, wrong data, wrong picture, stale state, data loss,
unescaped markup from file content), and what exactly is wrong. Do not assume the steps lead to a
problem; many such suspicions are wrong because of a check or reset elsewhere.
Call append_items once, list "verdict", one object: {{"id": {m['id']}, "verdict": "CONFIRMED" | "REFUTED" | "UNCERTAIN", "evidence": "the call path you traced with file:line for each step, and the decisive line", "repro": "for CONFIRMED: exact steps or a browser-console JavaScript snippet that shows it; else empty"}}.
Then call finish.
"""
        io.open(os.path.join(VD, f"v-{m['id']:03d}-b.md"), 'w', encoding='utf-8', newline='\n').write(trace)
    print(len(d['kept']), 'findings,', 2 * len(d['kept']), 'verify tasks')

elif MODE == 'merge':
    d = json.load(io.open(os.path.join(OUTD, 'va-global-reduced.json'), encoding='utf-8'))
    v = {}
    for f in glob.glob(os.path.join(OUTD, 'v-*.json')):
        if f.endswith('.transcript.json'): continue
        side = 'a' if f.endswith('-a.json') else 'b'
        try:
            for it in (json.load(io.open(f, encoding='utf-8')).get('verdict') or []): v[(int(it.get('id')), side)] = it
        except Exception: pass
    for m in d['kept']:
        a, b = v.get((m['id'], 'a')) or {}, v.get((m['id'], 'b')) or {}
        va, vb = (a.get('verdict') or 'MISSING').upper(), (b.get('verdict') or 'MISSING').upper()
        # both independent checks must confirm; a refutation with no confirmation = refuted; else uncertain
        if va == 'CONFIRMED' and vb == 'CONFIRMED': m['verdict'] = 'CONFIRMED'
        elif 'REFUTED' in (va, vb) and 'CONFIRMED' not in (va, vb): m['verdict'] = 'REFUTED'
        else: m['verdict'] = 'UNCERTAIN'
        m['verdicts'] = [va, vb]
        m['v_evidence'] = [a.get('evidence'), b.get('evidence')]
        m['repro'] = a.get('repro') or b.get('repro')
    json.dump(d, io.open(os.path.join(OUTD, 'va-global-final.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    tally = {}
    for m in d['kept']: k = (m.get('severity'), m['verdict']); tally[k] = tally.get(k, 0) + 1
    for k in sorted(tally, key=lambda k: (sev.get(k[0], 3), k[1])): print(k, tally[k])
