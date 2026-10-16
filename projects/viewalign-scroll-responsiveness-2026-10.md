---
title: ViewAlign scroll responsiveness at 600 x 600 (measured 2026-10-16)
type: notes
---

> Type: notes

What ViewAlign (C:\work\MSA-viewer-new, measured at 5d637fb / v1.33.2 and re-checked at 1a24de0 / v243 1.33.3, same results within noise) actually does when a
user scrolls a moderately sized alignment, measured on the owner's machine, and where
the time goes. Companion to the landscape and techniques write-ups in
[software/](../software/): [msa-viewer-rendering-landscape.md](../software/msa-viewer-rendering-landscape.md)
and [browser-grid-rendering-techniques.md](../software/browser-grid-rendering-techniques.md).

## Set-up

| item | value |
|---|---|
| machine | AMD Radeon 780M (D3D11 via ANGLE), 16 threads, Chrome 155, canvas + rasterization + compositing all "Hardware accelerated" in chrome://gpu (headless and headful) |
| viewport | 1600 x 900 window; alignment area about 1582 x 819 px |
| cell | 8 x 15 px at the default 13 px "Courier New"; about 197 x 54 = 10,600 cells visible |
| fixtures | `_perf_600x600.fas`: 600 rows of the oma SubFam plate (CURATE__oma__subfam600_23seeds), densest 600-column window (23 % residues, the rest gaps); `_perf_1200x1200.fas`: the same rows doubled, 1,200 columns |
| input | 60 wheel steps of 120 px at 16 ms cadence through CDP `Input.dispatchMouseEvent mouseWheel`, one run per direction |
| metrics | requestAnimationFrame deltas (p50, p95, max, frames > 34 ms), Long Animation Frames with script attribution, scroll distance actually reached, CPU profile and devtools trace of 20 steps |
| scripts | ViewAlign repo `dev-tools/scroll-perf-tests/measure_modes.js`, `trace_mode.js`, `inval_trace.js`, `canvas_draw_bench.html` + `run_canvas_draw_bench.js` (run with `NODE_PATH=<repo>/node_modules`) |

Headless Chrome gives the same ranking as headful on this machine; headful numbers are
quoted where both exist.

## Where 600 x 600 lands in the code

600 x 600 = 360,000 residues. That is below `ALIGN_WINDOWED_DOM_THRESHOLD` (500,000) and far
below `ALIGN_CRAZY_VOLUME` (5,000,000), so the default path is the classic full DOM build:
one `<span>` per residue (361,201 spans), `display: inline-block`, native document
scrolling, the span cache off (over 80,000 residues). Off-screen rows get
`content-visibility: auto` (`styles.css` rule on `.block-block.cv-block ... > .seq-data`,
commit 0b330db) so that the initial layout does not cost 2 s. Canvas mode exists but is
opt-in at this size; it auto-switches only at 5 M residues.

## Results, 600 x 600

Mode switch into the view: Full 2.4-2.8 s, Block 1.8-1.9 s, Canvas 0.18 s.

| mode | direction | frame p50 | frame p95 | frame max | frames > 34 ms (of ~150) | wall for 60 steps |
|---|---|---|---|---|---|---|
| Full (DOM) | vertical | 16.7 ms | 109-158 ms | 134-266 ms | 58-61 | 6.8-9.1 s |
| Full (DOM) | horizontal | 8.4 ms | 9-17 ms | 17 ms (one 208 ms spike) | 0-4 | 1.3 s |
| Block (DOM) | vertical | 83-100 ms | 150-167 ms | 175-217 ms | 95 | 9-11.5 s |
| Canvas | vertical | 8.4 ms | 33-42 ms | 42-100 ms | 6-13 | 2.3-2.7 s |
| Canvas | horizontal | 8.4 ms | 42 ms | 50-92 ms | 8-17 | 1.8-2.3 s |

A perfectly smooth run would take 1.0 s of wall time for 60 steps. Full mode vertical
scrolling is 7-9x slower than the input, Block mode 9-11x, Canvas 2-2.7x.

Vertical scrolling in Full and Block mode scrolls the **document**, not the alignment
container: the container is 7,859 px (Full) or 31,430 px (Block) tall and the page
scrolls as a whole. Horizontal scrolling is mirrored through the custom
`horizontal-scrollbar` element.

## Where the time goes

**DOM modes (Full, Block).** No script runs in the long frames (Long Animation Frames carry
no script attribution; the CPU profile is 90 % `(program)`). Per wheel step the trace shows
one `Layout` of 10,000-12,000 dirty objects out of 218,000-240,000 (about 12 ms), three
`Paint` events covering the whole 4,855 x 7,859 px container layer (about 31 ms each), and a
`Layerize` of about 14 ms. The invalidation trace names the cause: `.seq-data` rows report
`StyleRecalcInvalidationTracking ... DisplayLock` (the `content-visibility: auto` lock
opening as a row enters the viewport) and about 130,000 spans report "Added to layout"
over 12 steps. Scrolling is therefore an incremental layout + full-layer repaint of every
row that comes into view, roughly 8 rows x 600 spans per 120 px step. Disabling the sticky
name column (`position: static`) changes nothing.

**Canvas mode.** The draw is viewport-only and correct, but one frame costs about 58 ms
(FireAnimationFrame 1,153 ms over 20 steps). Self time per frame: `drawImage` 51 ms,
`_makeGlyph` 5 ms (string key concatenation + Map lookup per cell), the per-cell style
closure 4.7 ms, `draw` 5.4 ms. The glyph cache turns each cell into one `drawImage` from a
separate small `<canvas>` per (letter, background, foreground, variant), so a frame is
10,600 `drawImage` calls from hundreds of different source canvases. The GPU process shows
46 `GPUTask` entries totalling 1.3 s for the same 20 steps.

## Micro-benchmark: how to paint 10,638 cells per frame on this GPU

Standalone page (`canvas_draw_bench.html`), same cell size, random ACGT- data, 60 frames per
strategy scrolling 8 rows per frame. Numbers are JavaScript time of the draw call sequence;
the "to 2nd rAF" rows include presentation.

| strategy | headful p50 | headful p95 | headless p50 | headless p95 |
|---|---|---|---|---|
| A fillRect + fillText per cell | 67 ms | 203 ms | 51 ms | 271 ms |
| A2 rects batched by colour, then one text pass | 25 ms | 70 ms | 35 ms | 61 ms |
| B drawImage from per-glyph canvases (ViewAlign today) | 36 ms | 46 ms | 57 ms | 75 ms |
| C drawImage from one atlas canvas (9-arg) | 29 ms | 40 ms | 55 ms | 78 ms |
| C2 drawImage from one atlas ImageBitmap | 30 ms | 40 ms | 47 ms | 55 ms |
| D blit previous frame + draw the 8 exposed rows | 3.9 ms | 6.3 ms | 7.2 ms | 10 ms |
| E coloured rects only, no letters | 3.4 ms | 34 ms | 4.3 ms | 75 ms |
| F cached per-row strips (warm), 54 drawImage per frame | 7.0 ms | 10 ms | 11 ms | 14 ms |

Draw-to-second-rAF (headful): A 67 ms, B 42 ms, C 33 ms, D 16.7 ms, F 16.7 ms (16.7 ms is
one 60 Hz frame, the floor of this measurement).

Reading: on this machine a Canvas 2D `drawImage` costs 3-5 µs regardless of whether the
source is a tiny canvas, an atlas, or an ImageBitmap, so **any strategy that re-issues one
draw call per visible cell cannot fit a 16 ms frame at 10,600 cells**. The atlas saves
15-20 % over per-glyph canvases; it does not change the class of the problem. Only
drawing fewer things per frame does: blitting the unchanged part of the previous frame and
painting the newly exposed strip (D), or drawing cached row strips (F), both stay under
8 ms with no loss of exactness (every cell is still drawn from the same glyph bitmaps).
Rect-only drawing (E) has a fast median but large spikes, so it is not a free zoom-out path
either.

## What this means for ViewAlign (owner's call, nothing changed)

1. At 600 x 600 the default (DOM) view is the slow one; Canvas is already 4x smoother but
   is hidden behind a mode switch and still misses frames on this GPU.
2. The Canvas frame is bound by draw-call count, not by JavaScript logic. The two
   structural fixes that keep exactness are: (a) scroll by blitting the current canvas and
   drawing only the exposed rows/columns; (b) cache rendered row strips (or tiles) and
   compose the viewport from them, invalidating a strip when its row changes. Both are
   what desktop viewers and spreadsheet canvases do; see the techniques write-up.
3. Cheap wins inside the current draw: build the glyph key from small integers instead of
   five concatenated strings; resolve the cell style once per column (conservation) and
   once per base, not per cell; avoid two `ctx.font` switches per row for the names (draw
   all names in one pass after the residues).
4. The DOM path cannot be made smooth at this size by CSS alone: `content-visibility: auto`
   trades initial layout for per-scroll layout + repaint. Making Canvas the default
   above about 100,000 residues (and keeping DOM for editing-heavy small alignments) is
   the simplest route to "instant" at 600 x 600; the user-facing question is what Canvas
   mode still lacks (consensus row, hidden columns in var-sites mode, Show 2D overlay).

## Open points

- Measured on one machine (Radeon 780M). The per-call cost of `drawImage` differs per GPU
  and driver; the ranking of strategies should not, but the absolute budget does.
- The 600 x 600 fixture is 77 % gaps. A residue-dense alignment paints more distinct
  colours per frame; the DOM result would not improve, the Canvas glyph-cache misses
  would be slightly more frequent in the first frames.
- Headless and headful agreed in ranking but not in absolute numbers (headful faster on
  drawImage, slower on fillText).

## Done the same day (v244, commits 0abeb93 + 0adf9a6 on Toki-bio/MSA-viewer main)

Canvas mode now scrolls by copying the previous frame by the whole device-pixel delta
(through a cleared scratch canvas: an in-place `drawImage` blends antialiased pixels, which
the pixel-identity check caught) and drawing only the exposed strip; offsets are quantised to
device pixels; with sticky names a horizontal step leaves the name column in place; names are
drawn in one pass per frame. Full redraw on dirty state, diagonal moves and page jumps.

| Canvas mode, 600 x 600, headful | before (v243) | after (v244) |
|---|---|---|
| vertical frames p95 / max | 33-42 / 50-100 ms | 8.5-9 / 9-33 ms |
| horizontal frames p95 / max | 33-42 / 50-92 ms | 8.5-8.8 / 9 ms |
| frames over 34 ms (of ~160-190) | 5-17 | 0 |
| cost of one scroll frame (regression check) | 25-50 ms (full) | 1.3-3.9 ms (copy) |

Regression check "Canvas mode: wheel scroll copies the frame ..." in `tests/regression/run-all.js`
asserts the copy path is used, that pixels equal a full redraw at the same offset, and that a
copy frame is well under a full frame. Still open: Full/Block (DOM) modes at this size are
unchanged (the default view), so the user only gets the fast path after choosing Canvas; the
decision whether Canvas becomes the default above ~100k residues (it lacks the consensus
row, hidden columns in var-sites mode and the Show 2D overlay) is the owner's.

## Can the DOM view be made cheaper? (tested the same day)

The user asked for this to be tested before deciding. Full mode, 600 x 600, headful, same
wheel protocol; the only change is CSS injected after the render.

| variant | vertical frames p95 | frames > 34 ms | wall for 60 steps |
|---|---|---|---|
| as shipped (`content-visibility: auto` on rows) | 110-160 ms | 58-61 | 6.8-9.1 s |
| plus `contain: layout style paint` on `.seq-data` | 117-125 ms | 58 | 6.7-7.0 s |
| `content-visibility` rule off (whole page laid out once) | 17-25 ms | 0-2 | 1.3-1.4 s |

Turning the rule off makes vertical scrolling smooth (the page is laid out once and the
compositor scrolls it). The price is paid by every operation that lays the page out again:

| Full mode, operation | 600 x 600, rule on | 600 x 600, rule off | 700 x 700, rule on | 700 x 700, rule off |
|---|---|---|---|---|
| switch to Full | 2.6 s | 5.8 s | 3.2 s | 10.1 s |
| Highlight diffs on | 2.7 s | 5.5 s | 3.2 s | 7.9 s |
| switch to Block | 1.9 s | 4.6 s | 2.5 s | 8.7 s |
| zoom back to 100 % | 0.4 s | 0.01 s | 0.5 s | 0.2 s |

So the DOM view at this size can scroll smoothly or switch/redraw in 2-3 s, not both; with
the rule off a horizontal scroll also shows one 0.7 s hitch when the wide layer is first
painted. Neither variant approaches the Canvas numbers, and the DOM operations are already
slow in absolute terms (a 361k-span layout costs seconds). Conclusion recorded for the owner:
the DOM view stays as it is (editing and overlays need it below 100k residues), and Canvas
becomes the default above `ALIGN_CANVAS_DEFAULT_RESIDUES` = 100,000 (v245) with the consensus
row drawn in Canvas. A genuinely cheaper DOM would need fewer elements (one span per run of
equal colour, or windowed rows), which is a rewrite of selection and editing, not a CSS change.

Scripts: `dev-tools/scroll-perf-tests/measure_modes.js` with `INJECT_CSS=...`, and the
`op_cost_content_visibility.js` next to it.

## v246-v247: editing in Canvas, hardening

v246: editing in Canvas mode (row drag by the name, rename by double-click, the toolbar no
longer pushes the canvas off screen, static-name hit test, chip no longer "view-only"); the
tools themselves (typing, gaps, Move, Slide, column select) already worked and a guard check
now asserts they give the same sequences as in Full mode. v247: device-pixel cell geometry,
strip snapped to row/column edges, a state-signature safety net for the frame copy, live
conservation/consensus caches, PageUp/PageDown in Canvas; the copy is pixel-identical at 1x,
1.25x and 1.5x scaling and after edits (see techniques write-up, section 10).

## Beyond scrolling: where the time went after v247, and v248 (2026-10-16)

Profiled in the default (Canvas) view at 600 x 600, headful, real input where possible.

| action | v247 | v248 |
|---|---|---|
| load to first Canvas frame | 0.2 s (the 2.5 s seen earlier was page navigation + fetch) | same |
| click a residue / Ctrl+Alt+click a column / Ctrl+click a name | one full redraw each, 29-45 ms | no cell redraw, overlay only |
| drag a residue-range selection | one full redraw per mouse move | overlay only |
| type one residue | one full redraw, 27-29 ms | one row + one column band repainted |
| undo, gap insert, gap column | full re-render (model operation), 30-35 ms | same |
| zoom, consensus toggle, Highlight diffs, search | full re-render, 15-36 ms wall | same |
| switch to Full (DOM) | 2.5 s (361k spans) | same |

v248: a transparent **overlay canvas** above the content holds selections, the pending and
edit cell, and the row-drag indicator; selection paths request an overlay redraw only, and
the content signature no longer includes selection state. Typing calls `redrawCells(row, col)`,
which repaints the row band, the column band and the consensus band and then accepts the new
state as drawn, so the next scroll step can still copy the frame. Check: four selection
actions and seven keystrokes cause zero full redraws and the view equals a full redraw.

DOM views below the Canvas default (Full mode, same protocol): 150 x 150 (22k residues) smooth
(p95 8.8 ms); 200 x 200 (40k) p95 16.5 ms with 2-3 frames of 50 ms; 300 x 300 (90k) p95 67 ms,
16 frames over 34 ms. So the DOM view is only smooth up to roughly 30-40k residues; the
default threshold of 100k leaves 40k-100k alignments on the janky DOM path unless the user
picks Canvas. Open for the owner: lower `ALIGN_CANVAS_DEFAULT_RESIDUES` to 40,000, at the cost
of Show 2D and var-sites column hiding needing a manual switch for those sizes.

## v249 (2026-10-16): Canvas default lowered to 40,000 residues

Owner's decision after the DOM measurements above ("yes, small indicator should show the
limitation statement for like 10 seconds"). `ALIGN_CANVAS_DEFAULT_RESIDUES` is 40,000 (was
100,000). When the view is chosen automatically, a 10-second status note names what Canvas
does not do yet, and `showMessage` appends that note to the next message so the following
"File loaded successfully!" does not hide it. Commits 200c5ee + f002692 on Toki-bio/MSA-viewer
main; quick suite: only the two pre-existing failures (Highlight diffs in-place check in DOM
modes, optional server-security suite on this Windows account).

### Could Canvas serve every alignment?

Nothing in Canvas is impossible; what remains is DOM-only today because it was written
against DOM spans and has not been ported:

| feature | Canvas today |
|---|---|
| Show 2D (block mask overlay, `renderBlockMaskOverlay`) | **ported in v254** (2026-10-17): overlay layer, `_canvasBlockMaskRects` + `_canvasBlockMaskHit` for click/tooltip |
| var-sites mode: hiding invariant columns | **ported in v251** (2026-10-17): slot layout (shown columns + marker cells) used by drawing, overlay, hit-testing, tooltips and jumps |
| codon mode: amino-acid annotation rows below each sequence | **ported in v253** (2026-10-17): stop cells, syn/nonsyn underlines, codon boxes, frameshift marks, frame tags, AA only; colours read from the stylesheet via probe elements |
| snapshots | **v254**: record and restore Canvas |
| Block (wrapped) layout | Canvas is a single unwrapped grid |
| Highlight diffs (dim conserved columns to 40 %) | **ported in v250** (2026-10-16): per-cell `globalAlpha` 0.4 and bold glyph variant in `drawRegion`, consensus and ruler included; frame copy unaffected |

Each is a drawing routine on the same data; the hard part is the cell-count arithmetic in
hit-testing once columns are hidden or rows carry extra annotation lines. Porting them one at
a time would make Canvas the only view and remove the mode switch.

### Porting order (owner: "ok go", 2026-10-16)

1. Highlight diffs: v250, done. Lesson for checks on the Canvas: untouched canvas pixels are
   transparent black, so a pixel test must composite over white before judging darkness, and
   "bold" must be measured as letter pixels against the cell's own background (white-on-black
   cells have fewer dark pixels when bold).
2. Variable sites only: v251, done. Design: a `slots` Int32Array (column, or `-1 - column` for the
   marker before that column) and `slotOf[column]` (slot or -1) built once per render; every
   x = slot * cellW, so the uniform grid, the frame copy and the strip snapping are untouched.
   Found on the way: `generateScale()` drops a label that does not fit the window it is asked
   for, so a ruler text made per drawn strip differs from a full redraw (now one text for the
   whole alignment); the rule under the ruler was drawn beneath sticky names when scrolled.
   Harness lessons: at a 700 px window the FASTA textarea overlaps the canvas (wheel went
   there); the zoom slider's render after a tooltip hover was lost once, so checks that need
   scrollable content load a larger alignment instead of zooming.
3. Codon annotation rows, snapshot restore, Show 2D, Block layout: after that, one at a time.

### v252-v253 (2026-10-17): owner reports and the codon port

- **Name clicks selected a residue** (v252): with sticky names and a scrolled view, the Canvas
  hit test returned the cell under the name column; Ctrl+C then copied that residue. Fixed in
  the three hit tests (name column = names).
- **"Low resolution" Canvas text** (v252): Chrome antialiases text on a transparent canvas in
  greyscale, on an opaque one with subpixel ClearType like DOM text. Measured headful on
  Windows (colour-fringe pixels in a 400 x 150 crop): Canvas 0, Full 10,375; opaque Canvas
  9,704. Content canvas, scroll scratch and glyph canvases are opaque now; clears are white
  fills; translucent tints compose over white. Headless Chrome has no subpixel text, so the
  check verifies the contract (context attributes, no transparent pixels), not the fringes.
- A screen-scale change (devicePixelRatio) now re-renders.
- Codon port (v253): rounded box ends must be image blits, not path fills, or the frame copy
  differs by one unit at the antialiased corners under the strip clip.
- GLM use: aider/glm-5.3 cannot take script.js (416k tokens > 262k context); it did the
  manual/index text edits in a worktree (reviewed, cherry-picked); glm.js (read-only) audited
  the codon block against the DOM functions and confirmed it. Windows: run aider with
  PYTHONUTF8=1 or it crashes on the rich progress bar (cp1252).

Remaining DOM-only after v254: only the Block (wrapped) layout, by decision (see the Plan below).

## Plan (owner's decision 2026-10-17): make the DOM Block view faster, do not port wrapping to Canvas

Why not Canvas: wrapping is for few rows over many columns, read like a page; at 600 x 600 a
wrapped view is six 600-row blocks and nobody reads that. The only gap is a wrapped view of an
alignment above ~40k residues where DOM Block is too slow, which is rare, and the port would be
as large as v250-v254 together (every Canvas geometry becomes block + column-in-block, plus a
new round of frame-copy exactness).

Why Block is slow (measured 2026-10-16): each scroll step unlocks `content-visibility: auto`
rows, a layout pass plus a whole-layer repaint; not JavaScript. Block pays more than Full
because each block is its own stack of rows.

Steps, in order:
1. Measure with `dev-tools/scroll-perf-tests/measure_modes.js` in Block mode at 10 x 5,000,
   20 x 5,000 and 50 x 3,000 (50k, 100k, 150k residues): baseline; then with
   content-visibility off (the harness has INJECT_CSS); then with the windowed renderer
   forced on (`ALIGN_WINDOWED_DOM_THRESHOLD` lowered). Record p50/p95/max frame times and the
   wall time of a mode switch and a Highlight diffs toggle for each.
2. If content-visibility off gives p95 under ~25 ms with rebuilds under ~3 s up to ~120k
   residues: ship it as a size band (off between the Canvas default and that limit), keep it on
   above. One threshold constant, a regression check that scroll frames stay under the budget
   at 100k in Block, and that the rebuild time stays bounded.
3. If the windowed renderer is the better trade from ~120k up: lower its threshold for Block
   only, after checking the v230 windowed checks still pass.
4. Expected result: Block p95 from ~70-100 ms to ~20-30 ms in that band; it stays above the
   Canvas numbers (8-9 ms, flat), which is the ceiling for any DOM path.
Tested and not worth repeating: paint containment, will-change (no gain).
