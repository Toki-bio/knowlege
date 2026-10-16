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
