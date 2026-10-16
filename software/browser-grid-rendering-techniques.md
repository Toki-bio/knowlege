---
title: Lag-free rendering of a large coloured character grid in the browser
type: explanation
---

> Type: explanation

How to scroll a 600 x 600 to 5,000 x 5,000 grid of coloured letters (an alignment, a
terminal, a spreadsheet) in a browser with no visible lag and every cell drawn exactly.
Written 2026-10-16 from primary sources (browser vendor docs, specs, Chromium source, and the
source of xterm.js, Glide Data Grid, VS Code, igv.js, JBrowse 2) plus the viewer survey in
[msa-viewer-rendering-landscape.md](msa-viewer-rendering-landscape.md) and the measurements
in [projects/viewalign-scroll-responsiveness-2026-10.md](../projects/viewalign-scroll-responsiveness-2026-10.md).
Anything marked "unverified" was not confirmed from a primary source.

## 1. The three ways the problem has been solved

**DOM per cell** (one element per residue). Layout cost scales with the number of cells, not
the viewport. VS Code's terminal left it in 2017 because layout alone could exceed a 16.6 ms
frame and DOM churn caused GC pauses; the canvas renderer was "5 to 45 times faster"
(https://code.visualstudio.com/blogs/2017/10/03/terminal-renderer). AlignmentViewer 2.0's
README: "Inserting many elements into the HTML DOM at one time absolutely kills
performance". Virtualising the DOM (Monaco's `viewLayer.ts`, MSABrowser) keeps only viewport
rows and works for text lines, but every scroll step still builds and lays out new nodes.
`content-visibility: auto` on rows defers the cost but does not remove it: the rows entering
the viewport are laid out and painted during the scroll (measured in ViewAlign: 100-160 ms
frames at 600 x 600).

**Canvas 2D, viewport-sized.** One canvas the size of the view, redrawn from data on each
scroll. Everyone in the survey who is fast uses some form of this. The cost per frame is the
number of draw calls, and that is where implementations differ (section 3).

**WebGL / WebGPU instanced quads.** xterm.js's WebGL renderer draws text as one instanced
quad per cell (`drawElementsInstanced`, 11 floats per cell, two alternating Float32Arrays
uploaded with `STREAM_DRAW`) from a glyph atlas of up to 32 pages of at most 4,096 px; VS
Code 1.41 reports canvas frames of "1 ms to 40 ms" vs WebGL "0.2 ms to 2 ms"
(https://code.visualstudio.com/updates/v1_41). Costs: context-loss handling (xterm has a 3 s
restore timer), about sixteen WebGL contexts per page before trouble (JBrowse 2 design notes),
and a shader pipeline to maintain. Figma rejected Canvas 2D for being immediate-mode with no
guaranteed GPU acceleration (https://www.figma.com/blog/building-a-professional-design-tool-on-the-web/).

## 2. What a frame may cost

60 Hz leaves 16.7 ms per frame, of which about 10 ms is yours after browser overhead
(https://web.dev/articles/rendering-performance); 120 Hz displays leave 8.3 ms. On the
ViewAlign test machine (Radeon 780M, Chrome 155, canvas hardware-accelerated) one Canvas 2D
`drawImage` of an 8 x 15 px cell costs 3-5 µs whatever the source (small canvas, atlas,
ImageBitmap), and `fillRect` + `fillText` per cell costs 5-7 µs with spikes to 20-50 µs. A
1,600 x 900 viewport holds about 10,600 cells, so **any design that issues one draw call per
visible cell per frame costs 30-60 ms on that machine** and misses every frame budget. The
viewers that feel smooth never do this; they draw fewer things per frame.

## 3. Drawing fewer things per frame

Ranked by measured or documented effect; the first two are the ones that change the class of
the problem.

1. **Blit the previous frame, paint only the exposed strip.** `drawImage(canvas, ...)` onto
   itself (the spec says the source is copied first, so overlap is legal:
   https://html.spec.whatwg.org/multipage/canvas.html), then draw the newly exposed rows or
   columns. Jalview's `SeqCanvas.fastPaint` does this with `copyArea` on its back image; Glide
   Data Grid's `blitLastFrame` does it on canvas with integer device-pixel offsets,
   `imageSmoothingEnabled = false` during the copy, single-axis only (diagonal scroll, theme
   or selection change, or more than 100 columns changing force a full redraw)
   (https://raw.githubusercontent.com/glideapps/glide-data-grid/main/packages/core/src/internal/data-grid/render/data-grid-render.blit.ts).
   VS Code's minimap does the same on `ImageData` with a typed-array row copy and a dirty-band
   `putImageData`. Measured in the ViewAlign micro-benchmark: 3.9 ms p50 for a 120 px step
   versus 29-36 ms for a full redraw from an atlas. A page-sized jump falls back to a full
   redraw (Jalview: "scrolling too fast, more than a page at a time").
2. **Cache rendered strips or tiles and compose the view from them.** react-msa-viewer and
   Nightingale cache 10 x 10 residue tiles (200-entry LRU, 20 tiles of padding); react-msaview
   caches 512-cell raster tiles (LRU 64) for colours and 500 px canvas blocks for letters;
   igv.js paints 3x the viewport width and just moves the canvas with CSS `left` until the view
   leaves the painted range; AliView pre-renders 256 cell tiles per state and composes an
   `int[]` for the clip. Measured: cached per-row strips cost 7 ms per frame (54 `drawImage`s)
   vs 29-36 ms per-cell. Invalidate a strip when its row changes (edit, colour scheme).
3. **Glyph atlas instead of fillText.** Pre-render each (letter, background, foreground) cell
   once after `document.fonts.ready` and `drawImage` it; xterm.js keys atlas glyphs by
   background colour so antialiasing fringes match. Saves 15-20 % over per-glyph canvases
   and avoids fillText's spikes, but does not change the draw-call count. react-msaview's
   author found a letter texture cache "a minor part" of their speed-up (PR 79/80).
4. **Run-length and colour batching.** Jalview merges adjacent same-colour cells into one
   `fillRect`; SeaView and Jalview (monospaced fast path) draw a whole row's letters of one
   colour with one text call using space-masked strings. Cuts calls 3-5x on conserved DNA;
   depends on a true monospaced font whose advance equals the cell width.
5. **Drop letters below a cell size, never drop cells.** AliView bakes letters only above 6 px,
   UGENE needs a 7 px cell, Jalview requires `charWidth('M') <= charWidth`, igv.js 8 px rows,
   react-msaview 5 x 8 px. Every cell keeps its own rectangle down to 1 px (UGENE floors the
   column width at 1 px and lets the content be wider than the view). Below 1 px per cell
   AliView and Jalview's overview sample one residue per pixel, which drops residues; an exact
   zoom-out needs an aggregate per pixel column (none of the surveyed viewers computes one).
6. **Separate the content layer from the overlay.** UGENE draws its cached pixmap and then the
   selection frame on top, so selection changes never re-render cells; xterm.js uses four
   canvas layers (text, selection, link, cursor). Two stacked canvases in a browser.
7. **Coalesce input to one draw per frame.** Accumulate wheel/drag deltas and flush once in
   `requestAnimationFrame` (react-msaview `useWheelScroll.ts`); time-based throttles (BioJS msa
   30 ms, Wasabi 100 ms) only add latency. MDN: rAF-throttling `scroll` events does not reduce
   their rate, it is for correctness (https://developer.mozilla.org/en-US/docs/Web/API/Document/scroll_event).
8. **Build bitmaps off the main thread** when the composite is a pixel buffer: AlignmentViewer
   2.0 writes one RGBA pixel per residue into `ImageData` in a worker and transfers an
   `ImageBitmap` (`OffscreenCanvas.transferToImageBitmap`); JBrowse 2 renders tracks in
   workers and `drawImage`s the result. Adds message latency; only worth it when the draw
   itself does not fit the budget.

## 4. Scroll architecture

- **Where the scrollbar lives.** Chrome scrolls on the compositor thread and only then tells
  the main thread (https://www.chromium.org/developers/design-documents/compositor-thread-architecture/).
  A canvas *inside* the scrolled content moves with the thumb and then jumps when JS redraws
  it; a canvas *outside* the scroller (absolutely positioned underlay, repositioned from
  `scrollTop`) lags the thumb by one frame but never jumps. Glide Data Grid uses the second
  pattern: a `div` with `overflow: auto` filled with padder divs (each at most 5,000,000 px,
  total clamped to 33,554,400 px) and the canvas in a sibling underlay. react-msa-viewer /
  Nightingale's `FakeScroll` is the same idea (spacer-only scroller emitting deltas).
- **Maximum scroll extent.** Chrome 33,554,400 px, Safari 33,554,428, Firefox drops heights
  above 17,895,697 (https://meyerweb.com/eric/thoughts/2025/08/07/infinite-pixels/); causes:
  Chromium `LayoutUnit` is 1/64 px in an int32, Firefox `nscoord_MAX = (1 << 30) - 1` app units
  at 60 per px. A 5,000-row alignment at 14 px is 70,000 px, three orders of magnitude below.
- **Custom scrollbars** (AlignmentViewer 2.0, react-msaview, Wasabi, BioJS msa, ViewAlign
  Canvas mode) avoid the spacer but must reimplement thumb drag, keyboard, and accessibility,
  and the wheel listener must then be `{ passive: false }` to `preventDefault`, which blocks
  compositor scrolling of anything else. Chrome made document-level wheel listeners passive by
  default in Chrome 73 (https://developer.chrome.com/blog/scrolling-intervention-2/); a wheel
  listener on the element itself is not passive by default.
- **CSS transform vs redraw.** Moving an already painted canvas with `transform: translate`
  (react-msaview blocks, igv.js `style.left`) is compositor-only and free until the view leaves
  the painted area. `will-change: transform` promotes a layer at a memory cost; MDN calls it a
  "last resort" (https://developer.mozilla.org/en-US/docs/Web/CSS/will-change). `contain` and
  `content-visibility` act on DOM subtrees and do nothing for a single canvas.
- **Zoom.** react-msaview CSS-scales the stale bitmaps (max 2x) and redraws 150 ms after the
  wheel settles (`zoomSnapshot.ts`), so zoom never stalls.

## 5. Canvas 2D facts that matter

- Context attributes are fixed at the first `getContext`. `alpha: false` lets the browser
  treat the backdrop as opaque (and is a condition for LCD text antialiasing in Chrome and
  Firefox per a 2016 Chromium statement, https://www.w3.org/Bugs/Public/show_bug.cgi?id=23516,
  current status unverified). `willReadFrequently: true` forces a software canvas.
  `desynchronized: true` bypasses the compositor queue at the price of tearing and flicker
  (https://developer.chrome.com/blog/desynchronized). Do not set the last two.
- Chromium drops a canvas to the CPU after two `getImageData` readbacks
  (`kFallbackToCPUAfterReadbacks = 2`), leaves canvases under 128 x 129 px unaccelerated, and
  stops accelerating new canvases once a document has 100 canvases with 95 % of them
  unaccelerated (`html_canvas_element.cc`, `base_rendering_context_2d.cc`). Hit-test from data,
  never from pixels; do not make one canvas per row.
- Size limits (https://raw.githubusercontent.com/jhildenbiddle/canvas-size/master/docs/index.md):
  Chrome and Edge side 65,535 px, area 16,384^2; Firefox 122+ side 32,767, area 23,168^2;
  Safari area 16,384^2, iOS 4,096^2. A whole 600 x 600 alignment at 12 px is 7,200^2 =
  52 M px (legal, about 200 MB RGBA; 830 MB at DPR 2); 5,000 x 5,000 at 12 px is 60,000^2 and
  exceeds every limit. Never render the whole grid to one bitmap.
- Geometry in device pixels: backing store = CSS size x `devicePixelRatio`; xterm.js floors the
  cell width and ceils the height in device pixels and derives the CSS size from that, so no
  cell straddles a device pixel. MDN: avoid sub-pixel coordinates and scaling in `drawImage`
  (https://developer.mozilla.org/en-US/docs/Web/API/Canvas_API/Tutorial/Optimizing_canvas).
  Detect DPR changes with `matchMedia('(resolution: Ndppx)')`.
- Text: `document.fonts.ready` resolves when loading and layout finished; a font used only on
  canvas may need an explicit `document.fonts.load` first (unverified). `textBaseline = 'middle'`
  centres the em box, not the ink; xterm.js scans the glyph's alpha channel for the ink box.
  Grayscale vs LCD antialiasing differs across machines, so pixel-diff tests need a reference
  rendered on the same machine.
- OffscreenCanvas is Baseline since March 2023; `createImageBitmap` is asynchronous.

## 6. Data layout and garbage

- One `Uint8Array` per row (or one flat row-major array) of palette indices: 1 byte per cell,
  25 MB for 5,000 x 5,000. AliView uses `byte[]` per sequence; Jalview `char[]` (2 bytes);
  every browser viewer surveyed uses JavaScript strings. No primary benchmark of strings vs
  typed arrays for this access pattern was found (unverified); the known cost of strings is
  per-cell allocation when keys or substrings are built in the draw loop.
- Column-major mirrors (for per-column statistics) double memory; compute column statistics
  once at load and cache per column, not per cell per frame.
- GC: V8 scavenges the young generation often and runs idle-time GC before the next frame
  (https://v8.dev/blog/trash-talk); the VS Code renderer post attributes "noticeable delays"
  to DOM-churn garbage, xterm.js reuses two Float32Arrays. In a draw loop avoid per-cell
  string concatenation (cache keys), closures, and arrays.

## 7. Measuring honestly

- rAF timestamp deltas during a synthesized scroll (histogram, p95, frames over 34 ms) is the
  simplest truthful metric; JBrowse 2 measured frame rate the same way. Drive input with CDP
  `Input.dispatchMouseEvent` type `mouseWheel` (deltas in CSS px) at a fixed cadence, or
  `Input.synthesizeScrollGesture` (default speed 800 px/s), and throttle with
  `Emulation.setCPUThrottlingRate` because a developer machine is not the user's.
- DevTools Performance panel: Frames track (green rendered, yellow partially presented, red
  dropped), long tasks over 50 ms, Rendering tab "Frame Rendering Stats", "Paint flashing",
  "Scrolling performance issues" (https://developer.chrome.com/docs/devtools/performance/reference).
- Long Animation Frames report only frames over 50 ms with script attribution (Chrome 123+):
  good for stalls, blind to 16-34 ms jank. Event Timing excludes `wheel`, `mousemove` and
  `touchmove`, so wheel latency cannot be measured with it.
- A devtools trace splits the time into `Layout`, `Paint`, `Layerize`, `FireAnimationFrame`,
  `GPUTask`; the `invalidationTracking` category names the DOM node and reason for each
  style/layout invalidation (this is how the `content-visibility` cost in ViewAlign was found).
- For exactness, pixel-diff the viewport against a CPU-rendered reference of the same cells.

## 8. Recommended architecture for 600 x 600 to 5,000 x 5,000

1. Data: row-major `Uint8Array` of palette indices; per-column statistics cached once.
2. One viewport-sized canvas with `{ alpha: false }`, integer device-pixel cell geometry, no
   `willReadFrequently`, no `desynchronized`, no `getImageData`.
3. A cell atlas (letter x background x foreground x style) rendered once after fonts are
   ready, rebuilt on zoom or DPR change; `drawImage` with smoothing off.
4. Scroll by blitting the canvas onto itself and painting only the exposed strip; full redraw
   on diagonal moves, jumps over a page, zoom and colour-scheme changes. Optionally cache row
   strips so even the full redraw is 54 draws, not 10,600.
5. Input: native scroller with a spacer and the canvas in an underlay (or custom bars if the
   UI needs them), all sources coalesced into one dirty flag drawn per rAF.
6. Overlay (selection, cursor, search hits) on a second canvas so it never re-renders cells.
7. Measure with rAF deltas under CPU throttling; escalate to OffscreenCanvas in a worker and
   then to WebGL instanced quads only if the measured full redraw still misses the budget.

## 9. Pitfalls that look like optimisations

- Rendering the whole grid to one bitmap and translating it (area limits, hundreds of MB).
- `willReadFrequently: true` "to be safe" (software canvas); `getImageData` for hit testing
  (two readbacks drop to CPU); `desynchronized: true` for latency (tearing, flicker).
- `ctx.scale(dpr)` with fractional cell sizes (seams, resampling); blitting with smoothing on
  or non-integer offsets (blur accumulates).
- One canvas per row or tile under 128 x 129 px (unaccelerated, and can disable acceleration
  document-wide).
- A non-passive wheel listener on an element inside a native scroller.
- Time-based throttles on scroll or wheel (BioJS msa 30 ms, Wasabi 100 ms): latency, no savings.
- `will-change` left on permanently; `content-visibility: auto` on rows as a scroll fix (it
  moves layout and paint into the scroll).
- Letter atlas drawn on transparent and composited over colour (antialiasing fringes).
- Zoom-out by sampling one residue per pixel (drops residues) when exactness is required.
- Measuring wheel latency with Event Timing, 16 ms jank with Long Animation Frames, or only on
  an unthrottled developer machine.

## Sources
VS Code terminal renderer 2017 (https://code.visualstudio.com/blogs/2017/10/03/terminal-renderer), VS Code 1.41 notes (https://code.visualstudio.com/updates/v1_41), xterm.js addon-webgl `WebglRenderer.ts`, `GlyphRenderer.ts`, `TextureAtlas.ts`, README; Glide Data Grid README, `data-grid-render.blit.ts`, `infinite-scroller.tsx`; canvas-datagrid, Univer READMEs; Google Docs canvas announcement 2021; VS Code `viewLayer.ts`, `minimap.ts`, `minimapCharSheet.ts`; Figma blog 2015 and WebGPU post; Leaflet and MapLibre option docs; igv.js `trackViewport.js`; JBrowse 2 paper (PMC10108523), `offscreenCanvasPonyfill.ts`, "why not X" design notes; MDN Optimizing canvas, getContext, OffscreenCanvas, createImageBitmap, devicePixelRatio, scroll event, will-change, contain, FontFaceSet.ready, textBaseline, PerformanceEventTiming, Typed arrays; HTML spec canvas section; Chrome blogs desynchronized, scrolling interventions 1 and 2; Chromium `html_canvas_element.cc`, `base_rendering_context_2d.cc/.h`, `layout_unit.h`; Firefox `nsCoord.h`; canvas-size limits table; web.dev rendering-performance, canvas-performance, content-visibility; Chrome DevTools performance reference and rendering docs; Long Animation Frames docs; DevTools protocol `browser_protocol.json`; W3C bug 23516; V8 "Trash talk"; Eric Meyer "Infinite pixels" 2025. Viewer sources as listed in the landscape document.

## 10. Lessons from implementing the frame copy in ViewAlign (v244-v247, 2026-10-16)

- **Copy through a cleared scratch canvas, not onto yourself.** `drawImage(canvas, ...)` onto the
  same canvas composites source-over; where the destination still holds the previous frame,
  antialiased glyph pixels gain opacity on every copy (alpha 161 became 255 after six steps).
  Copy the frame into a cleared second canvas first, then clear the destination region and
  copy back. The pixel-identity check against a full redraw caught this at once.
- **Put cell edges on device pixels.** At 125 % or 150 % Windows scaling a 15 CSS px row is
  18.75 or 22.5 device px: every row edge is antialiased, the copy ends on a device pixel,
  and the row cut by that edge misses the neighbour's edge composite that a full redraw
  would apply. Fix: cell width and height are `ceil(css * dpr) / dpr` (xterm.js floors the
  width, ceils the height in device pixels; same idea), hairlines are `1 / dpr` CSS px, and
  the band heights (consensus rule + margin) are rounded the same way.
- **Widen the exposed strip to the nearest row or column edge** so every cell is either
  copied whole or drawn whole; the strip is cleared under a device-aligned clip before it is
  drawn. With both rules the copy is pixel-identical at 1x, 1.25x and 1.5x in every
  direction, including fractional trackpad deltas.
- **Keep the offsets on device pixels** (`round(offset * dpr) / dpr`), otherwise the copy
  delta and the strip's drawing phase disagree by a fraction of a pixel.
- **A safety net for state changes.** Compare a short signature of what can change a cell's
  look (sequence array identity and length, undo/redo stack lengths, conservation and
  consensus cache identities, selection sizes, trim previews, edit cell, colour scheme) with
  the previous frame; any difference forces a full redraw even if the caller forgot to mark
  the view dirty. Cost: one array compare per frame.
- **Read live caches each frame.** Edits reshade columns and recompute the consensus into
  state caches without a re-render; the draw must read those, not the arrays captured when
  the renderer was created, or copied frames keep stale shading.
- **Scroll-position inputs that bypass the scroller** (PageUp/PageDown, keyboard, scrollbar
  mirror) each need their own path into the same offset + draw request; the DOM
  `scrollTop` has no effect on a canvas viewport.
- **Measured result**: 600 x 600, Radeon 780M, headful Chrome 155: scroll frames p95 8.5-9 ms,
  none over 34 ms, in all four directions; a copy frame costs 1.3-4 ms against 25-50 ms for
  a full redraw.
