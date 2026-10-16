---
title: How MSA viewers render and scroll (landscape, 2026-10)
type: reference
---

> Type: reference

Source-level survey of how multiple-sequence-alignment viewers draw and scroll a plain
alignment: rendering technology, scroll mechanism, caches, in-memory layout, zoom-out
behaviour, and stated limits. Only the plain viewing experience; no analysis features.
Compiled 2026-10-16 for the ViewAlign responsiveness work
([projects/viewalign-scroll-responsiveness-2026-10.md](../projects/viewalign-scroll-responsiveness-2026-10.md));
the cross-cutting techniques and browser facts are in
[browser-grid-rendering-techniques.md](browser-grid-rendering-techniques.md).
Everything was read from the named source files unless marked "not verified".

## Browser viewers

### BioJS msa / MSAViewer (wilzbach/msa; archived 2022, README "This project is unmaintained")
- One `<canvas>` (`CanvasSeqBlock.js`). `CanvasSeqDrawer.drawLetters` makes two passes over
  the visible cells: `fillRect` per cell, then the letter by `drawImage` of a cached glyph
  (`CanvasCharCache`, keyed by letter only, one canvas per letter, wiped on cell-size change).
  Source comments: "REALLY expensive call on FF", `drawSeq` "TODO: very expensive method".
- Scroll: no native scroller; mousewheel + drag handlers move `_alignmentScrollLeft/Top` on
  the zoomer and trigger a full redraw (`this.el.width = this.el.width` to clear). No
  requestAnimationFrame; lodash `throttle(..., 30)` on non-WebKit.
- Data: Backbone model, sequence as a string. Zoom-out: letters skipped below
  `minLetterDrawSize` 11 px; the OverviewBox minimap is downsampled (modal colour per box),
  so not exact.
- Issues: #62 drawImage slow on Gecko (WebKit < 20 ms, Gecko > 45 ms per draw), #195 overview
  "horribly slow", #13 large alignments.

### react-msa-viewer (Plotly, 2018; repo now 404) and its port Nightingale MSA (EBI)
- Canvas, double-buffered: `DraggingComponent` keeps two stacked canvases and `swapContexts()`
  after each draw. The viewport is composed of **grid tiles** of 10 x 10 residues
  (`xGridSize`/`yGridSize`), each a cached canvas (`tileCache`, key `row-column`) built by
  `CanvasTilingGrid` from **residue tiles** (`residueTileCache`, key
  `${text}-${colorSchemeName}-${overlayFactor}`), blitted with `drawImage`.
  `cacheElements: 20` tiles of padding beyond the viewport; `CanvasCache` holds 200 entries,
  evicts oldest.
- Scroll: `FakeScroll` is an absolutely positioned div with a spacer of the full content size,
  so native scrollbars exist; its `scroll` event is converted to a movement delta. Drag is
  rAF-throttled (`utils/requestAnimation.js`, extra calls dropped). Plotly's version has
  `// TODO: handle wheel events`; Nightingale zooms through d3-zoom with a rAF-coalesced
  refresh. No `passive` listener options.
- Data: strings. Zoom-out: residue tile skips text below 4 px, fades to full at 10 px; no
  aggregation below 1 px (not verified). Claims: README "performant ... production-ready".
- Nightingale is used by UniProt/InterPro; issue #108 misalignment when zooming.

### AlignmentViewer 2.0 (Sander lab, alignmentviewer.org)
- Two layers. Colours: a worker (`MSAGenerationWorker.ts`) writes **one RGBA pixel per
  residue** into `ImageData`, per tile `new OffscreenCanvas` -> `transferToImageBitmap()`;
  tiles of 1,024 px for the WebGL path (PixiJS sprites, `scaleMode: NEAREST`) or one tile
  for the Canvas 2D path, drawn each frame with `imageSmoothingEnabled = false` and a scaled
  `drawImage`. Letters: `MSALetters.tsx` is **DOM text**, one `div` per distinct colour
  holding each visible row as a string with spaces where other colours sit; only the
  viewport slice is rendered.
- Scroll: no native scroller; React `onWheel` adds deltas to a Redux `worldOffsetPx`,
  content is shifted with `top/left`, custom scrollbar (DOM thumb + full-screen overlay while
  dragging). At most one partial cell per edge, no overscan, no rAF.
- Data: strings plus a per-column letter-count Map. Zoom-out: exact at >= 1 px per residue,
  NEAREST sampling below (inferred). README: "very large alignemnts (up to hundreds of
  millions of residues)", tested on about 23,000 x 250; "Inserting many elements into the
  HTML DOM at one time absolutely kills performance".

### react-msaview / JBrowse 2 MSA view (GMOD)
- A grid of `<canvas>` **blocks of 500 px** (`blockSize`, backing store scaled by
  `devicePixelRatio`), positioned absolutely; `renderMSABlock.ts` draws per cell `fillRect` +
  `fillText`, no glyph cache. With `cellFill === 'residue'` the background comes from
  **raster tiles** (`msaRaster.ts`): 512-cell tiles of 1 px per cell `ImageData`, drawn scaled
  with `imageSmoothingEnabled = false`, LRU of 64 tiles; when a cell is under one device pixel
  it averages 2^n cells per pixel.
- Scroll: all blocks sit in one wrapper moved by `transform: translate(...)` with
  `willChange: 'transform'`; only intersecting blocks exist (no overscan). `useWheelScroll.ts`
  registers `wheel` with `{ passive: false }`, **accumulates deltas and flushes once per
  requestAnimationFrame**, calls `preventDefault` only for consumed axes; drag likewise.
  Custom vertical scrollbar. Drawing runs in a MobX `autorun` (synchronous).
- Zoom: `zoomSnapshot.ts` CSS-scales the existing block bitmaps (max 2x) until the wheel
  settles (150 ms), dropping letters when a block has too many cells; letters need
  `minLetterRowHeight` 8 px and `minLetterColWidth` 5 px.
- Data: strings. Documented claim (gmod.org/JBrowseMSA/guide): "memory and draw time grow
  with the visible area, not with the size of the file"; 20,705 x 672 first paints in about
  4 s. Issues: #3 fast scrolling shows blank for seconds (17k sequences; closed 2024), #75
  slowdown at 5,000 x 5,000 fixed by smaller blocks (PR #80); PR #79 pre-rendered letter tiles
  closed unmerged, author: the letter texture cache was "a minor part" of the speed-up.

### Wasabi (veidenberg/wasabi)
- Many canvas tiles of 4,000 x 2,000 px appended as DOM elements (`canvaslist`, key
  `row|col`), each cell blitted with `drawImage` from per-symbol canvases made once
  (`makeCanvases()`, text only when font size > 7). Tiles are built in staggered
  `setTimeout`s (500 ms apart) and faded in.
- Scroll: custom scrollbars (malihu-based) set CSS `left` / `margin-top`; redraw debounced
  100 ms; no rAF. Zoom: discrete levels 1, 3, 5, 8, 11, 15, 20 px; at 1 px the tiles are
  colour-only (exact).

### MSABrowser (thekaplanlab)
- DOM: one absolutely positioned div per visible residue (`left: pos * 20px`), gaps skipped.
  Native scroll event re-windows the divs; prunes above 50,000 divs. Fixed 20 px cells, no zoom.

### igv.js (not an MSA viewer; tiling reference)
- The track canvas is painted **3x the viewport width** (and up to 3x the height); panning
  moves it with `canvas.style.left/top`; repaint only when leaving the painted extent or
  when the scale changes. Backing store scaled by devicePixelRatio. Letters via `strokeText`
  only at <= 0.1 bp per pixel and row height >= 8 px, else rectangles.

### NCBI Multiple Sequence Alignment Viewer (closed)
- ExtJS 5 front end. Release notes: 1.7 "Major rendering engine upgrade to support large
  size alignments", 1.8 smooth-animation scrolling and a fix for "back-end OpenGL (MESA)
  initialization" (suggests server-side rendering; not verified), 1.15 WebAssembly tooltips.
  Nothing public on canvas, tiles or limits.

### One-liners
MSAexplorer (2025): Python/Shiny + matplotlib/plotly, not a JS viewer. ProViz: rendering
not verified. Alvis: command-line, vector output. JSAV: jQuery DOM table. AlignmentViewer 1:
D3. pymsaviz, ggmsa: static images.

## Desktop viewers

### AliView (Java/Swing; Larsson 2014, doi:10.1093/bioinformatics/btu531)
- **Pixel-array compositing, not Graphics2D per cell.** For the current zoom level
  `CharPixelsContainer` holds `CharPixels[256]`, one pre-rendered ARGB tile per byte value
  with foreground and background baked in (`CharPixels.getRGBArray()` renders once into a
  `BufferedImage` and keeps the `int[]`). Separate containers per state (default, selected,
  consensus, translation; `CompoundCharPixelsContainer` = `HashMap<Color, ...>` for
  position-dependent colouring). Rebuilt when the cell size changes (`AlignmentPane`).
- `paintComponent` converts the clip rectangle to matrix coordinates (`paneCoordToMatrixCoord`,
  padded by one row/column), allocates an `int[]` for exactly that clip, fills it by copying
  tiles (`SequencePainterNucleotide.copyPixels` -> `ImageUtils.insertRGBArrayAt`), wraps it in a
  `MemoryImageSource` and issues **one `drawImage`**. No `drawString` in the paint path.
  Rows are painted in parallel (1-3 threads, one `SequencePainter` per visible row).
- Scroll: plain `JScrollPane`, default scroll mode (`BACKINGSTORE_SCROLL_MODE` commented out:
  "not working on Mac Retina-screen"), double-buffered; viewport-only repaint, no retained
  image, no `copyArea`. Unit increments 16 px vertical, 160 px horizontal.
- Storage: one `byte[]` per sequence (`DefaultBases`), dots turned into '-' on load; huge files
  are memory-mapped (`FileSequence`, `MemoryMappedSequencesFile`), painted single-threaded
  ("more threads make reading file slower").
- Letters are baked into the tile only when the cell is wider than `MAX_CHARSIZE_TO_DRAW = 6`
  px (times the HiDPI factor); below that the tile is a flat colour.
- Zoom-out below 1 px per cell: a loop per **screen pixel** samples one residue per pixel
  (`ySeq = y * seqPerPixY`, `seqXPos = x * seqPerPix`): nearest-neighbour subsampling, not
  exact.
- Claims (paper): "handles alignments of unlimited size"; loads 2-14x faster than Jalview /
  SeaView; 39,407 x 54,103 (2.1 GB) opened in 17 s with 2.2 GB RAM vs Jalview 73 s / 4.7 GB vs
  SeaView 51 s / 5.7 GB; 479,726 x 46,512 (22.4 GB) indexed in 88 MB RAM. Table 2 at 1920 x
  1200: visible at once AliView "Unlimited", Jalview 495 rows x 1,830 columns, SeaView 106 x
  305, ClustalX 120 x 345, Mesquite 68 x 1,650.

### Jalview 2.11.5 (Java/Swing)
- Graphics2D onto **one retained `BufferedImage`** sized to the viewport rounded down to whole
  cells (`SeqCanvas`). No glyph cache.
- **Fast scroll by self-copy.** `ViewportRanges` fires `STARTRES` / `STARTSEQ` changes;
  `SeqCanvas.fastPaint(dx, dy)` does `gg.copyArea(...)` on the back image, then `drawPanel`
  only for the exposed strip (the range narrowed to the delta, with `gg.translate`), guards
  against moving more than a page ("scrolling too fast, more than a page at a time" -> full
  repaint), sets `fastPaint = true`; `paintComponent` then just blits the image. Comments
  record a Swing pitfall (bundled repaints lost the flag) and a re-entrancy guard.
- `SequenceRenderer.drawSequence`: `drawBoxes` merges adjacent same-colour cells into one
  `fillRect` (run-length), then `drawText` only if `fm.charWidth('M') <= charWidth`. With a
  monospaced font, text shown, no groups and no per-residue text colour it draws a whole
  visible row with **one `drawString`**; otherwise one per character with glyph centring.
- Storage: `char[]` per sequence (2 bytes per residue), `List<SequenceI>`.
- Zoom-out: font shrinks until letters no longer fit, then boxes only (every cell still a
  rectangle of >= 1 px: exact). Whole-alignment view is the Overview window, which samples
  (`pixelsPerSeq`, `setRGB`), not exact.
- No numeric rendering claims in code or help; AliView's Table 2 measured 495 x 1,830 visible.

### UGENE (C++/Qt)
- `QPainter` into **one cached `QPixmap` of the viewport** (`MaEditorSequenceArea::drawAll`):
  on `completeRedraw` the pixmap is cleared and `drawVisibleContent` repaints every visible
  cell; otherwise the pixmap is just blitted and the selection frame / focus rect drawn on
  top, so selection changes never re-render cells. No glyph cache.
- Scroll: every scrollbar change sets `completeRedraw = true` and re-renders the whole
  visible region (`ScrollController` -> `si_completeRedraw`); no blit of the retained pixmap.
  Code comment: "// SANGER_TODO: optimize".
- `SequenceAreaRenderer::drawRow`: per cell `fillRect` then, only when text is visible,
  `drawText(QString(c))`: one QString and one drawText per residue.
- Storage: ungapped `DNASequence` plus a gap model `QVector<U2MsaGap>` per row; `charAt` goes
  through a 200-character gapped window cache with 10 % look-behind ("This access may be very
  slow for rows with a large gap models").
- Zoom: `resizeMode = zoomFactor < 1 ? OnlyContent : FontAndContent`; text needs a 7 px cell;
  column width floored at 1 px, so the alignment just gets wider than the view (exact).
- Docs: "Flexible zooming for large alignments"; no numbers.

### SeaView 5.1 (C++/FLTK)
- FLTK immediate drawing in `DNA_obj::draw()`, optional double-buffered window; no glyph
  cache, no retained image. Scrolling redraws all visible rows (`first_seq`, `first_site`);
  damage flags exist only for cursor / name / one-row edits.
- Text strategy: per row it builds up to 22 **colour-specific strings** (every residue not of
  colour i replaced by a space) and issues **one `fl_draw` per colour per row** (about 5 text
  calls per DNA row instead of one per residue); background `fl_rectf` per row.
- Storage: `char **sequence` plus a precomputed per-residue colour-class matrix
  `char **col_rank`. No zoom-out mode beyond font size. Measured by AliView at 106 x 305
  visible.

### CLC, Geneious, MEGA, BioEdit (closed)
Only documentation exists: CLC zooms "to any desired level"; Geneious and MEGA publish no
rendering limits; BioEdit is unmaintained (7.2.5, 2013). Rendering internals not verified for
any of them.

## Summary table (desktop)

| program | tech | caches | per-cell strategy | storage | claimed scale |
|---|---|---|---|---|---|
| AliView | int[] compositing -> MemoryImageSource | 256 pre-rendered tiles per state; rebuilt on zoom | tile blit, no drawString; letters only above 6 px; rows in parallel | byte[] per sequence; mmap for huge files | "unlimited"; 39k x 54k in 17 s / 2.2 GB |
| Jalview 2.11 | Java2D on retained BufferedImage | viewport back-image; copyArea on scroll, paint exposed strip only | run-length fillRect; text only if 'M' fits; one drawString per row when monospaced | char[] per sequence | none stated; 495 x 1,830 visible (AliView T2) |
| UGENE | QPainter into cached QPixmap | one viewport pixmap; overlay drawn without re-render; full re-render on scroll | fillRect + drawText(QString(c)) per cell; text off when zoom < 1 | ungapped sequence + gap model; 200-char window cache | "flexible zooming"; no numbers |
| SeaView | FLTK immediate | none | one fl_draw per colour class per row (space-masked strings) | char** + colour-class matrix | none; 106 x 305 visible |
| CLC / Geneious / MEGA / BioEdit | closed | not verified | not verified | not verified | none found |

## Summary table (browser)

| viewer | render | scroll | cache | data | rAF | zoom-out exact |
|---|---|---|---|---|---|---|
| BioJS msa | 1 canvas, fillRect + glyph drawImage | custom wheel/drag, full redraw, 30 ms throttle | letter glyphs | string | no | minimap = modal colour |
| react-msa-viewer / Nightingale | 2 canvases swapped; 10 x 10 grid tiles of residue tiles | drag (rAF), FakeScroll spacer | 200-entry tile + residue caches | string | drag only | no |
| AlignmentViewer 2.0 | 1 px/residue ImageBitmap from a worker + DOM letters per colour | onWheel -> Redux; top/left; custom bars | worker output by key | string | no | NEAREST below 1 px |
| react-msaview | 500 px canvas blocks; 512-cell raster tiles | wheel passive:false + rAF accumulate; wrapper translate; custom bar | 64 raster tiles; zoom snapshot | string | yes | averaged 2^n cells |
| Wasabi | 4,000 x 2,000 canvas tiles, symbol drawImage | custom bars set left/margin; 100 ms debounce | symbol canvases, tile list | not verified | no | 1 px colour |
| MSABrowser | div per residue | native scroll, windowed | div map | string | no | n/a |
| igv.js | 3x viewport canvas, CSS left/top | shift vs repaint | feature cache | n/a | no | n/a |
| ViewAlign (for comparison) | DOM span per residue below 500k residues; viewport canvas with per-glyph canvases above or on request | document scroll (DOM) / wheel -> offset + rAF (canvas) | glyph canvases keyed letter+colours+variant | string per row | canvas only | n/a |

Common gap: every browser viewer read stores rows as JavaScript strings; none uses typed
arrays. None blits the previous frame and paints only the exposed strip; the ones that scroll
smoothly do so by painting beyond the viewport (igv.js, react-msaview blocks) or by caching
tiles (react-msa-viewer, react-msaview raster).

## Sources (desktop section)
AliView/AliView master: src/main/java/aliview/gui/pane/AlignmentPane.java, CharPixels.java, CharPixelsContainer.java, CompoundCharPixelsContainer.java, RGBArray.java, ImageUtils.java, SequencePainter.java, SequencePainterNucleotide.java; aliview/sequences/DefaultBases.java, InMemorySequence.java, FileSequence.java; sequencelist/MemoryMappedSequencesFile.java; AliViewWindow.java. Paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC4221126/.
Jalview 2.11.5.0 (Debian source mirror sources.debian.org/data/main/j/jalview/2.11.5.0+dfsg-2/): src/jalview/gui/SeqCanvas.java, SeqPanel.java, AlignmentPanel.java, SequenceRenderer.java, viewmodel/ViewportRanges.java, datamodel/Sequence.java, Alignment.java, renderer/OverviewRenderer.java; https://www.jalview.org/help/html/memory.html.
ugeneunipro/ugene master: src/corelibs/U2View/src/ov_msa/MaEditorSequenceArea.cpp/.h, SequenceAreaRenderer.cpp/.h, MsaEditorSequenceArea.cpp, ScrollController.cpp, DrawHelper.cpp, MaEditor.cpp/.h; src/corelibs/U2Core/src/datatype/msa/Msa.h, MsaRow.h/.cpp; https://ugene.net/docs/alignment-editor/overview/alignment-editor-features/ and .../zooming-and-fonts/.
SeaView 5.1 (sources.debian.org/data/main/s/seaview/1:5.1+20251106-1/): seaview.cxx, seaview.h.
CLC: resources.qiagenbioinformatics.com manuals View_alignments / Zoom_functionality_in_View_Area. Geneious help centre (algorithm sizes only). BioEdit help PDF not readable (HTTP 403).

## Sources (browser section)
wilzbach/msa: src/views/canvas/CanvasSeqBlock.js, CanvasSeqDrawer.js, CanvasCharCache.js, src/model/Sequence.js, src/g/zoomer.js, src/views/OverviewBox.js, issues 62/195/13.
galaxyproject/react-msa-viewer (fork of plotly's): src/components/Canvas/SequenceViewer.js, CanvasComponent.js, DraggingComponent.js, CanvasCache.js, src/utils/requestAnimation.js; npm react-msa-viewer 1.0.1.
sanderlab/AlignmentViewer2.0: src/components/msa-blocks-and-letters/*, src/components/virtualization/*, src/webworkers/MSAGenerationWorker.ts, src/redux/VirtualizationReducers.ts, README; https://alignmentviewer.org/.
GMOD/react-msaview: packages/lib/src/components/msa/MSACanvas.tsx, MSACanvasBlock.tsx, renderMSABlock.ts, msaRaster.ts; packages/lib/src/model.ts, constants.ts, calculateBlocks.ts, zoomSnapshot.ts, useWheelScroll.ts, useDragScroll.ts; https://gmod.org/JBrowseMSA/guide; issues 3, 75, PRs 79, 80.
ebi-webcomponents/nightingale: packages/nightingale-msa/src/components/Canvas/CanvasTilingGrid.ts, CanvasCache.ts, SequenceViewer.ts, FakeScroll.ts, DraggingComponent.ts; nightingale-new-core withZoom.
veidenberg/wasabi: javascript/script.js, javascript/lib/scrollbars.js. thekaplanlab/MSABrowser: javascript/msabrowser.js. igvteam/igv.js: js/trackViewport.js, trackView.js, browser.js, bam/alignmentTrack.js.
NCBI: https://www.ncbi.nlm.nih.gov/tools/msaviewer/release-notes/, .../embedding-api/.
