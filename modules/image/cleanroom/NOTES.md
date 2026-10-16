# Analyze Skeleton clean-room implementation — notes and evidence log

Implementation: `cleanroom/askel_cr.js` (single file, no dependencies, Node + browser/worker; export shape of SPEC §10, `PRUNE: true`).
Sources used: `SPEC_analyze_skeleton.md`, the imagej.net prose documentation page (pruning options are described there only in
words), the oracle run as a black box (`SkelDriver` on `/tmp/ij.jar:/tmp/askel_out`, never inspected), `test/test_skeleton.js`,
`test/load_core.js`, `ref/askel/SkelDriver.java`, `src/img_io.js` (TIFF writer), `src/core/30_binary.js` (via load_core, used only
to generate skeletonized test inputs). No Fiji/AnalyzeSkeleton source, earlier port, built HTML, git history or forum code was opened.

All experiment scripts are in `cleanroom/exp/` (`oracle.js` runs the oracle; `cmp.js` compares full outputs; `gen.js` makes random
volumes with a proper PRNG, because the harness LCG has a short cycle in double arithmetic, so its random cases repeat).

## 1. Design summary

- **Indexing.** Foreground voxels get compact ids in the specified scan order (z, x, y). A padded `Int32Array` id map gives
  branch-free 26-neighbour access; the 26 (8 in 2D) offsets are stored in *pick priority* order
  (dx = 1,0,−1; dy = 1,0,−1; dz = −1,0,1), so `pick` from the spec = "first hit in that list".
- **Labelling.** Classification by neighbour count; 26-connected components by BFS; tree numbering by first end-point / junction /
  slab in scan order; junction clusters with the exploration order of SPEC §4.
- **Graph construction.** Per tree: explicit vertex list (clusters first, then vertices in creation order), phases A/B/C walks of
  SPEC §5. All slab lists live in one shared `Int32Array` (edges keep `[s0, s1)` ranges), so memory is linear.
- **Metrics.** Triples/quadruples, averages, edges as in the spec.
- **Pruning** (section 3): end pruning on the first graph; cycle detection by an edge-stack depth-first search; cuts; re-analysis.
- **Longest shortest path** (section 2): exact matrix procedure for n ≤ 600 vertices (`FW_MAX`); above that, the same quantities
  without the matrix.

## 2. Longest shortest path (SPEC §7)

### 2.1 Spec incomplete: the maximum can be a diagonal entry (self-loop)
The spec says that the start point falls back to the previous tree's start point "if the path list is empty (best stayed 0)". The path
list is also empty when the best entry lies on the diagonal (a = b), which happens when the last edge touching a vertex is a self-loop:
the diagonal starts as the self-loop length and the relaxation lowers it to min(loop, 2 × shortest incident edge).
- `exp/e01_diag.js`: loop on a junction cluster + 1-slab tail. Oracle: `T 0 2 1 1 4 10 8.5355… 1 0 13.6568… 6.82842712474619 0.0 0.0 0.0`,
  edges: tail 3.4142, loop 13.6568 (self-loop, created last). lsp = 6.828 = 2 × 3.414 (diagonal), start point (0,0,0), LSPIMG empty.
- `exp/e02_diag_prev.js`: same with a vertical line tree before it (cal 0.5, 0.7): tree 0 start (0, 2.8, 0); tree 1 lsp 3.7205
  (= 2 × 1.8602), start (0, 2.8, 0), the previous tree's start point.
Conclusion: when the best pair is (r, r), lsp = that diagonal value, no path voxels, and the start point is the previous tree's.
The spec's matrix semantics already produce this. Only the parenthetical remark in the spec is incomplete.

### 2.2 No size limit (improvement)
For n > 600: build the simple graph with "last edge wins" weights (matrix overwrite semantics). Diagonal values are
min(self-loop length if it is the last edge touching the vertex, else 0; 2 × smallest incident weight). The matrix is exactly
symmetric under in-place relaxation, so Fiji's choice (first strictly greater entry, row-major) is: a = smallest vertex whose
eccentricity (or diagonal) reaches the maximum M, and b = smallest partner of a. Eccentricities come from:
- the block-cut tree (Tarjan biconnected components, Dijkstra inside each block, rerooting DP over the block-cut tree), when
  Σ block² ≤ 4·10⁶. Trees and tree-like skeletons are linear time;
- otherwise (big biconnected blocks such as meshes), exact eccentricity bounding (Takes–Kosters style lower/upper bounds with
  alternating sweeps), followed by resolving the smallest-index candidate.
Then one Dijkstra from a gives b, the distance and the predecessor path. `warn` gets an entry when more than two vertices may reach
M within 1e-12 relative (the start point may differ from Fiji), or when an equal-length alternative exists along the chosen path
(path voxels may differ). These ties depend on Fiji's relaxation order and floating-point association, so they are reported
rather than reproduced.

Evidence:
- `exp/e03_lsp_unit.js`: 4000 random multigraphs (2–61 vertices, multi-edges, self-loops, shuffled creation order, 30 % integer
  weights to force ties), non-matrix vs matrix: 0 unflagged mismatches with both the block-cut method and with forced bounding
  (`BB=0`); 706 flagged, 648 of them were still identical.
- Harness forced onto the non-matrix path (`exp/force_large.js`, FW_MAX = 1, and also BLOCK_BUDGET = 0), FUZZ 3000, SEED 5 and
  SEED 8: 3048 / 3050 pass; the 2 differences (`figure 8`, `grid (loops)`: LSP image only, symmetric equal-length loops) both
  carried the tie warning.
- Oracle on large graphs (`exp/e04_large_oracle.js`, default thresholds, so the non-matrix path is used): random trees with 1232,
  1220, 1684, 2137, 3108, 3455 and **8494** vertices (the last one took the oracle 1212 s); skeletons of blobs with holes (loops)
  with 1240–1296 vertices per tree: **all identical** (lsp, start point, LSP image), no warnings.
- `exp/e41_mesh.js MODE=bounding`: two jittered-grid mesh skeletons (≈900 vertices, one big block; cal 0.73–0.83 × 1.137) via
  the bounding method: identical to the oracle.

## 3. Pruning (SPEC §8): derived behaviour

### 3.1 Overall pipeline (`e10`–`e13`, `e16`, `e39`)
1. Analyse the input (A0).
2. **Prune ends** (if `pruneEnds`): in vertex-list order, every end-point vertex that still has a branch loses it. The branch's
   slab voxels and the end-point voxel are removed; the opposite vertex stays. For a branch between two end-points, only the
   first one in the vertex list is removed: `line` → only (6,1) remains, tagged 30; `twopix` → (2,1) remains.
3. **Cycle pruning** (if prune ≠ none): cycles are searched on the graph left after step 2 (see 3.3).
4. If at least one voxel was cut: the whole image, minus the cut voxels **and minus the end-pruned voxels**, is analysed again from
   scratch. Every output (trees and their order, edges, tags, LSP) comes from this re-analysis, and ends are *not* pruned
   again (`theta` with prune 1 + ends: the end-points created by the cuts remain).
   - `e16_reanalysis.js`: for prune 1/2/3 × 300 random volumes (236 with cuts), the oracle output = plain analysis of
     (input − removed voxels): 0 mismatches.
   - `e39.js` (two trees: theta with tail + pure theta, prune 1 + ends): the left tree is re-analysed without its tail
     (3 branches, avg 8.319 = (9.30 + 6 + 9.66)/3, 0 end-points), the right one is cut, and the tree order changes
     (the cut tree now has end-points and becomes tree 0).
5. If nothing was cut, the output is A0 with the end pruning applied (no re-analysis):
   - removed voxels are 0 in the tag image; removed edges are not listed;
   - `branches`, `avgLength`, `maxLength` and `junctions` are unchanged (A0 values, which include the pruned branches);
   - `endPoints`, `slabs` (= A0 count − removed), `triples` and `quadruples` (neighbour counts ignoring removed voxels) are recomputed;
   - LSP is computed on the reduced graph (removed vertices and edges dropped, vertex order kept).
   Evidence (`e10`, `e11`): cross: `T 0 4 1 0 5 0 3.457… 0 0 5.0 0.0 4.0 4.0 0.0`, no E lines (before: 4 edges, quadruples 1);
   H: one bridge edge left, lsp 5.0 from the second cluster (8,3); loop + tail: only the self-loop remains, lsp 0 with start
   (6,2) (a single vertex now, whereas the unpruned graph had the diagonal case of 2.1).
   I could not separate "slab counter − removed slabs" from "remaining slab voxels": the two were equal in 12 000 random
   volumes (`e42.js`).

### 3.2 Why prune ends usually suppresses cycle pruning
The search starts from the tree's **root** = the first non-cluster vertex of the vertex list (the first end-point vertex; a
cycle component's vertex), or the first cluster when there is none. It does not restart from other vertices. End pruning always
removes that root when a tree has end-points, so such trees get no cycle cuts (`e12` loop + tail, `e13` thetatail/loop2j:
prune 1 + ends = ends only), while trees without end-points are still cut (`theta`, `e39` right tree).
Root variants on 1000-case sets (`e20.js`): "first cluster" 1 and 4 mismatches; "first end vertex" and "v1 of the first
edge" 0 and 0.

### 3.3 Cycle detection (`e17`, `e18`, `e20`)
Depth-first search with an explicit stack of (vertex, via edge, from vertex). The root is pushed first. On a pop:
- if the edge is already labelled, skip;
- if the vertex is already visited, the edge is a BACK edge (from → vertex);
- otherwise it is a TREE edge, and the vertex records its parent and parent edge.
A newly visited vertex pushes all its unlabelled incident edges in creation order.
For each back edge, in discovery order, the cycle list = [back edge, then the parent edges from `from` up to the back edge's
target]. Recursive DFS variants left 6 / 400 mismatches; this edge-stack DFS left 0 / 400 (seed 3) and 0 / 2000 (seeds 11, 12).
Decisive case (`e19.js`, seed 3 #126): three parallel edges with 3, 10 and 7 slabs between two clusters. The oracle cut the 3-slab
and the 7-slab edges, which requires the tree edge to be the 10- or 7-slab edge, not the first one.

### 3.4 Which voxel is cut
- **Cycle component** (pure slab ring): the start slab (the ring's vertex voxel) is always removed, by every method.
  Evidence: `e12` ring → (2,2) = index 0 of the slab list; `e15` lenVsSlab (16 slabs) → index 0; prune 2/3 mismatches on
  rings disappeared with this rule (`e36.js`: the 3-voxel ring at (4,3)… cut at its start slab (4,3)).
- **Shortest branch**: the branch with the fewest **slabs** (not the shortest length), first in cycle-list order on ties (strict <).
  It is cut at slab index ⌊n/2⌋ (`e14_mid.js`: self-loops with 8/10/12/14/18 slabs → 4/5/6/7/9; 3 → 1; 1 → 0).
  `e17.js`: slab count 0 / 400 mismatches; length 20–30 / 400. Example: ring 10,10,4 cut the 5.61-long 1-slab edge before the
  3.0-long 1-slab edge.
- **Intensity of a voxel**: the mean of the original image over a box centred on the voxel. Half-sizes:
  rx = round(R/pw), ry = round(R/ph), rz = round(R/pd), with R = max(pw, ph, pd) (pd included even in 2D; rz = 0 in 2D) and
  round = ⌊x + 0.5⌋. Out-of-bounds voxels count as 0, and the divisor is always the full box size.
  Evidence:
  - `e23` (single dark slab k → slab k − 1 removed = first of the three 3×3 windows that contain it);
  - `e24`/`e27` (bright/dark combinations reproduce the 3×3 *mean including the centre*);
  - `e25` (dark background voxel at distance 2 has no effect);
  - `e32`/`e33` (pixel sensitivity: with ph = 1.93, pixels 2 columns away matter);
  - `e34` (window extents for 9 calibrations: (1,1.93) → rx 2; (1,2.6) → 3; (1,3.2) → 3; (0.5,1) → 2; (0.3,1) → 3;
    (2,1) → ry 2; (1,0.4) → ry 3; (0.7,1.3) → rx 2);
  - `e35` (zero padding: 0 / 300; in-bounds-only mean: 84 / 300; R without pd in 2D: 6 / 300).
- **Lowest intensity voxel**: the darkest slab over the cycle (edges in cycle-list order, slabs in list order, strict <).
  Junction voxels are not candidates (including them gave 434 / 600 mismatches).
- **Lowest intensity branch**: the branch with the lowest mean of voxel intensities over [v1's voxels (all cluster voxels),
  its slabs, v2's voxels (if v2 ≠ v1)], strict <. It is cut at its darkest slab. Mean over slabs only: 6 / 300; with the
  vertex voxels: 0 / 300 (`e37.js`).
- Validation (`e38.js`, 1000 new cases per mode, seed 77): prune 1, 2, 3: 0 mismatches in the removed-voxel sets.
- Missing original image with prune 2/3: the oracle aborts with an exception (no output). This implementation adds a `warn`
  entry and skips cycle pruning (end pruning still applies).

### 3.5 Full-output validation of pruning
- `exp/e40_prune_fuzz.js` (all 4 prune modes × prune ends on/off, mixed 2D/3D, random calibrations and intensity images, compared
  in full, including tags and the LSP image): seeds 101–104 × 1500 cases, 0 mismatches; 201 × 800, 0; with the final file
  901 × 1500 (mixed) and 902 × 1000 (prune ends always on): 0 mismatches.
- `e43.js`: 600 cases with several end-points attached to the same cluster without slabs (zero-length-walk edges), prune modes
  with/without ends: 0 mismatches.
- Harness with PRUNE enabled: see section 5.

## 4. Spec items confirmed unchanged
SPEC §2–§6 and the §7 matrix semantics passed as written: 50 fixed cases and 3 × 3000 fuzz cases, without any change to the
first version of the walk code. Empty images, 1×1, 3×1, 1×1×3, full 3×3×3 and non-positive calibrations (→ 1) also match the
oracle with all prune modes.

## 5. Acceptance results (final file)
- `SKEL_IMPL=cleanroom/askel_cr.js node test/test_skeleton.js`: **102 / 102** pass (50 base + 52 pruning cases, which are
  enabled because PRUNE = true).
- FUZZ = 3000, SEED = 1, 2, 3 (and 11, 12, 13): **5852 / 5852** pass for each seed.

## 6. Timings (Node 22, this container, single thread)
| case | time |
|---|---|
| 1500×1500 fully foreground (one 2.25 M-voxel cluster), no LSP / with LSP | 0.76–0.94 s / 0.73–0.90 s |
| 2000×2000 serpentine single path (≈1 M voxels) + LSP | 0.40–0.60 s |
| 256³ sparse random walks (198 k voxels, 142 trees) + LSP | 0.71–0.98 s |
| 256³ 2 % random noise (335 k voxels, 255 k trees) + LSP / + shortest-branch + ends | 1.55 s / 2.73 s |
| comb skeleton, one tree with ≈52 k vertices + LSP | 0.59–0.93 s |
| grid mesh 1300² (≈21 k vertices in one biconnected block) + LSP | 0.65–0.72 s (block-cut-only version: 61 s) |
| same mesh with prune 1 / 2 + ends / 3 (+ LSP) | 1.04 / 0.60 / 1.15 s |
| exact matrix LSP, n = 300 / 600 | 59 ms / 358 ms |

## 7. Remaining known differences / limitations
- LSP for graphs above 600 vertices: the start point and path can differ from Fiji only when the maximum is (nearly) tied or
  an equal-length alternative path exists. Such cases always carry a `warn` entry. A small n always uses the exact matrix
  procedure.
- `lspImage` value 255 after end pruning without a re-analysis: removed voxels are left at 0. The driver prints only the
  96-valued voxels, so this cannot be checked against the oracle.
- Unverifiable details that never showed up in 12 000+ random volumes: end vertices with self-loops (impossible geometrically,
  as far as I can tell); slab counter ≠ slab voxels (see 3.1).
- Fiji's own limitation (nested cycles may survive a single prune pass, per the documentation) is reproduced on purpose; there
  is no iteration.
