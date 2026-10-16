# Analyze Skeleton (2D/3D) — behavioural specification for a clean-room implementation

Purpose: an independent JavaScript implementation whose **outputs are identical** to Fiji's *Analyze Skeleton (2D/3D)*
(Arganda-Carreras et al., *Microsc. Res. Tech.* 73:1019–1029, 2010), written **without** access to the Fiji source code.

This document describes *what the outputs are* (observable behaviour, confirmed by running the Fiji jar as a black box).
It deliberately says nothing about how Fiji's code is organised. Choose your own data structures and decomposition.

## 0. Clean-room rules (mandatory)

- Do **not** open, search for, fetch or decompile: the Fiji / `sc.fiji.analyzeSkeleton` source (GitHub, Maven, jar internals, any mirror),
  any earlier JavaScript port of it (including anything in this repository's history, built HTML files, or scratch directories),
  or forum/blog posts that quote its code.
- Allowed: this specification; the published paper; the user documentation page (https://imagej.net/plugins/analyze-skeleton) as prose;
  running the compiled oracle as a black box (`java -cp /tmp/ij.jar:/tmp/askel_out SkelDriver ...`, see §9) on any input you construct;
  the test harness `test/test_skeleton.js`; general algorithm knowledge (connected components, Dijkstra, tree diameter, …).
- If a behaviour is not described here, find it by black-box experiments (construct small inputs, run the oracle, compare). Record what
  you found in `cleanroom/NOTES.md` (inputs used and the observed outputs) — that file is the evidence that the rule was derived from behaviour.

## 1. Input / output

`analyzeSkeleton(vol, cal, opt)`
- `vol = {w, h, d, data: Uint8Array(w*h*d)}`, voxel index `i = x + w*(y + h*z)`; non-zero = foreground (skeleton voxel).
- `cal = {pw, ph, pd}` voxel size (non-positive or missing → 1).
- `opt = {lsp: bool, prune: "none"|"shortest branch"|"lowest intensity voxel"|"lowest intensity branch", pruneEnds: bool, orig: {data: Uint8Array}|null, progress?: f(0..1)}`.

Returns
```
{ trees: [ {branches, junctions, endPoints, junctionVoxels, slabs, avgLength, triples, quadruples, maxLength,
            lsp, spx, spy, spz   // only when opt.lsp
          } ...],                // one per skeleton, in skeleton order (§3)
  edges: [ {tree, length, v1:[x,y,z], v2:[x,y,z], euclid, lengthRa, nSlabs} ... ],  // in creation order (§5), tree = 0-based
  tag: Uint8Array(N),            // 30 end-point, 70 junction voxel, 127 slab, 0 background
  lspImage: Uint8Array(N)|null,  // when opt.lsp: 255 on every foreground voxel, 96 on path voxels (§7); else null
  warn: [string] }
```
Distances: `d(p,q) = sqrt(((xp-xq)pw)^2 + ((yp-yq)ph)^2 + ((zp-zq)pd)^2)` in double precision. `spx, spy, spz` are calibrated (x·pw, …).
`v1/v2` are voxel coordinates (uncalibrated) of the *first voxel* (§4) of each vertex; `euclid = d(first(v1), first(v2))`.

## 2. Voxel classes

Neighbourhood = 26 voxels (in 2D, d = 1, this reduces to 8). Out-of-bounds counts as background.
`n` = number of foreground neighbours: `n < 2` → **end-point** (30), `n > 2` → **junction voxel** (70), `n = 2` → **slab** (127).
An isolated voxel is an end-point.

**Scan order** (used for every ordered list below): z outermost, then **x**, then **y** innermost
(i.e. within a slice, column by column). Lists `E`, `J`, `S` = end-points, junction voxels, slabs in scan order.

**Neighbour pick rule** `pick(p, predicate)`: among the in-bounds foreground 26-neighbours q of p that satisfy the predicate,
consider offsets (dx, dy) in lexicographic order dx = −1..1 (outer), dy = −1..1 (inner); for each (dx, dy) take the hit with the
smallest dz (−1 before 0 before 1). The result is the hit belonging to the **last** (dx, dy) that has any hit; "none" if no hit.
Every "next voxel" choice in this spec uses this rule.

## 3. Skeletons (trees)

A skeleton = a 26-connected component of foreground voxels.
Numbering: first the components that contain an end-point, ordered by their first end-point in `E`; then components without
end-points that contain a junction voxel, ordered by their first junction voxel in `J`; then the remaining components (pure cycles
of slabs), ordered by their first slab in `S`. A component started from a slab in this way is a **cycle component**; its **start slab** is
that first slab.

Per tree: `endPoints` = number of end-points, `junctionVoxels` = number of junction voxels.

## 4. Junctions (vertices made of junction voxels)

Within a tree, junction voxels are grouped into **clusters** = connected components under 26-adjacency *through junction voxels only*.
`junctions` = number of clusters. Cluster order = order of their first voxel in `J`.

Each cluster has an ordered voxel list. Its first element (the **first voxel**) is its earliest voxel in `J`. The rest of the order is the
discovery order of this exploration (it matters for §5 phase B):
```
list = [seed]; mark seed; queue = [seed] (read from the front, never shrinks except by advancing a head index)
next = pick(seed, unmarked junction voxel)
repeat while next exists or the queue is not exhausted:
   if next exists: append next to list and queue; mark it; next = pick(next, unmarked junction voxel)
   else:           next = pick(queue[head], unmarked junction voxel); if none: advance head
```

An end-point vertex = a vertex consisting of one end-point voxel. Each time the rules below "make an end vertex" a *new* vertex object is created.

## 5. Branches (edges)

All voxels start *unvisited* for this section. Per tree, three phases in this order. `cluster(q)` = the cluster that contains junction voxel q.

**Walk W(s)** (s = the voxel where the walk starts; returns length L, running-average length R, a *terminal voxel*, and maybe a *far vertex*):
- mark s visited; history = [s]; L = R = 0; next = pick(s, unvisited).
- If no next: return L = 0, R = 0 with **no** terminal and no far vertex.
- While next is a slab: increment the tree's slab counter; append next to the branch's slab list; L += d(prev, next); push next to history
  and add the running-average step (below); mark it visited; prev = next; next = pick(prev, unvisited).
- If next exists (it is an end-point or junction voxel): L += d(prev, next); push to history and add the running-average step; mark it visited.
  If it is an end-point → far vertex = new end vertex {next}. If it is a junction voxel → far vertex = cluster(next) and
  **L += d(first(cluster), next) and R += the same amount**. Terminal = next.
- If there is no next: terminal = prev (a slab); far vertex undefined.
- Running-average step: after the k-th push (k = 1, 2, …; history[0] is s), with m = min(k, 5): R += d(v_k, v_{k−m}) / m, where v_j is the j-th pushed voxel.

**Junction neighbour search** `J(t, exclude)`: q = pick(t, *visited* junction voxel not in `exclude`).

**Phase A — from end-points.** For each end-point e of the tree in `E` order that is still unvisited:
create vertex A = {e} and add it to the tree's vertex list (§7); branch slab list = []; run W(e).
- If L = 0: q = J(e, {e}). If q exists: edge (A, cluster(q)) with length = d(e, q) + d(first(cluster(q)), e)  ← note: the second term uses e, not q;
  R = 0; slab list empty; counts as a branch; added to the total length; **not** considered for `maxLength`. If q is none: no edge.
  Continue with the next end-point.
- Else if the terminal is a slab t: q = J(t, {e}). If q exists: far vertex = cluster(q), L += d(q, t) + d(first(cluster(q)), q).
  If none: far vertex = A itself and L += d(e, t)  (i.e. d(first(A), t)).
- Add the far vertex to the vertex list if new; edge (A, far) with L, R, slab list; branch count +1; total += L; maxLength = max(maxLength, L).

**Phase B — from junctions.** For each cluster V (cluster order), for each voxel j of V in V's list order:
mark j visited; then repeatedly n = pick(j, unvisited) until none:
- If n is a junction voxel: mark it visited (no branch).
- Else: slab list = [n]; increment the tree's slab counter; L = d(j, n); run W(n) and add its L; R = W's R (the step j→n is not in R).
  If W returned no terminal, the terminal is n. If the terminal is a slab t: q = J(t, voxels of V); if q exists far = cluster(q) and L += d(q, t);
  if none far = V (L unchanged). (If W ended on a junction, far = that cluster as set by W.)
  Then L += d(first(V), j). Edge (V, far); branch count +1; total += L; maxLength = max(maxLength, L).

**Phase C — cycle component.** If the tree is a cycle component: vertex C = {start slab} (added to the vertex list); slab list = [start slab];
slab counter +1; run W(start slab). If L ≠ 0: branch count +1, total += L, maxLength = max. **Always** add edge (C, C) with L, R, slab list
(the closing step back to the start is not part of L).

Per tree: `branches` = branch count; `slabs` = the slab counter (not necessarily the number of slab voxels); `avgLength` = total / branches
(0 when no branches); `maxLength` as accumulated (0 if none). Edges are output in creation order, tree by tree;
`nSlabs` = slab-list length; `length` = L; `lengthRa` = R.

## 6. Triple / quadruple points

For each cluster: nb = Σ over its voxels p of the number of in-bounds 26-neighbours of p whose class is slab or end-point
(a neighbour adjacent to two cluster voxels counts twice). nb = 3 → triple point, nb = 4 → quadruple point. Classes are those of §2.

## 7. Longest shortest path (opt.lsp)

Graph per tree: vertex list order = all clusters in cluster order, then vertices in the order they are added in §5.
Edges = the tree's edges. Matrix semantics (n = number of vertices):
- If n = 1: lsp = 0, start point = first voxel of that vertex.
- Else: A = n×n all +∞, P = n×n all −1. For each edge in creation order with endpoints r, c: A[r][r] = 0; A[c][c] = 0; A[r][c] = A[c][r] = length
  (a later edge between the same pair overwrites; a self-loop edge therefore leaves A[r][r] = its length); P[r][r] = −1; P[c][c] = −1; P[r][c] = r; P[c][r] = c
  (in that order). Then the all-pairs relaxation with strict improvement: for k, for i, for j: if A[i][k] + A[k][j] < A[i][j] then A[i][j] = that sum,
  P[i][j] = P[k][j] (in place, k outer).
- best = 0; scan i ascending, j ascending; take entries with finite value **strictly greater** than best → (a, b).
- Path: walk from b back to a: while current ≠ a: pre = P[a][current]; among the tree's edges joining pre and current choose the one with the smallest
  length (first in creation order on ties). Orient it so that it goes from current to pre: the vertex on the current side contributes its voxels
  (those not already in the path list), then the edge's slab voxels in the direction current→pre (marked 96 in `lspImage`), then the pre-side vertex's
  voxels (not already listed). current = pre.
- lsp = best; start point = the first voxel in the path list (= first voxel of vertex b). If the path list is empty (best stayed 0) the start point is
  **the previous tree's start point** (0,0,0 before the first tree).
- `lspImage`: 255 on all foreground voxels, then 96 on every slab voxel placed on any tree's path. Junction and end voxels are not marked.

**Required improvement over Fiji (no size limit).** Fiji's all-pairs relaxation is O(n³). Your implementation must reproduce the semantics above
exactly, and stay fast for large graphs (target: 20 000 vertices well under 10 s). Reasonable approach: use the exact matrix procedure for small n
(you choose the threshold, ≤ ~800); above it compute the same quantities by other means (single-source shortest paths on the sparse graph,
tree-specific O(n) methods when the graph is a tree). Distances must agree with the matrix procedure to 1e-12 relative; the chosen pair, the path and
the start point must agree whenever the maximum is unique (no other pair within 1e-12 relative). Where ties make the result depend on the relaxation
order, either reproduce it or emit a `warn` entry saying the start point may differ from Fiji — document which in NOTES.md.

## 8. Pruning (opt.prune, opt.pruneEnds) — derive from documentation + black box

Not specified here. The documentation describes: *prune cycle method* (none / shortest branch / lowest intensity voxel / lowest intensity branch — the
last two use the original grey image `opt.orig`) and *prune ends* (remove end-point branches). Find the exact behaviour (what is removed, iteration,
what the outputs and the tagged image become after pruning, interaction with LSP) by experiments with the oracle (prune index 0..3, pruneEnds 0/1,
optional original image), record evidence in NOTES.md, and implement it. Export `PRUNE = true` once all pruning cases pass.

## 9. Oracle and tests

- Oracle: `java -cp /tmp/ij.jar:/tmp/askel_out SkelDriver <skeleton.tif> <prune 0..3> <pruneEnds 0/1> <lsp 0/1> <pw> <ph> <pd> [original.tif]`
  prints `TREES n`, one `T` line per tree (branches junctions endPoints junctionVoxels slabs avgLength triples quadruples maxLength [lsp spx spy spz]),
  one `E` line per edge (tree length v1x v1y v1z v2x v2y v2z lengthRa nSlabs), `TAG idx:value …`, and with lsp `LSPIMG idx …` (voxels = 96).
  `-list file` runs many cases (tab-separated argument lines). Write TIFFs with `IO.writeTIFF` from `src/img_io.js` (see the test).
- Harness: `SKEL_IMPL=cleanroom/askel_cr.js node test/test_skeleton.js` (50 fixed cases), plus `FUZZ=<n> SEED=<s>` for random cases.
  Required: 0 failures on the fixed cases and on FUZZ ≥ 3000 with several seeds; with PRUNE, the pruning cases too.
- Performance: linear-time traversals (no quadratic queues/arrays). Targets: 1500×1500 fully foreground (all junction voxels) < 5 s; a 2000×2000
  serpentine single path with LSP < 3 s; 3D 256³ sparse random skeleton < 10 s.

## 10. Deliverable

`cleanroom/askel_cr.js`: a single self-contained file, no dependencies, usable both as a browser/worker script and from Node:
```
(function(root){ "use strict"; ... function analyzeSkeleton(vol,cal,opt){...}
  const api={analyzeSkeleton, PRUNE:<bool>, TAGS:{END:30,JUNC:70,SLAB:127,SPATH:96}};
  if(typeof module!=="undefined"&&module.exports) module.exports=api;
  if(root&&root.IMG_CORE){ root.IMG_CORE.analyzeSkeleton=analyzeSkeleton; root.IMG_CORE.SKELETON_TAGS=api.TAGS; root.IMG_CORE.SKELETON_PRUNE=api.PRUNE; }
})(typeof self!=="undefined"?self:globalThis);
```
Header comment: state that it is an independent implementation from this behavioural specification and black-box testing, and cite the paper.
Plus `cleanroom/NOTES.md` (evidence log, decisions, any remaining known differences).
