Clipper package: /Users/robinthibaut/PycharmProjects/loop-cgal

• - P1 ‑ Preserve flags ignored: TriMesh::cutWithSurface accepts preserve_intersection,
    preserve_intersection_clipper, and use_exact_kernel but never uses them; it always calls PMP::clip
    on an inexact Surface_mesh and drops the intersection curve (mesh.cpp:300‑365). In fault_network.py
    these flags are set to True (fault_network.py:330‑334,394‑398,434) and the validator assumes the
    shared intersection edge remains. When two faults overlap, the intersection ring is deleted instead
    of constrained, leaving open borders that later remeshing or extract_feature_edges can mis-handle or
    crash on.
  - P1 ‑ Intersection miss → silent no‑op: cutWithSurface only clips when PMP::do_intersect returns
    true (mesh.cpp:341‑366). With the inexact kernel and open sheets (not closed volumes), near‑coplanar
    or just‑touching faults often return false, so the downstream fault is left unclipped while
    add_abutting_fault is still recorded (fault_network.py:423‑436). The pipeline then believes clipping
    occurred, but geometry and topology remain inconsistent.
  - P1 ‑ Orientation may flip randomly: The positive/negative side decision is based on PyVista implicit
    distance against an open TriMesh surface (fault_network.py:427‑433). Neither TriMesh creation nor
    cutWithSurface enforces consistent normals after remeshing or clipping (mesh.cpp:268‑288,300‑365).
    A slight numerical wobble can reverse the sign, causing the entire downstream surface to be removed
    instead of the intended portion, yielding empty or inverted faults.
  - P2 ‑ Open meshes fed to PMP::clip: CGAL’s clip is designed for closed, consistently oriented meshes;
    here both target and clipper are open sheets (fault panels and planes). With problematic borders, clip
    may return false or create self‑intersections even after ensure_valid_mesh repairs, leaving malformed
    meshes that later steps treat as valid because failures are logged only when LoopCGAL::verbose is on
    (mesh.cpp:345‑353).
  - P2 ‑ Trace‑trimming heuristic fragile: _trim_surfaces_to_trace assumes the highest‑Z boundary polyline
    is the active trace and trims laterally based on analytic endpoints (fault_network.py:188‑283). After a
    complex clip the true trace can be internal or split into multiple components, so the fallback may pick
    an arbitrary boundary loop and cut through the middle of the fault, producing tiny or empty meshes and
    derailing subsequent clipping/validation.