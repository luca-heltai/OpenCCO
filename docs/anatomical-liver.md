# Anatomical liver inputs

Tracking: [fork issue #5](https://github.com/luca-heltai/OpenCCO/issues/5), following
[ImmersX #343](https://github.com/luca-heltai/immersx/issues/343).

The import stage validates real vascular trees without growing them. A separate
experimental growth command is described below; neither stage establishes
hydraulic feasibility or a successful coupled tissue simulation.
The existing `generateTree3D` command still starts from an artificial segment;
it does not consume these anatomical inputs yet.

## Prepare the real dataset

Use the existing ImmersX preparation script, not synthetic feeder stubs:

```sh
PYTHON=/absolute/path/to/python-with-vmtk
IMMERSX="$HOME/c++/immersx"
"$PYTHON" "$IMMERSX/scripts/prepare_ircad_liver.py" \
  "$PWD/build/liver/ircad12-prepared" --with-portal
```

The output directory must be new. The script downloads the IRCAD12 mesh archive,
checks the source hashes, retains the largest intrahepatic vascular components,
converts millimetres to metres and writes `liver.stl`, `arterial.vtk`,
`venous.vtk`, `portal.vtk` and `provenance.json`. Its growth flags select a small
subset of existing leaves; the adapter never silently enables other endpoints.

The tested installed script SHA-256 is
`8008852f1c8a41e1368ca4beb32c1fee68183fab3404a0f44301fef49e4a2b1f`
(ImmersX checkout `272e94b45d0f7ecf106701ff93108c7e59fbdbbe`).
The original tutorial's VMTK environment uses Linux package versions that do not
resolve unchanged on macOS. On the tested Apple Silicon machine, VMTK 1.5.0 was
not available for `osx-arm64`. A project-local `osx-64` environment under Rosetta
worked with Python 3.10.13, VMTK 1.5.0, vtk-base 9.2.6 and NumPy 1.26.4:

```sh
micromamba create -y --platform osx-64 -p "$PWD/build/envs/ircad-prep" \
  -c conda-forge python=3.10.13 vmtk=1.5.0 vtk-base=9.2.6 numpy=1.26.4
PYTHON="$PWD/build/envs/ircad-prep/bin/python"
```

Do not install this into an existing user environment. Native VTK/NumPy suffice
for the conversion stage if the preparation directory already exists; VMTK is
needed only to run the original extraction script.

The data source is [3D-IRCADb-1](https://www.ircad.fr/research-and-development/data-sets/liver-segmentation-3d-ircadb-01/),
case `3Dircadb1.12`, under its stated **CC-BY-NC-ND-4.0** licence. Keep the source
attribution, licence and preparation manifest with data. Patient-derived files
are local build artifacts, not test fixtures or code-PR payloads. This workflow
does not grant permission to redistribute derivatives.

## Convert and audit

Build OpenCCO as described in [the macOS guide](build-macos.md), then run:

```sh
cmake --build build --target inspectAnatomy testAnatomicalTree
"$PYTHON" scripts/prepare_liver_inputs.py \
  --prepared "$PWD/build/liver/ircad12-prepared" \
  --preparation-script "$IMMERSX/scripts/prepare_ircad_liver.py" \
  --inspector "$PWD/build/tools/inspectAnatomy" \
  --output "$PWD/build/liver/validated-inputs"
ctest --test-dir build --output-on-failure --no-tests=error
```

Conversion creates a new directory and refuses to overwrite existing output.
For each role it writes a portable native file, executes the native inspector,
roundtrips the file through C++, and compares every coordinate, radius, flag and
edge with VTK readback using exact double equality. No node merging, radius
averaging, reordering or connectivity reorientation occurs. Original VTK and
STL files are copied unchanged. The output includes the original preparation
script, provenance, source/output SHA-256 hashes, runtime versions and an
explicit `growth_performed: false` flag. Failed conversion retains a failed
manifest rather than a successful result record.

The tested raw source meshes passed all four hashes in the ImmersX script.
Measured prepared inputs:

| Role | Nodes | Line cells | Root index | Enabled growth endpoints |
| --- | ---: | ---: | ---: | --- |
| Arterial | 177 | 176 | 83 | 10 |
| Hepatic venous | 1085 | 1084 | 0 | 460, 489 |
| Portal | 739 | 738 | 0 | 305, 452 |

All three passed native structural validation and exact VTK/native roundtrip.
The prepared surface SHA-256 is
`42491b71fda4b3f475c60cb5d1ce2d75bb44e8b7dd993c6f112c2f7891daa8e1`.
Preparation changes versus other platforms are measured in output hashes; no
cross-platform bitwise reproduction of VMTK extraction is assumed.

## Native input format

`src/AnatomicalTree.h` provides a small physical-data value type and stream
reader/writer. `inspectAnatomy INPUT [NEW_ROUNDTRIP_OUTPUT]` validates a file and
prints a JSON summary. It returns nonzero on invalid input and refuses to
replace an existing roundtrip output.

Version 1 is whitespace-delimited:

```text
OpenCCOAnatomy 1
arterial metres
N E
x y z radius root grow
... N node records ...
a b
... E edge records ...
```

Roles are `arterial`, `venous` and `portal`. Node identifiers are zero-based
record indices. Positions and radii are in metres. `root` and `grow` are binary
flags; there must be exactly one degree-one root, and growth flags may identify
only nonroot degree-one endpoints. The connected graph must have `E=N-1` with
finite coordinates, positive finite radii and positive finite edge lengths.
Malformed/truncated records, invalid indices, trailing data, unsupported versions
and non-metre units are rejected. Allocation is bounded to ten million nodes.
Supplied edge direction is retained but topology is validated as undirected;
root-relative traversal must use the explicit root, not assume node zero.

This representation is separate from `CoronaryArteryTree`: inserting imported
nodes directly into its mutable arrays is **not** an anatomy-preservation API.
Its current optimization updates global length factors, radii and resistances.
Domain queries now use the controller's sampled point type: continuous
controllers retain real-valued SI coordinates, while image masks retain integer
voxel queries. Initialization and `getDomainCenter()` likewise preserve the
configured center. This corrects integer truncation; it does not add a tissue
mesh controller or anatomical immutability.
Protecting fixed geometry and radius traces, constraining attachments, handling
other-role obstacles and exporting generated-vs-imported provenance are still
required before anatomical growth is accepted.

## Experimental growth from imported tips

`growAnatomicalTree` now attaches native CCO subtrees at existing `Grow` endpoints
without optimizing the imported vessels. It reuses `ImageMaskDomainCtrl<3>`;
no new tissue controller is required. Original coordinates, connectivity and
both endpoint-radius traces are retained exactly in the exported prefix.
The new subtree root is fixed at the original tip, not the legacy `--posInit`
location heuristic.

The command requires the `.anatomy` input, existing `.vol` mask, world-coordinate
origin/pitch, physical flow in m³/s and an explicit local pressure-drop budget in
Pa. Pass all three original trees through `--obstacles`, and previously grown
`segments.csv` files through `--grown-obstacles` when growing subsequent roles.
All output directories must be new. The executable writes `nodes.csv`,
`segments.csv`, attachment diagnostics and materialized inspection checkpoints
at 16-terminal increments. These are not restart checkpoints.

Generated-subtree computations use voxel coordinates with
`Q_native = Q_physical / pitch³`, viscosity in Pa·s and physical length factor 1.
The legacy executables retain their length-scaling policy. Imported vessels
are never rescaled. Generated branches are cylindrical native CCO outputs;
imported branches retain their two original endpoint radii. A shared junction
may therefore have different radius traces on its incident branches.

Before committing a candidate, the application checks distinct generated node
positions, nonincident self/source collisions and conservative voxel-domain wall
margins. Max-endpoint-radius capsules provide the initial collision bound.
When these overlap, the existing ImmersX support-plane predicate for flat-ended,
linearly tapered frusta can certify actual separation. The port preserves its
roundoff guards and conservative rejection when no witness is found; it does
not lower a physical clearance threshold. Exact source revision and licensing
are recorded in `src/helpers/FrustumGeometry.h` and `licenses/`.
Only the source edge incident to an attachment is exempted for the directly
attached generated branch, whose initial direction must remain outward.
Source vessels are not checked against each other or silently repaired.

`role_growth_complete` means the requested geometric tip count was reached.
`coupled_ready` remains false: cross-role terminal matching, explicit treatment
of unused anatomical leaves, hydraulic feasibility, original-surface clearance,
archive conversion and actual coupled-loader/solver checks are still required.
Computational sizing flows and local pressure-drop budgets must be recorded;
they are not automatically an adopted physiological configuration.

## Remaining scientific gates

The target from ImmersX #343 is nominally one terminal per global generation
background cell **per role**, with measured per-cell occupancy and at least 95%
coverage within a **maximum tissue-to-terminal distance of `2h`**. This is not a
minimum terminal separation. Cell counts, mesh-size definition and any count
calibration bounds must be recorded with each run.

Small geometry diagnostics, nominal-density three-role growth, preservation of
all imported anatomy, hydraulic/reference feasibility and coupled transient
convergence are separate gates. A binary tree or a successful executable exit
alone establishes none of these. The published ImmersX A4245 arterial result
and its partial companion results are reference evidence, not OpenCCO outputs.
