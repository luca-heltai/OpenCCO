#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Convert ImmersX-prepared IRCAD12 VTK trees to audited OpenCCO SI inputs.

Run ImmersX's prepare_ircad_liver.py --with-portal first. This adapter does not
extract centerlines, change growth flags, infer roots, or run CCO growth.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


ROLES = ("arterial", "venous", "portal")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_records(path):
    import vtk

    reader = vtk.vtkUnstructuredGridReader()
    reader.SetFileName(str(path))
    reader.Update()
    grid = reader.GetOutput()
    count = grid.GetNumberOfPoints()
    if reader.GetErrorCode() or count < 2:
        raise ValueError(f"Cannot read nonempty vascular tree: {path}")
    arrays = {}
    for name in ("Radius", "Root", "Grow"):
        array = grid.GetPointData().GetArray(name)
        if (array is None or array.GetNumberOfComponents() != 1 or
                array.GetNumberOfTuples() != count):
            raise ValueError(f"Missing or malformed {name} point array: {path}")
        arrays[name] = [array.GetTuple1(i) for i in range(count)]
    for name in ("Root", "Grow"):
        if any(value not in (0., 1.) for value in arrays[name]):
            raise ValueError(f"Nonbinary {name} flags: {path}")
    nodes = [(*grid.GetPoint(i), arrays["Radius"][i],
              int(arrays["Root"][i]), int(arrays["Grow"][i]))
             for i in range(count)]
    edges = []
    for i in range(grid.GetNumberOfCells()):
        cell = grid.GetCell(i)
        if cell.GetCellType() != vtk.VTK_LINE or cell.GetNumberOfPoints() != 2:
            raise ValueError(f"Expected VTK_LINE cell {i}: {path}")
        edges.append((cell.GetPointId(0), cell.GetPointId(1)))
    return nodes, edges


def write_records(path, role, nodes, edges):
    with path.open("x") as stream:
        stream.write(f"OpenCCOAnatomy 1\n{role} metres\n{len(nodes)} {len(edges)}\n")
        for node in nodes:
            stream.write(" ".join(repr(value) for value in node) + "\n")
        for a, b in edges:
            stream.write(f"{a} {b}\n")


def read_records(path):
    tokens = path.read_text().split()
    if (len(tokens) < 6 or tokens[:2] != ["OpenCCOAnatomy", "1"] or
            tokens[3] != "metres"):
        raise ValueError(f"Invalid native roundtrip header: {path}")
    n, e = map(int, tokens[4:6])
    if len(tokens) != 6 + 6 * n + 2 * e:
        raise ValueError(f"Invalid native roundtrip counts: {path}")
    nodes = [tuple(map(float, tokens[6 + 6*i:12 + 6*i])) for i in range(n)]
    offset = 6 + 6*n
    edges = [tuple(map(int, tokens[offset + 2*i:offset + 2*i + 2])) for i in range(e)]
    return tokens[2], nodes, edges


def prepare(args):
    import vtk

    prepared = args.prepared.resolve()
    inspector = args.inspector.resolve()
    preparation_script = args.preparation_script.resolve()
    inputs = [prepared / f"{role}.vtk" for role in ROLES]
    inputs += [prepared / "liver.stl", prepared / "provenance.json", preparation_script]
    for path in inputs:
        if not path.is_file():
            raise ValueError(f"Missing preparation input: {path}")
    if not inspector.is_file():
        raise ValueError(f"Missing native inspectAnatomy executable: {inspector}")
    provenance = json.loads((prepared / "provenance.json").read_text())
    if provenance.get("case") != "3Dircadb1.12" or provenance.get("output_units") != "m":
        raise ValueError("Expected ImmersX IRCAD12 preparation in metres")
    # Read all roles before creating output. Never silently infer missing flags.
    records = {role: source_records(prepared / f"{role}.vtk") for role in ROLES}
    hashes = {str(path): sha256(path) for path in inputs}
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    manifest = {
        "schema": 1, "status": "input-validation-in-progress", "case": "3Dircadb1.12",
        "units": "metres", "license": provenance.get("license"),
        "source_page": provenance.get("source_page"),
        "source_provenance": provenance, "input_sha256": hashes,
        "adapter_sha256": sha256(Path(__file__).resolve()),
        "inspector_sha256": sha256(inspector),
        "vtk_version": vtk.vtkVersion.GetVTKVersion(),
        "python_version": sys.version, "growth_performed": False, "roles": {},
    }
    manifest_path = output / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    try:
        for role, (nodes, edges) in records.items():
            path = output / f"{role}.anatomy"
            write_records(path, role, nodes, edges)
            roundtrip = output / f"{role}.roundtrip.anatomy"
            result = subprocess.run([str(inspector), str(path), str(roundtrip)],
                                    capture_output=True, text=True)
            if result.returncode:
                raise ValueError(f"Native inspector rejected {role}: {result.stderr.strip()}")
            summary = json.loads(result.stdout)
            actual_role, actual_nodes, actual_edges = read_records(roundtrip)
            if actual_role != role or actual_nodes != nodes or actual_edges != edges:
                raise ValueError(f"Native roundtrip changed {role} input data")
            summary["native_roundtrip_exact"] = True
            manifest["roles"][role] = summary
            shutil.copyfile(prepared / f"{role}.vtk", output / f"{role}.vtk")
        shutil.copyfile(prepared / "liver.stl", output / "liver.stl")
        shutil.copyfile(prepared / "provenance.json", output / "ircad-source-provenance.json")
        # Record the actual preparation script, without touching its repository.
        shutil.copyfile(preparation_script, output / "prepare_ircad_liver.py")
        if any(sha256(Path(path)) != digest for path, digest in hashes.items()):
            raise ValueError("Preparation inputs changed during conversion")
        manifest["output_sha256"] = {
            path.name: sha256(path) for path in sorted(output.iterdir())
            if path.is_file() and path != manifest_path}
        manifest["status"] = "input-validation-passed"
    except Exception as error:
        manifest["status"] = "input-validation-failed"
        manifest["error"] = str(error)
        raise
    finally:
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(manifest_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", type=Path, required=True,
                        help="Output of ImmersX prepare_ircad_liver.py --with-portal")
    parser.add_argument("--inspector", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="New output directory")
    parser.add_argument("--preparation-script", type=Path,
                        default=Path.home() / "c++/immersx/scripts/prepare_ircad_liver.py")
    prepare(parser.parse_args())


if __name__ == "__main__":
    main()
