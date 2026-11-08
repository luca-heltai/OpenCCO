#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

BINARY = Path(sys.argv.pop(1)).resolve()
WITH_SURFACE = "--with-surface" in sys.argv
if WITH_SURFACE:
    sys.argv.remove("--with-surface")


def cube_mask(directory):
    mask = directory / "cube.vol"
    n = 31
    header = ("Center-X: 15\nCenter-Y: 15\nCenter-Z: 15\nX: 31\nY: 31\nZ: 31\n"
              "Voxel-Size: 1\nAlpha-Color: 0\nVoxel-Endian: 0\nInt-Endian: 0123\nVersion: 2\n.\n")
    data = bytes(255 if min(x,y,z) > 0 and max(x,y,z) < n-1 else 0
                 for z in range(n) for y in range(n) for x in range(n))
    mask.write_bytes(header.encode() + data)
    return mask


class AnatomicalGrowthCLI(unittest.TestCase):
    def test_all_roles_preserve_sources_and_materialize_real_growth(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            mask = cube_mask(directory)
            coordinates = {
                "arterial": ((.005, .010, .010), (.010, .010, .010)),
                "portal": ((.010, .025, .015), (.010, .020, .015)),
                "venous": ((.025, .020, .020), (.020, .020, .020)),
            }
            files = []
            for role, (root, tip) in coordinates.items():
                path = directory / f"{role}.anatomy"
                path.write_text(f"OpenCCOAnatomy 1\n{role} metres\n2 1\n" +
                                " ".join(map(repr, root)) + " .0003 1 0\n" +
                                " ".join(map(repr, tip)) + " .0002 0 1\n0 1\n")
                files.append(path)
            original = {p: p.read_bytes() for p in files}
            prior = []
            for role in ["arterial", "portal", "venous"]:
                output = directory / role
                command = [str(BINARY), "--input", str(directory/f"{role}.anatomy"),
                           "--mask", str(mask), "--origin", "0", "0", "0", "--pitch", ".001",
                           "--flow", "1e-8", "--pressure-drop", "400", "--terminals", "8",
                           "--seed", "23", "--output", str(output), "--obstacles", *map(str, files)]
                if WITH_SURFACE and role == "portal":
                    command += ["--surface", str(Path(__file__).resolve().parent/"data/cube_wall.stl")]
                if prior:
                    command += ["--grown-obstacles", *map(str, prior)]
                result = subprocess.run(command, capture_output=True, text=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stderr[-2000:] + result.stdout[-2000:])
                summary = json.loads((output/"summary.json").read_text())
                self.assertEqual(summary["actual"], 8)
                self.assertTrue(summary["role_growth_complete"])
                self.assertFalse(summary["coupled_ready"])
                with (output/"nodes.csv").open() as stream:
                    rows = list(csv.DictReader(stream))
                with (output/"segments.csv").open() as stream:
                    edges = list(csv.DictReader(stream))
                self.assertEqual(sum(int(r["active_tissue"]) for r in rows), 8)
                for i, expected in enumerate(coordinates[role]):
                    self.assertEqual(tuple(float(rows[i][k]) for k in ["x", "y", "z"]), expected)
                    self.assertEqual(int(rows[i]["source_node"]), i)
                self.assertEqual((edges[0]["a"], edges[0]["b"], edges[0]["generated"]), ("0", "1", "0"))
                self.assertEqual((float(edges[0]["r0"]), float(edges[0]["r1"])), (.0003, .0002))
                self.assertEqual(len(edges), len(rows)-1)
                self.assertEqual(len(edges)-1, 15)
                prior.append(output/"segments.csv")
            for path, content in original.items():
                self.assertEqual(path.read_bytes(), content)


    @unittest.skipUnless(WITH_SURFACE, "Original-surface build required")
    def test_seed_inside_stl_but_outside_native_mask_is_not_committed(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            mask = cube_mask(directory)
            anatomy = directory/"near-wall.anatomy"
            anatomy.write_text("OpenCCOAnatomy 1\nvenous metres\n2 1\n"
                               ".009 .0025 .015 .0001 1 0\n"
                               ".010 .0015 .015 .0001 0 1\n0 1\n")
            output = directory/"grown"
            command = [str(BINARY), "--input", str(anatomy), "--mask", str(mask),
                       "--origin", "0", "0", "0", "--pitch", ".001", "--flow", "1e-8",
                       "--pressure-drop", "400", "--terminals", "1", "--seed", "23",
                       "--surface", str(Path(__file__).resolve().parent/"data/cube_wall.stl"),
                       "--match-attachment-radius", "--output", str(output)]
            result = subprocess.run(command, capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr[-2000:] + result.stdout[-2000:])
            with (output/"nodes.csv").open() as stream:
                for node in csv.DictReader(stream):
                    if node["source_node"] == "-1":
                        for coordinate in ["x", "y", "z"]:
                            self.assertGreaterEqual(float(node[coordinate]), .001)
                            self.assertLess(float(node[coordinate]), .030)

    def test_explicit_two_port_allocation(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            mask = cube_mask(directory)
            anatomy = directory / "two-port.anatomy"
            anatomy.write_text("OpenCCOAnatomy 1\nportal metres\n4 3\n"
                               ".005 .015 .015 .0003 1 0\n"
                               ".010 .015 .015 .0003 0 0\n"
                               ".015 .010 .015 .0002 0 1\n"
                               ".015 .020 .015 .0002 0 1\n0 1\n1 2\n1 3\n")
            allocation = directory / "allocation.csv"
            allocation.write_text("node,terminals,flow_fraction\n2,3,.25\n3,5,.75\n")
            output = directory / "grown"
            command = [str(BINARY), "--input", str(anatomy), "--mask", str(mask),
                       "--origin", "0", "0", "0", "--pitch", ".001", "--flow", "1e-8",
                       "--pressure-drop", "400", "--terminals", "8", "--seed", "23",
                       "--allocation", str(allocation), "--output", str(output)]
            result = subprocess.run(command, capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr[-2000:] + result.stdout[-2000:])
            for node, quota, fraction in [(2, 3, .25), (3, 5, .75)]:
                record = json.loads((output/f"attachment-{node}.json").read_text())
                self.assertEqual(record["actual"], quota)
                self.assertEqual(record["design_flow_fraction"], fraction)
            allocation.write_text("node,terminals,flow_fraction\n2,3,.25\n3,5,.5\n")
            command[-1] = str(directory/"invalid")
            result = subprocess.run(command, capture_output=True, text=True, timeout=30)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("unit flow fraction", result.stderr)
            allocation.write_text("node,terminals,flow_fraction\n2,3,.25\n3,5,.75\n")
            command[-1] = str(directory/"radius-matched")
            command += ["--match-attachment-radius", "--local-candidates", "64",
                        "--searches", "1", "--samples", "1"]
            result = subprocess.run(command, capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr[-2000:] + result.stdout[-2000:])
            summary = json.loads((directory/"radius-matched/summary.json").read_text())
            self.assertFalse(summary["pressure_drop_prescribed"])
            self.assertEqual(summary["local_candidate_budget"], 64)
            local_insertions = 0
            for node in [2, 3]:
                record = json.loads((directory/f"radius-matched/attachment-{node}.json").read_text())
                self.assertEqual(record["inlet_radius_m"], .0002)
                self.assertGreater(record["required_pressure_drop_pa"], 0.)
                self.assertNotEqual(record["required_pressure_drop_pa"], 400.)
                local_insertions += record["candidate_counts"].get("local_insertions", 0)
            self.assertGreater(local_insertions, 0)
            initial = directory/"radius-matched"
            allocation.write_text("node,terminals,flow_fraction\n2,4,.25\n3,8,.75\n")
            continued = command.copy()
            continued[continued.index("--terminals")+1] = "12"
            continued[continued.index("--output")+1] = str(directory/"continued")
            continued += ["--bootstrap", str(initial)]
            result = subprocess.run(continued, capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr[-2000:] + result.stdout[-2000:])
            summary = json.loads((directory/"continued/summary.json").read_text())
            self.assertEqual(summary["bootstrap_tip_count"], 8)
            self.assertEqual(summary["actual"], 12)
            self.assertTrue(summary["role_growth_complete"])
            corrupt = directory/"corrupt"
            corrupt.mkdir()
            lines = (initial/"source.anatomy").read_text().splitlines()
            first = lines[3].split()
            first[3] = format(float(first[3])*1.01, ".17g")
            lines[3] = " ".join(first)
            (corrupt/"source.anatomy").write_text("\n".join(lines)+"\n")
            invalid = continued.copy()
            invalid[invalid.index("--bootstrap")+1] = str(corrupt)
            invalid[invalid.index("--output")+1] = str(directory/"invalid-source")
            result = subprocess.run(invalid, capture_output=True, text=True, timeout=30)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("Bootstrap source anatomy changed", result.stderr)
            # The discrete mask cannot supply one million distinct tissue sites.
            # Force a finite partial continuation, retaining BOTH accepted parts.
            allocation.write_text("node,terminals,flow_fraction\n2,500000,.25\n3,500000,.75\n")
            partial = continued.copy()
            partial[partial.index("--terminals")+1] = "1000000"
            partial[partial.index("--output")+1] = str(directory/"partial")
            partial[partial.index("--local-candidates")+1] = "0"
            result = subprocess.run(partial, capture_output=True, text=True, timeout=30)
            self.assertNotEqual(result.returncode, 0)
            summary = json.loads((directory/"partial/summary.json").read_text())
            self.assertFalse(summary["role_growth_complete"])
            self.assertGreaterEqual(summary["actual"], 8)
            for node in [2, 3]:
                self.assertTrue((directory/f"partial/subtrees/attachment-{node}.cco").is_file())
            with (directory/"partial/nodes.csv").open() as stream:
                restored_nodes = list(csv.DictReader(stream))
            for node in [2, 3]:
                self.assertEqual(restored_nodes[node]["source_grow"], "1")
            adjacency = {i: [] for i in range(len(restored_nodes))}
            with (directory/"partial/segments.csv").open() as stream:
                for edge in csv.DictReader(stream):
                    if edge["generated"] == "1":
                        a, b = int(edge["a"]), int(edge["b"])
                        adjacency[a].append(b)
                        adjacency[b].append(a)
            seen, pending = {3}, [3]
            for node in pending:
                for neighbor in adjacency[node]:
                    if neighbor not in seen:
                        seen.add(neighbor)
                        pending.append(neighbor)
            self.assertEqual(sum(restored_nodes[node]["active_tissue"] == "1" for node in seen), 5)


if __name__ == "__main__":
    unittest.main()
