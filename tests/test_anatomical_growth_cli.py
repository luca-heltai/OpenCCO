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


if __name__ == "__main__":
    unittest.main()
