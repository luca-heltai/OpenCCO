#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Exercise the Python/native anatomy boundary without downloading patient data."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


sys.dont_write_bytecode = True
INSPECTOR = Path(sys.argv.pop(1)).resolve()
MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts/prepare_liver_inputs.py"
spec = importlib.util.spec_from_file_location("liver_inputs", MODULE_PATH)
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)


class AnatomyBoundaryTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.nodes = [(0., 0., 0., .0012345678901234567, 0, 0),
                      (0., .02, 0., .002, 1, 0),
                      (-.02, -.01, 0., .0005, 0, 1),
                      (.02, -.01, 0., .0006, 0, 0)]
        self.edges = [(1, 0), (0, 2), (0, 3)]

    def test_native_roundtrip_all_roles_and_unrelated_cwd(self):
        for role in adapter.ROLES:
            path = self.root / f"{role}.anatomy"
            out = self.root / f"{role}.roundtrip.anatomy"
            adapter.write_records(path, role, self.nodes, self.edges)
            result = subprocess.run([str(INSPECTOR), str(path), str(out)],
                                    cwd=self.root, check=True, capture_output=True, text=True)
            report = json.loads(result.stdout)
            self.assertEqual(report, dict(role=role, units="metres", nodes=4, edges=3,
                                          root=1, growth_endpoints=1,
                                          connected_tree=True, growth_performed=False))
            self.assertEqual(adapter.read_records(out), (role, self.nodes, self.edges))
            original = path.read_bytes()
            # Both the Python exporter and the native inspector refuse overwrite.
            with self.assertRaises(FileExistsError):
                adapter.write_records(path, role, self.nodes, self.edges)
            result = subprocess.run([str(INSPECTOR), str(path), str(path)],
                                    capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(path.read_bytes(), original)

    def test_invalid_radius_fails_in_native_inspector(self):
        for value in [-1., 0., float("inf"), float("nan")]:
            path = self.root / f"bad-{value}.anatomy"
            nodes = [(*self.nodes[0][:3], value, 0, 0), *self.nodes[1:]]
            adapter.write_records(path, "arterial", nodes, self.edges)
            result = subprocess.run([str(INSPECTOR), str(path)], capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(result.stdout)

    def test_truncated_roundtrip_is_rejected(self):
        path = self.root / "truncated.anatomy"
        adapter.write_records(path, "portal", self.nodes, self.edges)
        path.write_text(path.read_text().rsplit("\n", 2)[0] + "\n")
        with self.assertRaises(ValueError):
            adapter.read_records(path)


if __name__ == "__main__":
    unittest.main()
