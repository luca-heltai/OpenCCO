#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Check requested/actual terminal counts, finite budgets and seeded replay."""
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import xml.etree.ElementTree as ET


BINARIES = [Path(arg).resolve() for arg in sys.argv[1:3]]
del sys.argv[1:3]


class GenerationCLITest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.sequence = 0

    def run_tree(self, binary, terminals, extra=()):
        self.sequence += 1
        out = self.root / str(self.sequence)
        out.mkdir()
        with (out / "run.log").open("w") as log:
            result = subprocess.run([str(binary), "-n", str(terminals), "--seed", "23",
                                     "-x", "tree.xml", *extra], cwd=out,
                                    stdout=log, stderr=subprocess.STDOUT, timeout=30)
        return result.returncode, out

    def actual_terminals(self, path):
        document = ET.parse(path)
        nodes, edges = document.findall(".//node"), document.findall(".//edge")
        ids = {node.attrib["id"] for node in nodes}
        self.assertEqual(len(ids), len(nodes))
        self.assertEqual(len(edges), len(nodes) - 1)
        adjacency = {node: set() for node in ids}
        for edge in edges:
            a, b = edge.attrib["from"], edge.attrib["to"]
            self.assertIn(a, ids)
            self.assertIn(b, ids)
            self.assertNotEqual(a, b)
            self.assertNotIn(b, adjacency[a])
            adjacency[a].add(b)
            adjacency[b].add(a)
        seen, pending = {"n0"}, ["n0"]
        for node in pending:
            for other in adjacency[node] - seen:
                seen.add(other)
                pending.append(other)
        self.assertEqual(seen, ids)
        self.assertEqual(len(adjacency["n0"]), 1)
        self.assertTrue(all(len(neighbors) in (1, 3) for neighbors in adjacency.values()))
        self.assertTrue(all(math.isfinite(float(value.text))
                            for value in document.findall(".//float")))
        for edge in edges:
            radii = [attr.find("float") for attr in edge.findall("attr")
                     if attr.attrib["name"].strip() == "radius"]
            self.assertEqual(len(radii), 1)
            self.assertGreater(float(radii[0].text), 0.)
        return sum(len(neighbors) == 1 for node, neighbors in adjacency.items() if node != "n0")

    def test_full_counts_and_exact_seeded_replay(self):
        for binary in BINARIES:
            for extra in [[], ["-s"]]:
                status, out = self.run_tree(binary, 8, extra)
                self.assertEqual(status, 0, (out / "run.log").read_text()[-2000:])
                self.assertEqual(self.actual_terminals(out / "tree.xml"), 8)
                again, replay = self.run_tree(binary, 8, extra)
                self.assertEqual(again, 0, (replay / "run.log").read_text()[-2000:])
                self.assertEqual((out / "tree.xml").read_bytes(), (replay / "tree.xml").read_bytes())

    def test_last_allowed_search_counts_as_success(self):
        for binary in BINARIES:
            status, out = self.run_tree(binary, 2, ["--searchTrials", "1"])
            self.assertEqual(status, 0, (out / "run.log").read_text()[-2000:])
            self.assertEqual(self.actual_terminals(out / "tree.xml"), 2)

    def test_partial_tree_is_not_success(self):
        for binary in BINARIES:
            status, out = self.run_tree(binary, 40,
                                        ["--searchTrials", "1", "--samplingTrials", "1"])
            self.assertEqual(status, 1)
            self.assertLess(self.actual_terminals(out / "tree.xml"), 40)

    def test_invalid_budgets_and_counts(self):
        for binary in BINARIES:
            for extra in [["--searchTrials", "0"], ["--samplingTrials", "-1"],
                          ["-a", "0"], ["-g", "-1"], ["-m", "-1"]]:
                status, out = self.run_tree(binary, 8, extra)
                self.assertEqual(status, 1)
                self.assertFalse((out / "tree.xml").exists())
            status, out = self.run_tree(binary, 0)
            self.assertEqual(status, 1)
            self.assertFalse((out / "tree.xml").exists())


if __name__ == "__main__":
    unittest.main()
