#!/usr/bin/env python3
"""Small standalone tests; all writes are restricted to TemporaryDirectory."""
import copy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

import stage_group as stage


class StageGroupTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="star-recovery-stage-")
        self.addCleanup(self.temp.cleanup)
        root = Path(self.temp.name)
        self.cache, self.scratch = root / "cache", root / "scratch"
        self.cache.mkdir()
        self.scratch.mkdir()
        self.manifest = self.scratch / "group_00012.json"
        self.tokens = self.scratch / "tokens.list"
        self.output = self.scratch / "production_input.list"
        self.data = dict(schemaVersion=1, sourceJobid=stage.SOURCE_JOBID,
                         originalSubjob=12, inputs=[])
        for i in range(3):
            name = "st_physics_22121033_raw_300000{}.picoDst.root".format(i)
            blob = b"root" + bytes([i]) * 50
            path = self.cache / name
            path.write_bytes(blob)
            uri = "root://node//production_3p85GeV_fixedTarget_2021/P24iy/2021/121/" + name
            self.data["inputs"].append(dict(path=str(path), sha256=hashlib.sha256(blob).hexdigest(),
                                            entries=10 + i, uri=uri))

    def run_stage(self, data=None):
        self.manifest.write_text(json.dumps(self.data if data is None else data))
        self.tokens.write_text(str(self.manifest) + "\n")
        return stage.stage_group(self.tokens, self.scratch, self.cache, self.output)

    def test_order_and_checksums(self):
        report = self.run_stage()
        self.assertEqual(report["expectedEntries"], 33)
        paths = [Path(p) for p in self.output.read_text().splitlines()]
        self.assertEqual([p.name for p in paths], [Path(i["path"]).name for i in self.data["inputs"]])
        self.assertEqual([stage.digest(p) for p in paths], [i["sha256"] for i in self.data["inputs"]])

    def test_empty_recorded_and_omitted(self):
        self.data["inputs"][1]["entries"] = 0
        report = self.run_stage()
        self.assertEqual(report["expectedEntries"], 22)
        self.assertEqual(len(report["stagedInputs"]), 2)
        self.assertEqual(report["omittedEmptyInputs"][0]["position"], 1)

    def test_hash_mismatch(self):
        self.data["inputs"][0]["sha256"] = "0" * 64
        with self.assertRaises(ValueError):
            self.run_stage()
        self.assertFalse(self.output.exists())

    def test_empty_hash_is_checked(self):
        self.data["inputs"][0].update(entries=0, sha256="0" * 64)
        with self.assertRaises(ValueError):
            self.run_stage()
        self.assertFalse(self.output.exists())

    def test_all_empty_fails(self):
        for entry in self.data["inputs"]:
            entry["entries"] = 0
        with self.assertRaises(ValueError):
            self.run_stage()

    def test_wrong_dataset_fails(self):
        self.data["inputs"][0]["uri"] = self.data["inputs"][0]["uri"].replace("2021", "2019")
        with self.assertRaises(ValueError):
            self.run_stage()

    def test_duplicate_fails(self):
        self.data["inputs"][1] = copy.deepcopy(self.data["inputs"][0])
        with self.assertRaises(ValueError):
            self.run_stage()

    def test_cache_symlink_fails(self):
        original = Path(self.data["inputs"][0]["path"])
        link = self.cache / "link.root"
        link.symlink_to(original)
        self.data["inputs"][0]["path"] = str(link)
        with self.assertRaises(ValueError):
            self.run_stage()

    def test_no_clobber(self):
        self.output.write_text("existing\n")
        with self.assertRaises(ValueError):
            self.run_stage()
        self.assertEqual(self.output.read_text(), "existing\n")


if __name__ == "__main__":
    unittest.main()
