#!/usr/bin/env python3
"""Synthetic, no-ROOT tests for the saved-candidate purity study.

env -u PYTHONPATH /usr/bin/python3 tests/test_lambda_purity_scan.py
"""
import copy
import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import yaml

REPO = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("purity_scan", REPO / "common/macro/studyLambdaKfPurity.py")
SCAN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SCAN)
CONFIG = REPO / "config/qa/lambda_kf_purity_study_20260927.yaml"


def fixture(folder):
    """Independent simple spectra, with explicit invalid-final-cut rows."""
    c = SCAN.load_config(CONFIG)
    c.update(inputEvents=100, reconstructedEvents=100, bootstrapReplicates=200,
             minWindowCountForSelection=10, minSidebandCountForSelection=1)
    c["cuts"] = {"nSigmaProton": copy.deepcopy(c["cuts"]["nSigmaProton"])}
    c["cuts"]["nSigmaProton"]["values"] = [2.5]
    c["tripleFamilies"] = []
    c["eventIdsFile"] = str(folder / "events.tsv")
    c["signalRetentionFloors"] = [.95, .8]
    (folder / "events.tsv").write_text("".join("{}\t1\t{}\n".format(i, i+1) for i in range(100)))
    cfg = folder / "config.yaml"
    cfg.write_text(yaml.safe_dump(c, sort_keys=False))
    export = folder / "export"
    export.mkdir()
    effective = {
        "schemaVersion": 2, "backend": "finder_topo", "selectionProfile": "lambda_imp5",
        "pidProfile": "pico_nsigma", "nSigmaProton": 3., "nSigmaPion": 3.,
        "imp5MinDCAProton": .7, "imp5MinDCAPion": 1., "imp5MaxPathLength": 100.,
        "minMass": 1.05, "maxMass": 1.25, "maxMassError": -1., "maxChi2Ndf": -1.,
        "maxTopoChi2Ndf": -1., "maxDaughterDistance": .5, "maxDistanceToPv": .5,
        "minDecayLength": -1., "minDecayLengthSignificance": -1., "minVertexLineSignificance": -1.,
        "minCosPointing": .998, "useTof": False, "strictTofPid": False, "cleanKaonsWithTof": False,
        "useHftTracksOnly": False, "reconstructAntiLambda": False, "applyLambdaGeometryCuts": True,
    }
    prototype = {
        "runId": 1, "eventId": 1, "protonId": 1, "pionId": 2, "selected": 1,
        "mass": 1.1155, "massError": .001, "pdg": 3122, "imp5PathValid": 1,
        "protonHelixPathLength": 20., "pionHelixPathLength": 20., "protonHasTof": 0, "pionHasTof": 0,
        "protonPidPull": 2., "pionPidPull": 2., "protonDcaToPv": 2., "pionDcaToPv": 2.,
        "daughterDistance": .1, "distanceToPv": .1, "cosPointing": .9999, "chi2Ndf": .5,
        "topoChi2Ndf": .5, "decayLength": 10., "decayLengthSignificance": 20.,
        "vertexLineLengthSignificance": 20.,
    }
    rows = []
    # Last ten input events have no candidate at all; they must still appear
    # in the bootstrap event axis and split counters.
    for event in range(1, 91):
        for candidate in range(11):
            row = dict(prototype, eventId=event, protonId=candidate*2+1, pionId=candidate*2+2)
            if candidate in (8, 9):
                row["mass"] = 1.1025 if candidate == 8 else 1.1305
                row["protonPidPull"] = 3. if event % 2 else 2.
            if candidate == 10:
                row["selected"] = 0; row["daughterDistance"] = .8
            rows.append(row)
    SCAN.write_csv(export / "candidates.csv", list(prototype), rows)
    masses = np.array([row["mass"] for row in rows if row["selected"]])
    counts = np.bincount(SCAN.mass_bins(masses, c), minlength=202)
    native = [{"bin": i, "low": 1.05+(i-1)*.001, "high": 1.05+i*.001,
               "count": int(counts[i]), "error": float(np.sqrt(counts[i]))} for i in range(202)]
    SCAN.write_csv(export / "selected_histogram.csv", list(native[0]), native)
    head = {"inputRoot": "synthetic.root", "inputEvents": 100, "reconstructedEvents": 100, "rawEntries": len(rows)}
    (export / "metadata.txt").write_text(yaml.safe_dump(head)+"\nKFParticleEffectiveConfiguration\n"+
        yaml.safe_dump(effective)+"\nKFEventSelectionConfiguration\nKFEventSelection.policy: legacy_lambda_imp5\n")
    datasets = folder / "datasets.tsv"
    datasets.write_text("id\texport_dir\nbaseline\texport\n")
    return c, cfg, export, datasets


class LambdaPurityScanTest(unittest.TestCase):
    def test_grid_bounded_unique(self):
        c = SCAN.load_config(CONFIG)
        points = SCAN.enumerate_points(c)
        self.assertEqual(len(points), 1020)
        self.assertEqual(len(points), len(set(points)))
        self.assertEqual(points[0], ())
        self.assertLessEqual(max(map(len, points)), 3)

    def test_event_hash_is_candidate_independent(self):
        c = SCAN.load_config(CONFIG)
        values = [SCAN.event_train(1, event, c["splitSeed"], c["trainFraction"]) for event in range(100)]
        self.assertTrue(any(values)); self.assertFalse(all(values))
        self.assertEqual(values, [SCAN.event_train(1, event, c["splitSeed"], c["trainFraction"]) for event in range(100)])

    def test_inclusive_extra_cut_boundaries(self):
        c = {"cuts": {"nSigmaProton": {"field": "protonPidPull", "operator": "abs_le", "values": [2.]},
                      "minDecayLength": {"field": "decayLength", "operator": "ge", "values": [3.]}}}
        d = {"data": {"protonPidPull": np.array([-2., 2., 2.001]), "decayLength": np.array([3., 2.99, 3.])},
             "selected": np.array([True, True, True])}
        d["cut_masks"] = SCAN.cut_masks(d, c)
        mask = SCAN.point_mask(d, (("nSigmaProton", 2.), ("minDecayLength", 3.)))
        self.assertEqual(mask.tolist(), [True, False, False])

    def test_purity_definitions_and_no_clipping(self):
        result = SCAN.measure(100, 10, 10, .75)
        self.assertEqual(result["B"], 15)
        self.assertEqual(result["S"], 85)
        self.assertEqual(result["purity"], .85)
        self.assertAlmostEqual(result["SoverB"], 85/15)
        self.assertLess(SCAN.measure(1, 10, 10, .75)["purity"], 0)
        self.assertTrue(np.isnan(SCAN.measure(0, 0, 0, .75)["purity"]))

    def test_low_sideband_point_cannot_win(self):
        c = SCAN.load_config(CONFIG)
        base = {"signalRetention": 1., "S": 100., "Nwin": 120., "Nleft": 10., "Nright": 10.,
                "score": .8, "changedParameters": 0, "dataset_id": "baseline", "point_id": 0}
        zero = dict(base, Nleft=0., Nright=0., score=1., point_id=1)
        self.assertTrue(all(row["point_id"] == 0 for floor, row in SCAN.choose_points([base, zero], c)))

    def test_export_replay_and_fail_closed_bits(self):
        with tempfile.TemporaryDirectory(prefix="lambda-purity-tests-") as temp:
            c, cfg, export, datasets = fixture(Path(temp))
            data = SCAN.load_dataset("baseline", export, c)
            self.assertEqual(int(data["selected"].sum()), 900)
            table = SCAN.read_numeric_csv(export / "candidates.csv")
            table["selected"][0] = 0
            rows = [dict(zip(table, values)) for values in zip(*table.values())]
            SCAN.write_csv(export / "candidates.csv", list(table), rows)
            with self.assertRaisesRegex(ValueError, "selected bits"):
                SCAN.load_dataset("baseline", export, c)

    def test_nonfinite_csv_rejected(self):
        with tempfile.TemporaryDirectory(prefix="lambda-purity-tests-") as temp:
            path = Path(temp) / "bad.csv"
            path.write_text("a,b\n1,nan\n")
            with self.assertRaisesRegex(ValueError, "nonfinite"):
                SCAN.read_numeric_csv(path)

    def test_complete_synthetic_integration_and_no_overwrite(self):
        with tempfile.TemporaryDirectory(prefix="lambda-purity-tests-") as temp:
            folder = Path(temp)
            c, cfg, export, datasets = fixture(folder)
            output = folder / "result"
            summary = SCAN.run(cfg, datasets, output)
            self.assertEqual(summary["pointsPerDataset"], 2)
            self.assertEqual(sum(summary["inputEventSplitCounts"].values()), 100)
            self.assertTrue(json.loads((output / "status.json").read_text())["complete"])
            events = np.load(output / "chosen_event_counts.npz")
            self.assertEqual(events["event_ids"].shape, (100, 2))
            self.assertTrue((events["counts"][:, 90:, :] == 0).all())
            choices = json.loads((output / "chosen_points.json").read_text())
            self.assertEqual(choices["selectionSplit"], "train")
            self.assertEqual(choices["selectionSideband"], "default")
            self.assertTrue((output / "chosen_cuts.tsv").read_text().startswith("selection\tkey\toperator\tvalue"))
            with (output / "chosen_histograms.csv").open() as stream:
                hist = list(csv.DictReader(stream))
            for choice in choices["choices"]:
                bins = [int(row["bin"]) for row in hist if row["selection"] == choice["choice_id"]]
                self.assertEqual(bins, list(range(202)))
            with self.assertRaisesRegex(ValueError, "refusing overwrite"):
                SCAN.run(cfg, datasets, output)

    def test_wrong_ordered_event_id_rejected(self):
        with tempfile.TemporaryDirectory(prefix="lambda-purity-tests-") as temp:
            folder = Path(temp)
            c, cfg, export, datasets = fixture(folder)
            (folder / "events.tsv").write_text("0\t1\t999\n")
            with self.assertRaisesRegex(ValueError, "exactly inputEvents"):
                SCAN.run(cfg, datasets, folder / "result")


if __name__ == "__main__":
    unittest.main(verbosity=2)
