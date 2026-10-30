"""Offline tests requiring only Python's standard library."""
import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

import audit_assets as assets
from audit_external_cohort import candidate_attributes


class AssetTests(unittest.TestCase):
    def test_inventory_missing_duplicate_and_overlapping_roots(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "one").mkdir(); (root / "two").mkdir()
            (root / "one/a.svs").write_bytes(b"one")
            (root / "one/b.svs").write_bytes(b"same")
            (root / "two/b.svs").write_bytes(b"other")
            (root / "unrelated.svs").write_bytes(b"ignore")
            found, errors = assets.index_slides([root, root / "one"], {"a.svs", "b.svs", "missing.svs"})
            self.assertFalse(errors)
            self.assertEqual(len(found["a.svs"]), 1)
            self.assertEqual(len(found["b.svs"]), 2)
            self.assertEqual(found["missing.svs"], [])

    def test_input_hash_tampering_and_inventory_zip(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); package = root / "package"; data = package / "data"; data.mkdir(parents=True)
            slides = data / "slide_manifest.csv"
            slides.write_text("source,patient_id,file_name,patch_count\nbatch01,p,a.svs,2\n")
            coordinates = data / "patch_coordinates.csv.gz"; coordinates.write_bytes(b"placeholder")
            lock = {"files": {slides.name: assets.sha(slides), coordinates.name: assets.sha(coordinates)},
                    "patients": 1, "slides": 1, "patches": 2}
            (data / "input_lock.json").write_text(json.dumps(lock))
            images = root / "images"; images.mkdir(); (images / "a.svs").write_bytes(b"test")
            with contextlib.redirect_stdout(io.StringIO()):
                result = assets.audit([images], root / "audit", package, probe_hardware=False)
            self.assertEqual(result["status"], "READY_FOR_PILOT")
            self.assertTrue((root / "audit.zip").is_file())
            coordinates.write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "input changed"):
                assets.verify_inputs(package)

    def test_receptor_candidate_detection_does_not_treat_erbb2_as_er(self):
        attrs = [{"clinicalAttributeId": "ERBB2_STATUS", "displayName": "HER2 status"},
                 {"clinicalAttributeId": "ER_STATUS", "displayName": "ER clinical status"},
                 {"clinicalAttributeId": "PAM50", "displayName": "PAM50"}]
        found = candidate_attributes(attrs)
        self.assertEqual([x["clinicalAttributeId"] for x in found["ER_binary"]], ["ER_STATUS"])
        self.assertEqual([x["clinicalAttributeId"] for x in found["HER2_binary"]], ["ERBB2_STATUS"])


if __name__ == "__main__": unittest.main(verbosity=2)
