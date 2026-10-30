"""Regression tests for original UNI edge reads and guarded checkpoint migration."""
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from audit_assets import sha
from extract_hibou import stream_slide_mean as legacy_mean
from extract_hibou_edgefix import stream_slide_mean, slide_coordinate_hash, atomic_json
from boundary_migration import LEGACY_CODE_SHA256, fingerprint, edge_summary, migrate_partial


class Slide:
    dimensions = (512, 512)
    calls = []

    def __init__(self, *args):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass

    def read_region(self, position, level, size):
        if level != 0:
            raise AssertionError("Level changed")
        self.calls.append((position, level, size))
        x, y = position
        w, h = size
        yy, xx = np.mgrid[:h, :w]
        valid = (xx + x < 512) & (yy + y < 512)
        a = np.zeros((h, w, 4), dtype=np.uint8)
        a[valid] = [120, 70, 190, 255]
        return Image.fromarray(a)


def rows(xs=(1, 200), ys=(1, 200), size=128):
    return pd.DataFrame({"patient_id": ["p"] * len(xs), "file_name": ["a.svs"] * len(xs),
                         "patch_id": range(len(xs)), "x_level0": xs, "y_level0": ys,
                         "read_size_level0": [size] * len(xs)})


def encode(images):
    return np.array([[np.asarray(im).mean(), np.asarray(im).std(), 1.] for im in images], dtype=np.float32)


class EdgeTests(unittest.TestCase):
    def test_inbounds_matches_legacy_at_multiple_batch_sizes(self):
        d = rows()
        for batch in (1, 2, 8):
            np.testing.assert_array_equal(stream_slide_mean(Slide(), d, encode, batch, 3),
                                          legacy_mean(Slide(), d, encode, batch, 3))

    def test_right_bottom_corner_and_outside_match_original_uni_pixel_reads(self):
        d = rows(xs=(480, 20, 480, 600), ys=(20, 480, 480, 600))
        images = []
        for r in d.itertuples(index=False):
            # Exact Stage9B original read_patch path; no clipping or replacement.
            im = Slide().read_region((r.x_level0, r.y_level0), 0,
                                    (r.read_size_level0, r.read_size_level0)).convert("RGB")
            if im.size != (256, 256):
                im = im.resize((256, 256), Image.Resampling.LANCZOS)
            images.append(im)
        expected = encode(images).astype(np.float64).mean(axis=0).astype(np.float32)
        for batch in (1, 3):
            np.testing.assert_array_equal(stream_slide_mean(Slide(), d, encode, batch, 3), expected)
        summary, detail = edge_summary(Slide(), d, "a.svs")
        self.assertEqual(summary["edge_crossing_patches"], 4)
        self.assertEqual(summary["fully_outside_patches"], 1)
        self.assertEqual(len(detail), 4)

    def test_invalid_coordinates_and_bad_encoder_outputs_still_rejected(self):
        for key, value in (("x_level0", -1), ("y_level0", -1), ("read_size_level0", 0)):
            d = rows(); d.loc[0, key] = value
            with self.assertRaisesRegex(ValueError, "Invalid frozen"):
                stream_slide_mean(Slide(), d, encode, 2, 3)
        with self.assertRaisesRegex(ValueError, "Invalid encoder"):
            stream_slide_mean(Slide(), rows(), lambda images: np.full((len(images), 3), np.nan), 2, 3)


class MigrationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.out = Path(self.temp.name)
        self.checkpoints = self.out / "slide_checkpoints"
        self.checkpoints.mkdir()
        self.d = rows()
        self.template = pd.DataFrame({"patient_id": ["p"], "file_name": ["a.svs"], "patch_count": [2]})
        self.wsi = self.out / "a.svs"; self.wsi.write_bytes(b"source-image-bytes")
        self.paths = {"a.svs": str(self.wsi)}
        self.new = {"version": "brca_revision2_hibou_v2_edgefix", "code_sha256": "new-code",
                    "edge_handling": "full read", "boundary_fix_code_sha256": "migration-code",
                    "model": {"revision": "locked-model"}, "coordinate_sha256": "locked-manifest",
                    "pilot_only": False, "pilot_patches": 0}
        self.old = {k: v for k, v in self.new.items() if k not in {"edge_handling", "boundary_fix_code_sha256"}}
        self.old.update(version="brca_revision2_hibou_v1", code_sha256=LEGACY_CODE_SHA256)
        atomic_json(self.out / "extraction_protocol.json", self.old)
        token = hashlib.sha256(b"a.svs").hexdigest()[:20]
        self.array = self.checkpoints / (token + ".npy")
        np.save(self.array, np.ones(768, dtype=np.float32), allow_pickle=False)
        self.meta_path = self.array.with_suffix(".json")
        self.meta = {"file_name": "a.svs", "patient_id": "p", "patches": 2,
                     "protocol_sha256": fingerprint(self.old), "coordinates_sha256": slide_coordinate_hash(self.d),
                     "mean_sha256": sha(self.array), "wsi_sha256": sha(self.wsi), "elapsed_seconds": 1.}
        atomic_json(self.meta_path, self.meta)

    def tearDown(self):
        self.temp.cleanup()

    def migrate(self, writer=atomic_json):
        migrate_partial(self.out, self.new, self.template, self.d, self.paths,
                        Slide, sha, slide_coordinate_hash, writer)

    def test_migration_retains_vector_bytes_and_backs_up_lineage(self):
        vector = self.array.read_bytes()
        self.migrate()
        self.assertEqual(vector, self.array.read_bytes())
        self.assertEqual(json.loads(self.meta_path.read_text())["protocol_sha256"], fingerprint(self.new))
        self.assertEqual(json.loads((self.out / "boundary_v1_backup" / self.meta_path.name).read_text()), self.meta)
        journal = json.loads((self.out / "boundary_checkpoint_migration.json").read_text())
        self.assertEqual(journal["status"], "COMPLETE")
        self.assertEqual(journal["reused_slide_count"], 1)
        self.migrate()
        self.assertEqual(vector, self.array.read_bytes())

    def test_changed_vector_is_rejected_before_mutation(self):
        np.save(self.array, np.full(768, 2., dtype=np.float32), allow_pickle=False)
        with self.assertRaisesRegex(ValueError, "validation"):
            self.migrate()
        self.assertFalse((self.out / "boundary_v1_backup").exists())
        self.assertEqual(json.loads(self.meta_path.read_text()), self.meta)

    def test_changed_source_image_is_rejected_before_mutation(self):
        self.wsi.write_bytes(b"different")
        with self.assertRaisesRegex(ValueError, "Source slide changed"):
            self.migrate()
        self.assertFalse((self.out / "boundary_v1_backup").exists())

    def test_changed_model_settings_are_rejected(self):
        self.new["model"] = {"revision": "different"}
        with self.assertRaisesRegex(ValueError, "settings differ"):
            self.migrate()
        self.assertFalse((self.out / "boundary_v1_backup").exists())

    def test_inherited_edge_checkpoint_is_rejected(self):
        self.d.loc[0, "x_level0"] = 500
        self.meta["coordinates_sha256"] = slide_coordinate_hash(self.d)
        atomic_json(self.meta_path, self.meta)
        with self.assertRaisesRegex(ValueError, "cannot certify"):
            self.migrate()
        self.assertFalse((self.out / "boundary_v1_backup").exists())

    def test_interrupted_metadata_migration_can_resume(self):
        def interrupted(path, payload):
            if Path(path) == self.meta_path:
                raise RuntimeError("simulated interruption")
            atomic_json(path, payload)
        with self.assertRaisesRegex(RuntimeError, "interruption"):
            self.migrate(interrupted)
        self.assertEqual(json.loads((self.out / "boundary_checkpoint_migration.json").read_text())["status"], "PREPARED")
        self.migrate()
        self.assertEqual(json.loads((self.out / "boundary_checkpoint_migration.json").read_text())["status"], "COMPLETE")

    def test_orphan_checkpoint_is_rejected_before_mutation(self):
        np.save(self.checkpoints / "orphan.npy", np.ones(768, dtype=np.float32))
        with self.assertRaisesRegex(ValueError, "Unpaired"):
            self.migrate()
        self.assertFalse((self.out / "boundary_v1_backup").exists())

    def test_resume_after_checkpoint_rewrite_before_protocol_transition(self):
        def interrupted(path, payload):
            if Path(path) == self.out / "extraction_protocol.json" and payload == self.new:
                raise RuntimeError("simulated protocol interruption")
            atomic_json(path, payload)
        with self.assertRaisesRegex(RuntimeError, "interruption"):
            self.migrate(interrupted)
        self.assertEqual(json.loads(self.meta_path.read_text())["protocol_sha256"], fingerprint(self.new))
        self.assertEqual(json.loads((self.out / "extraction_protocol.json").read_text()), self.old)
        self.migrate()
        self.assertEqual(json.loads((self.out / "boundary_checkpoint_migration.json").read_text())["status"], "COMPLETE")

    def test_resume_after_protocol_transition_before_journal_completion(self):
        def interrupted(path, payload):
            if Path(path) == self.out / "boundary_checkpoint_migration.json" and payload["status"] == "COMPLETE":
                raise RuntimeError("simulated final journal interruption")
            atomic_json(path, payload)
        with self.assertRaisesRegex(RuntimeError, "interruption"):
            self.migrate(interrupted)
        self.assertEqual(json.loads((self.out / "extraction_protocol.json").read_text()), self.new)
        self.migrate()
        self.assertEqual(json.loads((self.out / "boundary_checkpoint_migration.json").read_text())["status"], "COMPLETE")


if __name__ == "__main__":
    unittest.main(verbosity=2)
