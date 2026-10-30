"""Core tests use synthetic image/feature adapters; identity checks use UNI data.

These tests do not claim to validate the live Hibou checkpoint or GPU adapter.
"""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from threadpoolctl import threadpool_limits

from audit_assets import PACKAGE, verify_inputs
from extract_hibou import stream_slide_mean, pool_patients
from run_encoder_benchmark import fit_site, feature_input, compare_group
import run_controls as controls


class SyntheticSlide:
    dimensions = (2048, 2048)
    def read_region(self, position, level, size):
        if level != 0: raise AssertionError("Changed coordinate level")
        return Image.new("RGB", size, (position[0], position[1], 10))


class PipelineTests(unittest.TestCase):
    def test_slide_mean_streaming_and_batch_boundary_parity(self):
        rows = pd.DataFrame({"patient_id": ["p"] * 5, "file_name": ["a"] * 5,
            "patch_id": range(5), "x_level0": [1, 2, 3, 4, 5], "y_level0": [10] * 5,
            "read_size_level0": [300] * 5})
        def fake(images):
            return np.array([[im.getpixel((0, 0))[0], 2] for im in images], dtype=np.float32)
        a = stream_slide_mean(SyntheticSlide(), rows, fake, 2, 2)
        b = stream_slide_mean(SyntheticSlide(), rows, fake, 4, 2)
        np.testing.assert_array_equal(a, [3., 2.])
        np.testing.assert_array_equal(a, b)
        bad = rows.copy(); bad.loc[0, "x_level0"] = 2000
        with self.assertRaisesRegex(ValueError, "boundary"):
            stream_slide_mean(SyntheticSlide(), bad, fake, 2, 2)

    def test_patient_pooling_is_equal_slide_weight_and_preserves_order(self):
        means = np.array([[1., 2.], [5., 6.], [10., 20.]], dtype=np.float32)
        slide_table = pd.DataFrame({"patient_id": ["p1", "p1", "p2"], "patches": [1000, 1, 4]})
        patient = pd.DataFrame({"patient_id": ["p2", "p1"]})
        result = pool_patients(means, slide_table, patient)
        np.testing.assert_array_equal(result, [[10., 20.], [3., 4.]])

    def test_pilot_not_accepted_as_complete_feature_freeze(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "feature_freeze.json").write_text(json.dumps({"status": "PILOT_ONLY", "model": "histai/hibou-b"}))
            with self.assertRaisesRegex(ValueError, "complete real"):
                feature_input(root)

    def test_identical_encoder_paired_effects_are_exactly_zero(self):
        frame = pd.DataFrame({"patient_id": ["p1", "p2", "p3", "p4"], "label": [0, 1, 0, 1],
                              "p_0": [.8, .2, .7, .1], "p_1": [.2, .8, .3, .9]})
        _, effects = compare_group(frame.iloc[::-1], frame, 100, {"endpoint": "test"})
        for row in effects:
            self.assertEqual(row["delta_hibou_minus_uni"], 0.)
            self.assertEqual(row["ci_lower_95"], 0.)
            self.assertEqual(row["ci_upper_95"], 0.)

    def test_internal_and_site_identity_for_every_endpoint(self):
        # Run the new fitting adapters on UNI features, solely to verify recipe parity.
        # No output from this identity test is labelled as a real second encoder.
        verify_inputs()
        data = PACKAGE / "data"
        X = np.load(data / "stage10a_uni_patient_features_v1.npy", allow_pickle=False)
        internal = pd.read_csv(data / "stage10a_internal_5fold_manifest_v1.csv")
        sites = pd.read_csv(data / "stage10a_site_holdout_manifest_v1.csv")
        old_internal = pd.read_csv(data / "uni_control_predictions.csv", float_precision="round_trip")
        old_site = pd.read_csv(data / "stage10c_site_holdout_predictions_v1.csv", float_precision="round_trip")
        with threadpool_limits(limits=1):
            for endpoint in controls.ENDPOINTS:
                d = internal[internal.endpoint.eq(endpoint)]
                train, test = d[d.internal_fold.ne(0)], d[d.internal_fold.eq(0)]
                prediction, _, _ = controls.fit_control(X[train.feature_row], train.label.to_numpy(int),
                    X[test.feature_row], "std_balanced", controls.SEED)
                for state in prediction:
                    old = old_internal[old_internal.endpoint.eq(endpoint) & old_internal.variant.eq("std_balanced")
                                       & old_internal.state.eq(state)].set_index("patient_id").loc[test.patient_id]
                    expected = old[[f"p_{c}" for c in range(d.label.nunique())]].to_numpy()
                    np.testing.assert_allclose(prediction[state], expected, rtol=0, atol=5e-5)
                d = sites[sites.endpoint.eq(endpoint) & sites.holdout_site.eq("A2")]
                train, val, test = [d[d.partition.eq(x)] for x in ["train", "validation", "test"]]
                prediction, _ = fit_site(X[train.feature_row], train.label.to_numpy(int), X[val.feature_row],
                                        val.label.to_numpy(int), X[test.feature_row])
                old = old_site[old_site.endpoint.eq(endpoint) & old_site.holdout_site.eq("A2")].set_index("patient_id").loc[test.patient_id]
                for state, prefix in [("full_raw", "prob_raw_class_"), ("validation_calibrated", "prob_calibrated_class_")]:
                    expected = old[[f"{prefix}{c}" for c in range(d.label.nunique())]].to_numpy()
                    np.testing.assert_allclose(prediction[state], expected, rtol=0, atol=5e-5)


if __name__ == "__main__": unittest.main(verbosity=2)
