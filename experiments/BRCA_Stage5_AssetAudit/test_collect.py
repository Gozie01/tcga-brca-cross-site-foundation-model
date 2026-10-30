import importlib.util,json,subprocess,sys,tempfile,unittest,zipfile
from pathlib import Path
from collect import wanted

class Collection(unittest.TestCase):
    def test_selection(self):
        self.assertTrue(wanted(Path('stage10c/ER_binary_A2_model.joblib')))
        self.assertTrue(wanted(Path('stage10a_uni_patient_features_v1.npy')))
        self.assertTrue(wanted(Path('cptac_patient_clinical_data.json')))
        self.assertFalse(wanted(Path('patch_embeddings.npy')))
        self.assertFalse(wanted(Path('patient_slide.svs')))
        self.assertFalse(wanted(Path('credentials.json')))

    def test_no_unpickling_and_no_unrelated_collection(self):
        with tempfile.TemporaryDirectory() as d:
            base=Path(d);root=base/'TCGA_BRCA_Mutation_Audit';root.mkdir()
            (root/'endpoint_model.joblib').write_bytes(b'not a pickle; copied as bytes')
            (root/'stage10a_uni_patient_features_v1.npy').write_bytes(b'feature fixture')
            other=base/'unrelated';other.mkdir();(other/'model.joblib').write_bytes(b'private')
            # Run a copy so test reports stay inside the temporary directory.
            script=base/'collect.py';script.write_bytes(Path(__file__).with_name('collect.py').read_bytes())
            subprocess.run([sys.executable,str(script),'--project-root',str(base)],check=True,capture_output=True)
            with zipfile.ZipFile(next(base.glob('brca_asset_audit_*.zip'))) as z:
                r=json.loads(z.read('asset_inventory.json'));self.assertEqual(r['model_file_count'],1);self.assertEqual(r['file_count'],2)
                self.assertFalse(r['models_loaded']);self.assertEqual(z.read('assets/TCGA_BRCA_Mutation_Audit/endpoint_model.joblib'),b'not a pickle; copied as bytes')

if __name__=='__main__':unittest.main()
