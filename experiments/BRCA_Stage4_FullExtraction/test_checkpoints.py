import io,json,tempfile,unittest,zipfile
from pathlib import Path
import numpy as np
from core import identity,commit,checkpoint,pool,read_rows

class Checkpoints(unittest.TestCase):
    def test_roundtrip_rejects_changed_input_and_corruption(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);a=np.ones((3,4),dtype=np.float32)
            commit(root,'a',a,'original',{'file_name':'slide','patient_id':'p','token':'a'})
            self.assertEqual(checkpoint(root,'a','original',4)['rows'],3)
            with self.assertRaises(ValueError):checkpoint(root,'a','changed',4)
            data=bytearray((root/'a.npy').read_bytes());data[-1]^=1;(root/'a.npy').write_bytes(data)
            with self.assertRaises(ValueError):checkpoint(root,'a','original',4)

    def test_partial_output_without_commit_is_not_reused(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);np.save(root/'a.npy',np.ones((2,4),dtype=np.float32))
            self.assertIsNone(checkpoint(root,'a','x',4))

    def test_invalid_values_rejected_even_with_matching_file_hash(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);commit(root,'a',np.full((2,4),np.nan,dtype=np.float32),'x',{'token':'a'})
            with self.assertRaises(ValueError):checkpoint(root,'a','x',4)

    def test_patient_pooling_is_equal_slide_not_patch_weighted(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);np.save(root/'a.npy',np.ones((100,2),dtype=np.float32));np.save(root/'b.npy',np.full((1,2),3,dtype=np.float32))
            rows=[{'token':'a','patient_id':'p','file_name':'a','rows':100},{'token':'b','patient_id':'p','file_name':'b','rows':1}]
            _,_,a,p=pool(root,rows,2)
            np.testing.assert_array_equal(a,np.array([[2,2]],dtype=np.float32));self.assertEqual(p,['p'])

if __name__=='__main__':unittest.main()
