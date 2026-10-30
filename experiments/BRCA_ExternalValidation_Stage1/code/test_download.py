"""Integrity and interruption tests using an in-process HTTP stand-in."""
import hashlib
from email.message import Message
import io
import json
from pathlib import Path
import tempfile
import unittest

from download_cptac import download_one
from prepare_inventory import specimen_matches, canonical_patient

PAYLOAD=b'II*\x00'+bytes(range(256))*700


class Response(io.BytesIO):
    def __init__(self, content, start, end, total, modified='Mon, 14 Sep 2020 17:57:47 GMT', status=206):
        super().__init__(content);self.status=status;self.headers=Message()
        self.headers['Content-Range']=f'bytes {start}-{end}/{total}'
        self.headers['Last-Modified']=modified


class Server:
    def __init__(self,payload=PAYLOAD):self.payload=payload;self.requests=[];self.bad_range=False;self.modified=None
    def __call__(self,request,timeout=None):
        value=request.get_header('Range');self.requests.append(value)
        start,end=value.removeprefix('bytes=').split('-');start=int(start);end=int(end) if end else len(self.payload)-1
        return Response(self.payload[start:end+1],start+1 if self.bad_range else start,end,len(self.payload),
                        self.modified or 'Mon, 14 Sep 2020 17:57:47 GMT')


class DownloadTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.dest=Path(self.tmp.name)
        self.row={'url':'https://pathdb.cancerimagingarchive.net/example.svs','file_name':'01BR001-test.svs',
                  'patient_id':'01BR001','bytes':str(len(PAYLOAD)),'prefix_sha256':hashlib.sha256(PAYLOAD[:4096]).hexdigest(),
                  'last_modified':'Mon, 14 Sep 2020 17:57:47 GMT','etag':''}
    def tearDown(self):self.tmp.cleanup()
    def run_download(self,server):return download_one(self.row,self.dest,server,lambda _:None)
    def test_full_transfer_and_existing_hash_verification(self):
        first=self.run_download(Server());self.assertEqual(first['sha256'],hashlib.sha256(PAYLOAD).hexdigest())
        self.assertEqual((self.dest/self.row['file_name']).read_bytes(),PAYLOAD)
        self.run_download(Server())
        path=self.dest/self.row['file_name'];path.write_bytes(PAYLOAD[:-1]+b'x')
        with self.assertRaisesRegex(ValueError,'receipt mismatch'):self.run_download(Server())
    def test_resume_from_prior_saved_prefix(self):
        part=self.dest/(self.row['file_name']+'.part');part.write_bytes(PAYLOAD[:17000])
        state={k:self.row.get(k,'') for k in ['url','bytes','prefix_sha256','last_modified','etag']}
        part.with_name(part.name+'.json').write_text(json.dumps(state))
        server=Server();self.run_download(server)
        self.assertIn('bytes=17000-',server.requests)
        self.assertEqual((self.dest/self.row['file_name']).read_bytes(),PAYLOAD)
    def test_remote_object_change_rejected_before_local_mutation(self):
        server=Server(b'II*\x00'+b'x'*(len(PAYLOAD)-4))
        with self.assertRaisesRegex(ValueError,'prefix changed'):self.run_download(server)
        self.assertEqual(list(self.dest.iterdir()),[])
    def test_changed_modification_time_rejected(self):
        server=Server();server.modified='Tue, 15 Sep 2020 17:57:47 GMT'
        with self.assertRaisesRegex(ValueError,'modification time changed'):self.run_download(server)
    def test_bad_resume_range_rejected(self):
        server=Server();n=0
        def open_response(request,timeout=None):
            nonlocal n
            n+=1
            if n>1:server.bad_range=True
            return server(request,timeout)
        with self.assertRaisesRegex(ValueError,'incorrect resume range'):self.run_download(open_response)
        self.assertFalse((self.dest/self.row['file_name']).exists())
    def test_partial_identity_mismatch_is_preserved_and_rejected(self):
        part=self.dest/(self.row['file_name']+'.part');part.write_bytes(PAYLOAD[:8000])
        part.with_name(part.name+'.json').write_text('{}')
        with self.assertRaisesRegex(ValueError,'Partial file identity'):self.run_download(Server())
        self.assertEqual(part.read_bytes(),PAYLOAD[:8000])
    def test_final_file_without_receipt_is_not_overwritten(self):
        path=self.dest/self.row['file_name'];path.write_bytes(b'keep me')
        with self.assertRaisesRegex(ValueError,'no overwrite'):self.run_download(Server())
        self.assertEqual(path.read_bytes(),b'keep me')
    def test_completed_part_can_finish_after_interruption(self):
        part=self.dest/(self.row['file_name']+'.part');part.write_bytes(PAYLOAD)
        part.with_name(part.name+'.json').write_text(json.dumps({k:self.row.get(k,'') for k in ['url','bytes','prefix_sha256','last_modified','etag']}))
        server=Server();self.run_download(server)
        self.assertEqual(server.requests,['bytes=0-4095'])
    def test_mid_transfer_network_failure_resumes_without_duplicate_bytes(self):
        server=Server();interrupted=False
        class Interrupted(Response):
            def read(self,n=-1):
                if self.tell()>=21000:raise OSError('Simulated network interruption')
                return super().read(min(n,7000))
        def open_response(request,timeout=None):
            nonlocal interrupted
            if request.get_header('Range')=='bytes=0-' and not interrupted:
                interrupted=True
                return Interrupted(PAYLOAD,0,len(PAYLOAD)-1,len(PAYLOAD))
            return server(request,timeout)
        self.run_download(open_response)
        self.assertIn('bytes=21000-',server.requests)
        self.assertEqual((self.dest/self.row['file_name']).read_bytes(),PAYLOAD)


class LinkageTests(unittest.TestCase):
    def test_only_expected_patient_prefix_is_normalized(self):
        self.assertEqual(canonical_patient('X01BR001'),'01BR001')
        self.assertEqual(canonical_patient('CPT000814'),'CPT000814')
        with self.assertRaises(ValueError):canonical_patient('XX01BR001')
    def test_exact_label_and_section_suffix(self):
        s={'stripped_label':'433c2eb6-560f-4387-93af-6c2e1a','specimen_id':'x'}
        self.assertEqual(specimen_matches('433c2eb6-560f-4387-93af-6c2e1a',s),'exact')
        self.assertIsNone(specimen_matches('433c2eb6-560f-4387-93af-6c2e1a_2',s))
    def test_abbreviated_alias_cannot_match_wrong_section_or_identifier(self):
        s={'specimen_id':'[43]-6c2e1a_D1_D1_2'}
        self.assertEqual(specimen_matches('433c2eb6-560f-4387-93af-6c2e1a_2',s),'unique_abbreviated_specimen_alias')
        self.assertIsNone(specimen_matches('433c2eb6-560f-4387-93af-6c2e1a_3',s))
        self.assertIsNone(specimen_matches('443c2eb6-560f-4387-93af-6c2e1a_2',s))


if __name__=='__main__':unittest.main(verbosity=2)
