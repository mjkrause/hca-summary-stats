#!/usr/bin/env python3

import os
import json
import requests
import responses
import unittest
import shutil
from src.matrix_summary_stats import MatrixSummaryStats
from more_itertools import first


class TestMatrixSummaryStats(unittest.TestCase):

    def setUp(self) -> None:
        self.path = os.path.dirname(os.path.abspath(os.curdir))
        #self.mtxdir = os.path.join(self.path,
        #                           'test/data/5ff7733b-916f-4cca-b466-99b1df250e25.mtx.zip')
                                   #'test/data/04b7e4ff-a4fd-4a77-9453-16cd0f7ca030.mtx.zip')
        #self.mss = MatrixSummaryStats(self.mtxdir)
        self.mss = MatrixSummaryStats()

    # def tearDown(self) -> None:
    #     shutil.rmtree(os.path.join(self.path, 'test/data/04b7e4ff-a4fd-4a77-9453-16cd0f7ca030.mtx'))
    #
    # # def test_eliminate_dupes(self):
    # #     print(self.mss.eliminate_dupes())
    #
    # def test_unzip_files(self):
    #     self.mss.unzip_files(os.path.join(self.path, 'test/data'))
    #     observed = os.listdir(first(os.path.splitext(self.mtxdir)))
    #     expected = ['barcodes.tsv', 'matrix.mtx', 'genes.tsv']
    #     self.assertEqual(observed, expected)

    @responses.activate
    def test_request_and_unzip_matrix(self):
        def request_callback(request):
            payload = json.loads(request.body)
            resp_body = {'value': payload['numbers']}
            headers = {'request-id': '728d329e-0e86-11e4-a748-0c84dc037c13'}

            return (200, headers, json.dumps(resp_body))

        responses.add_callback(
            responses.POST, 'http://hca.org/matrix',
            callback = request_callback,
            content_type = 'application/json',
        )

        resp = requests.post(
            'http://hca.org/matrix',
            json.dumps({'numbers': [1, 2, 3]}),
            headers = {'content-type': 'application/json'},
        )
        self.assertTrue(resp.json() == {'value': [1, 2, 3]})
        assert len(responses.calls) == 1
        assert responses.calls[0].request.url == 'http://hca.org/matrix'
        assert responses.calls[0].response.text == '{"value": [1, 2, 3]}'
        assert (
            responses.calls[0].response.headers['request-id'] == '728d329e-0e86-11e4-a748-0c84dc037c13'
        )

    # def test_upload_figs_to_s3(self):
    #     os.chdir(self.path)
    #     self.mss.upload_figs_to_s3()

