#!/usr/bin/env python3

import os
import unittest
import shutil
from src.matrix_summary_stats import MatrixSummaryStats
from more_itertools import first


class TestMatrixSummaryStats(unittest.TestCase):

    def setUp(self) -> None:
        self.path = os.path.dirname(os.path.abspath(os.curdir))
        self.mtxdir = os.path.join(self.path,
                                   'test/data/5ff7733b-916f-4cca-b466-99b1df250e25.mtx.zip')
                                   #'test/data/04b7e4ff-a4fd-4a77-9453-16cd0f7ca030.mtx.zip')
        self.mss = MatrixSummaryStats(self.mtxdir)

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

    def test_request_matrix(self):

        status_response = self.mss.request_matrix()
        print(status_response)
