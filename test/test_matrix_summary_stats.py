#!/usr/bin/env python3

import os
import unittest
from src.matrix_summary_stats import MatrixSummaryStats

class TestMatrixSummaryStats(unittest.TestCase):

    def setUp(self) -> None:
        path = os.path.dirname(os.path.abspath(os.curdir))
        self.mtxfile = os.path.join(path, 'test/data/mtxfile')
        self.mss = MatrixSummaryStats(self.mtxfile)

    def test_eliminate_dupes(self):
        print(self.mss.eliminate_dupes())
        self.assert