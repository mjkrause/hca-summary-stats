#!/usr/bin/env python3

import unittest
import loompy
from src.matrix_summary_stats import MatrixSummaryStats

class TestMatrixSummaryStats(unittest.TestCase):

    def setUp(self) -> None:
        self.loomfile = 'test/data/58eefffe-f0af-490b-ae49-bce09602b8e6.loom'

    def test_gene_max_count(self):
        pass