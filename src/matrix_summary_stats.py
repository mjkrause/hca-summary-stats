#!/usr/bin/env python3

import loompy
import numpy as np


class MatrixSummaryStats:

    def __init__(self, loomfile):
        self.loomfile = loomfile

    def get_max_gene_count(self):
        with loompy.connect(self.loomfile) as ds:
            gene_count = ds.ca['genes_detected']

        pass

