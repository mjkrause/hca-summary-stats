#!/usr/bin/env python3

import os
import csv


class MatrixSummaryStats:

    def __init__(self, mtxdir):
        self.genes = os.path.join(mtxdir, 'genes.tsv')

    def eliminate_dupes(self):
        # TODO: this method can probably be made more efficient
        genes = []
        with open(self.genes, 'r') as fh:
            reader = csv.reader(fh, delimiter='\t')
            for row in reader:
                genes.append(row)
        # Genes is a list of lists. For each list, first index: Ensemble ID, second: gene symbol.
        symbols = [L[1] for L in genes]
        seen = set()
        for symbol in symbols:
            if symbol not in seen:
                seen.add(symbol)
                # Get indices of duplicates in array symbols.
                dupes = [idx for idx, val in enumerate(symbols) if val == symbol]
                if len(dupes) > 1:
                    for idx in range(len(dupes)):
                        genes[dupes[idx]][1] = genes[dupes[idx]][1] + '.' + str(idx)

        return genes
