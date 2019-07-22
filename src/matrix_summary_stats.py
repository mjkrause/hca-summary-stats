#!/usr/bin/env python3

import os
import csv
import shutil
import gzip
import pandas as pd
from zipfile import ZipFile
from more_itertools import first


class MatrixSummaryStats:

    def __init__(self, mtxdir):
        self.mtxdir = mtxdir
        self.genes = None
        self.barcodes = None

    def unzip_files(self, path=None) -> None:
        with ZipFile(self.mtxdir) as archive:
            archive.extractall(path=path)
        matrix_path = first(os.path.splitext(self.mtxdir))  # remove extension ".zip"
        original_dir = os.getcwd()
        os.chdir(matrix_path)
        files = os.listdir('.')
        assert ['cells.tsv.gz', 'genes.tsv.gz', 'matrix.mtx.gz'] == files

        for file in files:
            with gzip.open(file, 'rb') as f_in:
                if f_in.name == 'genes.tsv.gz' or f_in.name == 'cells.tsv.gz':
                    self.preprocessing(f_in)  # writes files to disk
                elif f_in.name == 'matrix.mtx.gz':
                    outfile = first(os.path.splitext(f_in.name))
                    with open(outfile, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
            os.remove(file)
        files = os.listdir('.')
        assert 'barcodes.tsv' in files
        assert 'genes.tsv' in files
        self.genes = os.path.abspath('genes.tsv')
        self.barcodes = os.path.abspath('barcodes.tsv')
        os.chdir(original_dir)

    @staticmethod
    def preprocessing(fileobj: gzip.GzipFile) -> None:
        """Preprocessing TSV files from Matrix Service in order to use scanpy methods on matrix."""
        f = pd.read_table(fileobj, sep='\t')  # Pandas dataframe
        col_to_keep = []
        if fileobj.name == 'genes.tsv.gz':
            col_to_keep = ['featurekey', 'featurename']
            assert col_to_keep[0] in f.columns
            assert col_to_keep[1] in f.columns
        elif fileobj.name == 'cells.tsv.gz':
            fileobj.name = 'barcodes.tsv.gz'
            col_to_keep = 'cellkey'
            assert col_to_keep in f.columns
        f_new = f[col_to_keep]
        # Write to file without column or row headers.
        f_new.to_csv(first(os.path.splitext(fileobj.name)), index=False, header=False, sep='\t')

    def eliminate_dupes(self) -> list:
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
