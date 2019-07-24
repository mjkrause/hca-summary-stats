#!/usr/bin/env python3

import os
import time
import csv
import shutil
import gzip
import pandas as pd
import scanpy as sc
import config
import requests
import logging
from zipfile import ZipFile
from more_itertools import first

logger = logging.getLogger(__name__)


class MatrixSummaryStats:

    def __init__(self):
        self.matrix_zipfile_name = None
        self.matrix_path = None
        self.matrix_response = None
        self.genes = None
        self.barcodes = None

    def get_expression_matrix(self) -> None:
        status_response = self._request_matrix()
        s3_download_url = status_response.json()['matrix_url']
        self.matrix_response = requests.get(s3_download_url, stream=True)
        self.matrix_zipfile_name = os.path.basename(s3_download_url)

    def unzip_files(self, path=None) -> None:
        root_dir = os.getcwd()
        if path:
            os.chdir(path)
        with open(self.matrix_zipfile_name, 'wb') as matrix_zip_file:
            shutil.copyfileobj(self.matrix_response.raw, matrix_zip_file)
        ZipFile(self.matrix_zipfile_name).extractall()
        self.matrix_path = first(os.path.splitext(self.matrix_zipfile_name))  # remove extension ".zip"
        os.chdir(self.matrix_path)
        files = os.listdir('.')

        assert len(files) == 3
        assert 'cells.tsv.gz' in files
        assert 'genes.tsv.gz' in files
        assert 'matrix.mtx.gz' in files

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
        os.chdir(root_dir)

    def push_to_s3(self):
        pass

    def create_images(self) -> None:
        # Highest-expressing genes.
        adata = sc.read_10x_mtx(self.matrix_path, var_names='gene_symbols', cache=True)
        adata.var_names_make_unique()
        sc.pl.highest_expr_genes(adata, n_top=20, save='.png', show=False)

        # Highest-variable genes:
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e3)
        sc.pp.log1p(adata)  # logarithmize
        adata.raw = adata  # save raw data for later use
        sc.pp.log1p(adata)
        adata.raw = adata
        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
        sc.pl.highly_variable_genes(adata, save='.png')

    @staticmethod
    def _request_matrix() -> requests.models.Response:
        # Parameters.
        feature = 'gene'
        format_ = 'mtx'
        project = 'Single cell transcriptome analysis of human pancreas'  # get single project
        project_field = 'project.project_core.project_short_name'
        min_cell_count = 300
        min_cell_count_field = 'genes_detected'

        hca_matrix_service_url = config.endpoints['hca_matrix_service_url']

        list_projects_url = hca_matrix_service_url + '/filters/' + project_field
        assert project in list(requests.get(list_projects_url).json()['cell_counts'])

        payload = {'feature': feature,
                   'format': format_,
                   'filter': {
                       'op': 'and',
                       'value': [
                           {'op': '=',
                            'value': project,
                            'field': project_field},
                           {'op': '>=',
                            'value': min_cell_count,
                            'field': min_cell_count_field}
                           ]
                       }
                   }
        logger.info(f'Requesting expression matrix for project {project}')
        response = requests.post(hca_matrix_service_url + '/matrix', json=payload)

        while True:
            status_response = requests.get(hca_matrix_service_url + '/matrix/' +
                                           response.json()['request_id'])
            if status_response.json()['status'] == 'Complete':
                break
            logger.info(f'{status_response.json()["status"]} ...')
            time.sleep(30)

        return status_response

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
