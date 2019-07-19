import loompy
import numpy as np

from azul import config
from more_itertools import first

class ProjSummaryStats:

    def __init__(self, loomfile: str, var_top_hits: int):
        self.loomfile = loomfile
        self.tophits = var_top_hits
        self.high_expressing_genes = []
        #self.adata = sc.read(loomfile)

    def find_gene_with_max_var(self) -> list: # in each gene (=row)
        with loompy.connect(self.loomfile) as ds:
            # Create list of tuples, where the first tuple element denotes the gene name,
            # and the second tuple element denotes the variance of that gene across all cells.
            high_variance_genes = [(gene, ds[ds.ra.Gene == gene, :].var()) for gene in ds.ra.Gene]
            high_variance_genes.sort(key=lambda tup: tup[1], reverse=True)  # sort in place
            return high_variance_genes[:self.tophits]

    # def find_gene_with_max_expression(self) -> list: # in each cell
    #     high_exp_genes = []
    #     with loompy.connect(self.loomfile) as ds:
    #         arr = ds[:, ds.ca.CellID == ds.ca.CellID[cell_id]]
    #         idx = np.where(arr == arr.max())  # find maximum-expressing gene in cell
    #         return high_exp_genes.append((cell_id, set(ds.ra.Gene[first(idx)])))
    #         #return [cell_id for cell_id in range(ds.ca['CellID'].size)])

    def _load_loomfile(self) -> tuple:
        with loompy.connect(self.loomfile) as ds:
            col_headers = ds.ca.CellID
            row_headers = ds.ra.Gene
            return (col_headers, row_headers, ds[:,:])

    def compute(self):
        col_headers, row_headers, gene_by_cell_matrix = self._load_loomfile()

        expression_variance = self._find_gene_with_max_var(gene_by_cell_matrix, row_headers)
        gene_with_max_expression_per_cell = self._gene_with_max_expression_per_cell(gene_by_cell_matrix, row_headers)

        return (expression_variance, gene_with_max_expression_per_cell)

    def _find_gene_with_max_var(self, M: np.ndarray, row_headers: np.ndarray) -> list:
        expression_variance = np.var(M, axis=1)  # across all genes
        index = np.flip(np.argsort(expression_variance, axis=0))  # ordered desc.
        return list(row_headers[index][:self.tophits])

    def _gene_with_max_expression_per_cell(self, M: np.ndarray, row_headers: np.ndarray) -> list:
        index = np.argsort(M, axis=0)[-1, :]
        return list(row_headers[index])
