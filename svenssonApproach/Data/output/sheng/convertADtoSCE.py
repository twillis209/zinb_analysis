import anndata
import anndata2ri
import argparse
import os
from rpy2 import robjects

parser = argparse.ArgumentParser(description="Reads AnnotatedDataFrame object from file, converts to R object of the SingleCellExperiment class, and writes the R object to an RData file.")
parser.add_argument("-if", "--h5adFilePath", type=str, help="Path to h5ad file containing AnnotatedDataFrame")
parser.add_argument("-of", "--rdataFilePath", type=str, help="Path to RData file to contain converted SingleCellExperiment objects")

args = parser.parse_args()

anndata2ri.activate()

ad = anndata.read(args.h5adFilePath)

name = os.path.splitext(os.path.basename(args.h5adFilePath))[0]

sexp=anndata2ri.py2rpy(ad)

robjects.r.assign(name, sexp)

robjects.r("save({},file=\'{}\')".format(name,args.rdataFilePath))

anndata2ri.deactivate()
