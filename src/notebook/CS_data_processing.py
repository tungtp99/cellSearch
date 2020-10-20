import numpy as np

from anndata import AnnData
import scanpy as sc
import math

def get_highly_genes(profile):
    return sc.pp.highly_variable_genes(AnnData(profile), inplace = False)

def do_standardization(X):
    
    # Setup standardization config
    adata = AnnData(X)
    
    ###return adata.X
    
    # Perform standardization
    sc.pp.normalize_total(adata, target_sum=1e6, exclude_highly_expressed = True, inplace=True)

    # Apply logarithm
    sc.pp.log1p(adata, copy = False)
    
    # Return result
    return adata.X

def do_pca(X, n_components):
    from sklearn.decomposition import IncrementalPCA
    from sklearn.utils import resample

    transformer = IncrementalPCA(n_components = n_components, batch_size = None)
    
    cell_counts = int(X.shape[0])
    sample_fraction = 1.0 - (math.e ** (-10000 / cell_counts))
    
    if (sample_fraction < 1.0):
        sample_size = int(cell_counts * sample_fraction)
        sample_data = resample(X, n_samples = sample_size, replace = False, random_state = 42)
    transformer.partial_fit(sample_data)
    
    return transformer.transform(X)


