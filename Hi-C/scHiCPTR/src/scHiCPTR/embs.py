import os
import time
import logging
import numpy as np
from . import args
import umap
import anndata as ad
# import torch
# import torch.nn.functional as F
from multiprocessing import Pool
from scipy.sparse import csr_matrix
from typing import Optional
from pathlib import Path
from sklearn.decomposition import PCA

def neighbor_ave_cpu(A, pad):
	if pad==0:
		return A
	ngene, _ = A.shape
	ll = pad * 2 + 1
	B, C, D, E = [np.zeros((ngene + ll, ngene + ll)) for i in range(4)]
	B[(pad + 1):(pad + ngene + 1), (pad + 1):(pad + ngene + 1)] = A[:]
	F = B.cumsum(axis = 0).cumsum(axis = 1)
	C[ll :, ll:] = F[:-ll, :-ll]
	D[ll:, :] = F[:-ll, :]
	E[:, ll:] = F[:, :-ll]
	return (np.around(F + C - D - E, decimals=8)[ll:, ll:] / float(ll * ll))

def random_walk_cpu(A, rp):
	ngene, _ = A.shape
	A = A - np.diag(np.diag(A))
	A = A + np.diag(np.sum(A, axis=0) == 0)
	P = np.divide(A, np.sum(A, axis = 0))
	Q = np.eye(ngene)
	I = np.eye(ngene)
	for i in range(30):
		Q_new = (1 - rp) * I + rp * np.dot(Q, P)
		delta = np.linalg.norm(Q - Q_new)
		Q = Q_new.copy()
		if delta < 1e-6:
			break
	return Q

def impute_cpu(args):
	cell, c, ngene, pad, rp = args
	D = np.loadtxt(cell + '_chr' + c + '.txt')
	A = csr_matrix((D[:, 2], (D[:, 0], D[:, 1])), shape = (ngene, ngene)).toarray()
	A = np.log2(A + A.T + 1)
	A = neighbor_ave_cpu(A, pad)
	if rp==-1:
		Q = A[:]
	else:
		Q = random_walk_cpu(A, rp)
	return [cell, Q.reshape(ngene*ngene)]

def impute_embedding(
    network: np.ndarray,
    chromsize: dict,
    res: int = 1000000,
    impute_args: dict = None,
    pca_args: dict = None,
    umap_args: dict = None,
):
    """\

    """
    if network is None:
        logging.error("Missing key parameter '\network'\!")

    if chromsize is None:
        logging.error("Missing key parameter '\chromsize'\!")

    if impute_args is None:
        logging.warning("Using default imputaion parameters")
        impute_args = args.impute_args
    pad, rp, prct, ncpus = impute_args["pad"], impute_args["rp"], impute_args["prct"], impute_args["ncpus"]

    if pca_args is None:
        logging.warning("Using default pca parameters")
        pca_args = args.pca_args

    if umap_args is None:
        logging.warning("Using default umap parameters")
        umap_args = args.umap_args
    
    # The code references  "https://github.com/zhoujt1994/scHiCluster"
    matrix = []
    for i, c in enumerate(chromsize):
        ngene = int(chromsize[c] / res) + 1
        start_time = time.time()
        paras = [[cell, c, ngene, pad, rp] for cell in network]
        p = Pool(ncpus)
        result = p.map(impute_cpu, paras)
        p.close()
        index = {x[0]: j for j, x in enumerate(result)}
        Q_concat = np.array([result[index[x]][1] for x in network])
        if prct > -1:
            thres = np.percentile(Q_concat, 100 - prct, axis=1)
            Q_concat = (Q_concat > thres[:, None])
        end_time = time.time()
        print('Load and impute chromosome', c, 'take', end_time - start_time, 'seconds')
        ndim = int(min(Q_concat.shape) * 0.2) - 1
        pca = PCA(n_components=ndim)
        R_reduce = pca.fit_transform(Q_concat)
        matrix.append(R_reduce)
        print(c)
    matrix = np.concatenate(matrix, axis=1)
    pca = PCA(n_components=min(matrix.shape) - 1)
    data = pca.fit_transform(matrix)
    data = data[:, 0:pca_args["n_components"]]
    fit = umap.UMAP(**umap_args)
    mapper = fit.fit_transform(data) 

    return mapper
