import numpy as np
import platform

from scipy.sparse import data

from numpy.lib.npyio import savetxt
from scHiCPTR import scHiCPTR, eva, vis

if platform.system().lower() == 'windows':
    import umap.umap_ as umap
if platform.system().lower() == 'linux':
    import umap

## The input data is raw matrices.
# if __name__ == '__main__':

#     network = np.loadtxt('../data/SampleList.txt' , dtype=np.str)
#     mm9dim = [197195432,181748087,159599783,155630120,152537259,149517037,152524553,131738871,124076172,129993255,121843856,121257530,120284312,125194864,103494974,98319150,95272651,90772031,61342430,166650296]
#     chrom =  [str(i+1) for i in range(19)] + ['X']
#     chromsize = {chrom[i]:mm9dim[i] for i in range(len(chrom))}

#     ## KNN args
#     knn_args = {"n_neighbors": 100, "metric":'euclidean'}

#     ## Pruning args
#     prune_args = {"lambda_": 0.5, "threshold": 0.7, "louvain_res": 0.5}

#     ## Calculate pseudotime
#     trajectory = scHiCPTR.scHiCPTR(data=network, embedding=None, chromsize=chromsize, resolution=1000000, knn_args=knn_args, prune_args=prune_args, start=270)
#     result = trajectory._trajectory 


## The input data is the low-dimension representations.
if __name__ == '__main__':

    embeds = np.loadtxt('../data/embeds/cellcycle_250k.txt' , dtype=np.str)
    mm9dim = [197195432,181748087,159599783,155630120,152537259,149517037,152524553,131738871,124076172,129993255,121843856,121257530,120284312,125194864,103494974,98319150,95272651,90772031,61342430,166650296]
    chrom =  [str(i+1) for i in range(19)] + ['X']
    chromsize = {chrom[i]:mm9dim[i] for i in range(len(chrom))}

    ## KNN args
    knn_args = {"n_neighbors": 100, "metric":'euclidean'}

    ## Pruning args
    prune_args = {"lambda_": 0.5, "threshold": 0.7, "louvain_res": 0.5}

    ## Calculate pseudotime
    trajectory = scHiCPTR.scHiCPTR(data=None, embedding=embeds, chromsize=chromsize, knn_args=knn_args, prune_args=prune_args, start=270)
    result = trajectory._trajectory 

    ## Evaluation results
    result['cell_cycle_stages'] = ['G1']*280+['ES']*303+['MS']*262+['G2']*326
    result.to_csv('../data/results/pseudotime.csv')
    stages = ['G1', 'ES', 'MS', 'G2']
    eva.evaluation_ranks_AUC(stages, result, 'cellcycle', key='pseudotime')
    eva.kendall_rank_coefficient(predict_order=result, dataset='cellcycle', key='pseudotime')

    ## Visualize pseudotime
    vis.vis_trajectory(trajectory, save_path='../data/results')