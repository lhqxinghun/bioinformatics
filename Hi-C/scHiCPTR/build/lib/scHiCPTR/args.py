
mm9dim = [197195432,181748087,159599783,155630120,152537259,149517037,152524553,131738871,124076172,129993255,121843856,121257530,120284312,125194864,103494974,98319150,95272651,90772031,61342430,166650296]
chrom =  [str(i+1) for i in range(19)] + ['X']
chromsize = {chrom[i]:mm9dim[i] for i in range(len(chrom))}

## default args
## feature args
impute_args = {"pad": 1, "rp": 0.5, "prct": 20, "ncpus": 10}
pca_args = {"n_components": 120}
umap_args = {"n_neighbors":100, "min_dist":0.5, "n_components":20, "metric":'euclidean'}

## KNN args
knn_args = {"n_neighbors": 100, "metric":'euclidean'}

## pruning args
prune_args = {"lambda_": 0.5, "threshold": 0.7, "louvain_res": 0.5}