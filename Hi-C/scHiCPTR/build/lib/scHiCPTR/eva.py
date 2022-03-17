import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import igraph
from sklearn.metrics import roc_curve, auc
from scipy.stats import kendalltau
from scipy.sparse import csr_matrix
from collections import Counter
from typing import Optional, Union
from sklearn.metrics.cluster import normalized_mutual_info_score as nmi_score
from sklearn.metrics import adjusted_rand_score as ari_score
from sklearn import metrics
from munkres import Munkres, print_matrix

def computing_AUC(rankList, key:str):
    y_true = rankList['bench']
    y_score = np.max(rankList[key])-rankList[key]
    fpr, tpr, threshold = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)
    if roc_auc < 0.5:
        roc_auc = 1- roc_auc
    # plt.plot(fpr, tpr)
    return roc_auc

def evaluation_ranks_AUC(stages:list, predict_order:pd.DataFrame, dataset:str, key: str) -> pd.DataFrame:
    roc_auc_DF = list()
    if dataset == 'cellcycle':
        stage_ind = 'cell_cycle_stages'
    if dataset == 'gse':
        stage_ind = 'ESC_stages'

    for k, stage in enumerate(stages):
        stage1 = stage
        stage2 = stages[(k+1)%len(stages)]

        rank_list = predict_order.loc[(predict_order[stage_ind]==stage1)|(predict_order[stage_ind]==stage2)]
        rank_list['bench'] = 0
        rank_list.loc[rank_list[stage_ind]==stage2, 'bench'] = 1
        roc_auc = computing_AUC(rank_list, key)
        # print(stage1, stage2, roc_auc)
        roc_auc_DF.append(roc_auc)
    for k, stage in enumerate(stages):
        print(stage,stages[(k+1)%len(stages)], roc_auc_DF[k])
    roc_auc_DF=pd.DataFrame(roc_auc_DF,index=stages)    
    return roc_auc_DF

def compute_tau_b(x, y):
    coef, p = kendalltau(x, y)
    return coef, p

def kendall_rank_coefficient(predict_order: pd.DataFrame, true_order: Optional[list]=None, dataset: Optional[str]='cellcycle', key: str='pseudotime'):
    
    if true_order is None:
        if dataset == 'cellcycle':
            true_order = [0]*280+[1]*303+[2]*262+[3]*326
        if dataset == 'gse':
            true_order = [0]*55+[1]*94+[2]*87+[3]*101+[4]*171    
    true_order_sorted = sorted(true_order)
    # rank = dict(Counter(true_order_sorted))
    predict_order['index'] = np.array(range(len(predict_order)))
    predict_order = predict_order.sort_values(by=[key])
    predict_order['rank'] = np.array(true_order_sorted)
    predict_order = predict_order.sort_values(by='index')
    coef, p = compute_tau_b(true_order, predict_order['rank'].values)
    print(coef, p)
    return coef, p

def modularity_cal(
    adj: Union[csr_matrix, np.ndarray],
    groups: Optional[pd.DataFrame] = None,
) -> float:
    
    if isinstance(adj, csr_matrix):
        adj = adj.todense()
    
    g = igraph.Graph.Weighted_Adjacency(matrix=adj)
    Q = g.modularity(membership=groups['louvain'].cat.codes.values, weights='weight')
    print(Q)
    return Q

def cluster_acc(y_true, y_pred):
    y_true = y_true - np.min(y_true)

    l1 = list(set(y_true))
    numclass1 = len(l1)

    l2 = list(set(y_pred))
    numclass2 = len(l2)

    ind = 0
    if numclass1 != numclass2:
        for i in l1:
            if i in l2:
                pass
            else:
                y_pred[ind] = i
                ind += 1

    l2 = list(set(y_pred))
    numclass2 = len(l2)

    if numclass1 != numclass2:
        print('error')
        return

    cost = np.zeros((numclass1, numclass2), dtype=int)
    for i, c1 in enumerate(l1):
        mps = [i1 for i1, e1 in enumerate(y_true) if e1 == c1]
        for j, c2 in enumerate(l2):
            mps_d = [i1 for i1 in mps if y_pred[i1] == c2]
            cost[i][j] = len(mps_d)

    # match two clustering results by Munkres algorithm
    m = Munkres()
    cost = cost.__neg__().tolist()
    indexes = m.compute(cost)

    # get the match results
    new_predict = np.zeros(len(y_pred))
    for i, c in enumerate(l1):
        # correponding label in l2:
        c2 = l2[indexes[i][1]]

        # ai is the index with label==c2 in the pred_label list
        ai = [ind for ind, elm in enumerate(y_pred) if elm == c2]
        new_predict[ai] = c

    acc = metrics.accuracy_score(y_true, new_predict)
    f1_macro = metrics.f1_score(y_true, new_predict, average='macro')
    precision_macro = metrics.precision_score(y_true, new_predict, average='macro')
    recall_macro = metrics.recall_score(y_true, new_predict, average='macro')
    f1_micro = metrics.f1_score(y_true, new_predict, average='micro')
    precision_micro = metrics.precision_score(y_true, new_predict, average='micro')
    recall_micro = metrics.recall_score(y_true, new_predict, average='micro')
    return acc, f1_macro

def cluster_eva(y_true, y_pred, epoch=0):
    # accuracy, F1
    acc, f1 = cluster_acc(y_true, y_pred)
    # NMI
    nmi = nmi_score(y_true, y_pred, average_method='arithmetic')
    # ARI
    ari = ari_score(y_true, y_pred)
    print(epoch, ':acc {:.4f}'.format(acc), ', nmi {:.4f}'.format(nmi), ', ari {:.4f}'.format(ari),
            ', f1 {:.4f}'.format(f1))



    
