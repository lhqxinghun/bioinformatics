function [trPCA preserve preserveN evects mdata sdata]=trkpca(kernel,T,NP)
N=size(kernel,1);

unit = ones(N, N)/N;
% centering in feature space!
K_n = kernel - unit*kernel - kernel*unit + unit*kernel*unit;

[cvv,cdd] = eig(K_n);
[zz,ii] = sort(diag(cdd));
ii = flipud(ii);
evects = cvv(:,ii);
cdd = diag(cdd);
evals = cdd(ii);

if NP==-1
    NP=min(find(cumsum(evals)>=(T*sum(evals))));
end
preserve=sum(evals(1:NP))/sum(evals);
preserveN=NP;

trPCA=K_n*evects(:,1:NP);
[trPCA mdata sdata]=autosc(trPCA);
