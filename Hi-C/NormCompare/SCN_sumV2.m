function [DataFN] = SCN_sumV2(DataF)

%add by Hongqiang Lyu
load('/media/biology/2b52abd3-7feb-4bc9-8c84-a8408d93f7ec/data/temp/tempmatrix-preSCN.mat');

DataFN = spdiags(spfun(@(x) 1./x,sum(DataF,2)),0,size(DataF,1),size(DataF,1))*DataF;
DataFN = spdiags(spfun(@(x) 1./x,sum(DataFN',2)),0,size(DataFN',1),size(DataFN',1))*DataFN';

DataFN = spdiags(spfun(@(x) 1./x,sum(DataFN',2)),0,size(DataFN',1),size(DataFN',1))*DataFN';
DataFN = spdiags(spfun(@(x) 1./x,sum(DataFN',2)),0,size(DataFN',1),size(DataFN',1))*DataFN';

DataFN = spdiags(spfun(@(x) 1./x,sum(DataFN',2)),0,size(DataFN',1),size(DataFN',1))*DataFN';
DataFN = spdiags(spfun(@(x) 1./x,sum(DataFN',2)),0,size(DataFN',1),size(DataFN',1))*DataFN';


%Remove nan and infs
DataFN(isnan(DataFN))=0;
DataFN(isinf(DataFN))=0;

%RE-SYMETRIZED

DataFN=(DataFN+DataFN')/2;

%add by Hongqiang Lyu
save('/media/biology/2b52abd3-7feb-4bc9-8c84-a8408d93f7ec/data/temp/tempmatrix-aftSCN.mat', 'DataFN');

end
