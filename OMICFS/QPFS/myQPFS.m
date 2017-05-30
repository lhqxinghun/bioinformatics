%Nguyen X. Vinh, Jeffrey Chan, Simone Romano and James Bailey, "Effective Global Approaches for Mutual Information based Feature Selection". 
%To appear in Proceeedings of the 20th ACM SIGKDD Conference on Knowledge Discovery and Data Mining (KDD'14), August 24-27, New York City, 2014. 
% (C) 2014 Nguyen Xuan Vinh   
% Email: vinh.nguyen@unimelb.edu.au, vinh.nguyenx@gmail.com
% Performing the following mutual informattion based feature selection
% approaches:
% - Maximum relevance (maxRel)
% - Minimum redundancy maximum relevance (MRMR)
% - Minimum redundancy (minRed)
% - Quadratic programming feature selection (QPFS)
% - Mutual information quotient (MIQ)
function QPFSFS=myQPFS(a,C)

[n dim]=size(a);
H=zeros(dim,dim);
f=zeros(dim,1);
maxE=log(max(max(a)));

fprintf('QPFS started, calculating the MI matrix...')

H=computeMImatrix_4([a C]);
f=H(1:end-1,end);
H=H(1:end-1,1:end-1);

%-doing QPFS
%H=(H+lambda*eye(dim));
tic;
evalue=eig(H);
if evalue(1)<-10^-3
    fprintf('QPFS: non-positive Q\n');
    for i=1:dim H(i,i)=H(i,i)+3;end;
end

mq=sum(sum(H))/dim^2;
mf=sum(f)/dim;
alpha=mq/(mq+mf);
f=-alpha*f;
H=(1-alpha)*H;

A=-eye(dim);  % x_i >=0 contraints
b=zeros(dim,1);
Aeq=ones(1,dim);
beq=1;

options = optimset('Algorithm','interior-point-convex');
fprintf('Solving the QPFS formulation...\n');
x = quadprog(H,f,A,b,Aeq,beq,[],[],[],options);

f=0.5*(1-alpha)*x'*H*x -alpha*f'*x;
y=zeros(dim,2);
y(:,1)=-x;
for i=1:dim y(i,2)=i;end;
y=sortrows(y);
y=[y(:,2) -y(:,1)];
QPFSFS=y(:,1);
end
