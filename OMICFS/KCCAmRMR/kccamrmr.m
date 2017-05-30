%INPUT
%mydata: dataset containing CONTINUOUS FEATURES (all features are assumed to be continuous and normalized)
%classlabels: vector that contains the corresponding class labels
%mbins: number of bins that will be used for the discretization of the continuous features

%OUTPUT
%kccamrmr_ranking: Ranking of the features by KCCAmRMR algorithm

%For more details, refer to:
% Sakar, C.O., Kursun, O., Gurgen F., “A Feature Selection Method Based on Kernel Canonical
% Correlation Analysis and the Minimum Redundancy Maximum Relevance Filter Method? <br>
% Expert Sytems with Applications, vol. 39(3), pp. 3333-3344, 2012.
                
function [kccamrmr_ranking] = kccamrmr(mydata, classlabels, mbins,K)
    
[CC kk] = size(unique(classlabels));

[N n]=size(mydata);   %number of samples and number of columns
out=classlabels;
if CC == 2
    out = out + 1;
end

for i=1:n
    nanvector = isnan(mydata(:,i));
    sumCol = 0;
    for j=1:N
       if nanvector(j) == 0
           sumCol = sumCol + mydata(j,i);
       end
    end
    
    if sum(nanvector) == N
        mydata(:,i) = 0;
    else
        meanCol = sumCol / sum(nanvector == 0);
        for j=1:N
            if nanvector(j) == 1
                mydata(j,i) = meanCol;
            end
        end
    end
end

T=out;
out=zeros(N,CC);
for samp=1:N
    out(samp,T(samp))=1;
end

dme=mean(mydata);
dst=std(mydata);

n_mydata=zeros(N,n);
for i=1:n
    if dst(i)~=0
        n_mydata(:,i)=autosc(mydata(:,i));
        mydata(:,i)=floor(((mydata(:,i)-dme(i))/dst(i))+0.5);
    end

    f=find(mydata(:,i)<-mbins);
    if ~isempty(f)
        mydata(f,i)=-mbins;
    end

    f=find(mydata(:,i)>mbins);
    if ~isempty(f)
        mydata(f,i)=mbins;
    end
end

go=1;

begintime=cputime;

okernel=exp(-go*(dist(out,out').^2));

[otrPCA opreserve opreserveN oevects omdata osdata]=trkpca(okernel,.9,-1); 
var=0;
for svar=1:n    %for each feature in the dataset
%     disp(svar)
    d=n_mydata(:,svar);
    if length(unique(d))>2
        gd=10;              
        dkernel=exp(-gd*(dist(d,d').^2));
        [trPCA preserve preserveN evects mdata sdata]=trkpca(dkernel,.9,-1);   

        [wx wy cors]=cca(trPCA',otrPCA');
        f=find(cors>0.1);
        lf=length(f);
        R(:,var+1:var+lf)=trPCA*wx(:,f);
        subs{svar}=[var+1:var+lf];
        trcors{svar}=cors(f);
        var=var+lf;
    else
        lf=1;
        R(:,var+1:var+lf)=d;
        subs{svar}=[var+1:var+lf];
        if length(unique(d))==2
            ccc=corrcoef([d out]);
            trcors{svar}=max(ccc(2:end));
        else
            trcors{svar}=0;
        end
        var=var+lf;
    end
end

runtime=cputime-begintime;

me=mean(R);
st=std(R);
dim=size(R,2);

for i=1:dim
    n_mydata(:,i)=floor(((R(:,i)-me(i))/st(i))+0.5);

    f=find(n_mydata(:,i)<-mbins);
    if ~isempty(f)
        n_mydata(f,i)=-mbins;
    end

    f=find(n_mydata(:,i)>mbins);
    if ~isempty(f)
        n_mydata(f,i)=mbins;
    end
end

begintime=cputime;

res=zeros(n);
for i=1:n
    misorig(i)=mutualinfo(mydata(:,i),T);
end
runtime=cputime-begintime;
mypow=1;
mypow2=mypow*2;

ii=0;
for i=1:n
    mis(i)=0;
    thecor=trcors{i};
    for k=1:length(subs{i})
        mis(i)=mis(i)+mutualinfo(n_mydata(:,ii+k),T)*(thecor(k)^mypow2);
    end
    jj=0;
    for j=1:i-1
        thecorj=trcors{j};
        for k=1:length(subs{i})
            for l=1:length(subs{j})
                res(i,j)=res(i,j)+mutualinfo(n_mydata(:,ii+k),n_mydata(:,jj+l))*(thecor(k)^mypow)*(thecorj(l)^mypow);
            end
        end
        jj=jj+length(subs{j});
    end
    ii=ii+length(subs{i});
end
res=max(res,res');
resorig=res;
res=res-(eye(n)*1000);


clear kccamrmr_ranking t_mi c_mi mi_array
[tmp, idxstart] = max(misorig);
kccamrmr_ranking(1) = idxstart;
mis(idxstart)=-1000;
[tmp, idxleft] = sort(-mis);
idxleft=idxleft(1:K-1);
for k=2:K,
    clear t_mi clear mi_array;
    ncand = length(idxleft);
    curlastfea = length(kccamrmr_ranking);
    for i=1:ncand,
        t_mi(i) = mis(idxleft(i));
        mi_array(i) = mean(res(idxleft(i),kccamrmr_ranking(1:curlastfea)));
    end
    [tmpval, tmpidx] = max(t_mi - mi_array);
    kccamrmr_ranking(k) = idxleft(tmpidx); idxleft(tmpidx) = [];
end

