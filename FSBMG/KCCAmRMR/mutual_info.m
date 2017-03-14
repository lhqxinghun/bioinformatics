%INPUT
%u: discrete variable 1
%v: discrete variable 2

%OUTPUT
%mi: mutual information between u and v

% Mutual information computation
%compute the joint distribution

function [mi]=mutual_info(u,v)

if (size(u,1)~=size(v,1))
    disp('Both vectors must be column vectors of the same size');
    return;
end
N=size(u,1);
us=unique(u');
vs=unique(v');
maxU=numel(us);
maxV=numel(vs);

a=zeros(maxU,1);
b=zeros(maxV,1);
i=0;
for k=us
    i=i+1;
    a(i)=sum(u==k);
end
j=0;
for l=vs
    j=j+1;
    b(j)=sum(v==l);
end
mi=0;
for i=1:maxU
    k=us(i);
    for j=1:maxV
        l=vs(j);

        cm=sum(u'==k & v'==l);
        if (cm>0)
            rt= cm/N;
            mi=mi+rt*log2(rt/(a(i)*b(j)/N^2));
        end
    end
end

end