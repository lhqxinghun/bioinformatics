function RankedFea =  QPFS(inputfile)
%Description: Read matrices containing feature and class information of samples, output a specific number of ranked features according to QPFS criterion.
%inputfile    - Input/Sample files in UCI format.
%RankedFea    - Output/Ranked features.
%Example:
%RankedFea =  QPFS('diabetic.data');

%Configuration
maxbufsize=405900;
magnification=10;
%Read files
cd DATA;
fidinputfile=fopen(inputfile,'r');
inputline1=textscan(fidinputfile,'%s','BufSize',maxbufsize,'Delimiter','\n');
fclose(fidinputfile);
cd ..;
%Initialization
[row,col]=size(inputline1{1,1});
classindex={};
classnum=0;
vectornum=0;
rawfeavector={};
%Read feature vectors corresponding to different classes
for i=1:row
    inputline2=textscan(inputline1{1,1}{i,1},'%s','BufSize',maxbufsize,'Delimiter',',');
    [row1,col1]=size(inputline2{1,1});
    bpriorclass=ismember(classindex,inputline2{1,1}{row1,1});
    if(any(bpriorclass) == 0)
        classnum=classnum+1;
        classindex{classnum}=inputline2{1,1}{row1,1};
        vectornum(classnum)=1;
        fealen=row1-1;
        for j=1:fealen
            rawfeavector{classnum}{vectornum(classnum),j}=str2num(inputline2{1,1}{j,1});
        end    
    else        
        [bpcrow,bpcol]=size(bpriorclass);
        for j=1:bpcol
            if(bpriorclass(j)~=0)
                break;
            end
        end
        tempclassnum=j;
        vectornum(tempclassnum)=vectornum(tempclassnum)+1;
        fealen=row1-1;
        for j=1:fealen
            rawfeavector{tempclassnum}{vectornum(tempclassnum),j}=str2num(inputline2{1,1}{j,1});
        end
    end
end
%Data prepartion for QPFS
entrynum=0;
for i=1:classnum
    for j=1:vectornum(i)
        entrynum=entrynum+1;
        for m=1:fealen
            fearray(entrynum,m)=rawfeavector{i}{j,m};
        end
        classflag(entrynum)=i;
    end
end
fearray=QuantileDiscretize(fearray,magnification);
clear inputline1 inputline2 rawfeavector;
%QPFS feature selection
cd QPFS;

%Paramater for UCI data(medium-dimensional)
begintime=cputime;
QPFSFS=myQPFS(fearray,classflag');
runtime=cputime-begintime
expfeanum=fealen;

% %Paramater for GEMS and GEO data(high-dimensional, P>>N)
% begintime=cputime;
% QPFSFS=myQPFS(fearray,classflag');
% runtime=cputime-begintime
% expfeanum=200;

cd ..;
RankedFea=QPFSFS';
RankedFea=RankedFea(1:expfeanum);
clear;
end

function b=QuantileDiscretize(a,d)
if nargin<2 d=3;end;

[n dim]=size(a);
b=zeros(n,dim);
for i=1:dim
   b(:,i)=doDiscretize(a(:,i),d);
end
b=b+1;
end

% ----------------------------------------
function y_discretized= doDiscretize(y,d)
% ----------------------------------------
% discretize a vector
ys=sort(y);
y_discretized=y;

pos=ys(round(length(y)/d *[1:d]));
for j=1:length(y)
    y_discretized(j)=sum(y(j)>pos);
end
end