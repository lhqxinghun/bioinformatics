function RankedFea =  KCCAmRMR(inputfile)
%Description: Read matrices containing feature and class information of samples, output a specific number of ranked features according to KCCAmRMR criterion.
%inputfile    - Input/Sample files in UCI format.
%RankedFea    - Output/Ranked features.
%Example:
%RankedFea =  KCCAmRMR('diabetic.data');

%Configuration
maxbufsize=405900;
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
%Data prepartion for KCCAmRMR
entrynum=0;
for i=1:classnum
    for j=1:vectornum(i)
        entrynum=entrynum+1;
        for m=1:fealen
            fearray(entrynum,m)=rawfeavector{i}{j,m};
        end
        classflag(entrynum,1)=i;
    end
end
clear inputline1 inputline2 rawfeavector;
%KCCAmRMR feature selection
cd KCCAmRMR;

%Paramater for UCI data(medium-dimensional)
begintime=cputime;
mbins=100;
RankedFea=kccamrmr(fearray,classflag',mbins,fealen)
runtime=cputime-begintime

% %Paramater for GEMS and GEO data(high-dimensional, P>>N)
% begintime=cputime;
% mbins=100;
% expfeanum=200;
% RankedFea=kccamrmr(fearray,classflag',mbins,expfeanum)
% runtime=cputime-begintime

cd ..;
clear;
end