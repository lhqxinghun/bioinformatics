function RankedFea = OMICFS(fearray,classflag,psfeanum,expfeanum)
%Description: Read matrices containing feature and class information of samples, output a specific number of ranked features according to OMICFS criterion.
%fearray      - Input/m*n feature array, in which m is the number of samples and n is the total number of candidate features.
%classflag    - Input/m*1 vector, each element in the vector is the flag number of the class it belongs to.
%psfeanum     - Input/The number of pre-screening features (psfeanum = expfeanum = the number of all features: means all the features are expected to be ranked).
%expfeanum    - Input/The number of ranked features expected.
%RankedFea    - Output/Ranked features.
%Example:
%RankedFea = OMICFS(fearray,classflag,psfeanum,expfeanum);

%Check the input parameters 
if (psfeanum<expfeanum)
    disp('Warning: the number of pre-screening features is smaller than expected features.');
    return; 
end
%Computer Max-Relevance
[MaxRelFeaNum,MaxRelMICValue] = MIC(fearray,classflag',psfeanum);
%Cut off completely noneffective features 
cutindex=find(MaxRelMICValue==0);
if (isempty(cutindex)==0)
    cutfea=MaxRelFeaNum(cutindex:length(MaxRelFeaNum));
    warningline=sprintf('Warning: a total of %d completely noneffective features exist.',length(cutfea));
    disp(warningline);
    MaxRelFeaNum=MaxRelFeaNum(1:cutindex-1);
    MaxRelMICValue=MaxRelMICValue(1:cutindex-1);
end
%Initialization
len=length(MaxRelFeaNum);
candfeaflag=ones(1,len-1);
actexpfeanum=min(expfeanum,len);
RankedFea(1)=MaxRelFeaNum(1);
temporthvector=fearray(:,MaxRelFeaNum(1));
orthvector(:,1)=temporthvector./norm(temporthvector);
for i=1:actexpfeanum-1
    index1=0;
    temporthvector=[];
    tempFSMICscore=[];
    tempFSMICnum=[];
    for j=1:len-1
        if(candfeaflag(j)==1)
            index1=index1+1;
            temporthvector(:,index1)=fearray(:,MaxRelFeaNum(j+1));
            for m=1:i
                temporthvector(:,index1)=temporthvector(:,index1)-fearray(:,MaxRelFeaNum(j+1))'*orthvector(:,m)./(orthvector(:,m)'*orthvector(:,m))*orthvector(:,m);
            end
            temporthvector(:,index1)=temporthvector(:,index1)./norm(temporthvector(:,index1));
        end     
    end
    %Computer Max-Relevance and Min-Redundancy   
    [MaxRelMinRedFeaNum,MaxRelMinRedMICValue] = MIC(temporthvector,classflag',index1);
    index2=0;
    for j=1:len-1
        if(candfeaflag(j)==1)
          index2=index2+1;
          tempFSMICscore(index2)=MaxRelMinRedMICValue(find(MaxRelMinRedFeaNum==index2));
          tempFSMICnum(index2)=MaxRelFeaNum(j+1);
        end
    end
    [maxvualue,maxindex]=max(tempFSMICscore);
    candfeaflag(find(MaxRelFeaNum==tempFSMICnum(maxindex))-1)=0;
    orthvector(:,i+1)=temporthvector(:,maxindex);
    RankedFea(i+1)=tempFSMICnum(maxindex);
end
% Splice necessary noneffective features at the tail
if (actexpfeanum<expfeanum && isempty(cutindex)==0)
    splfeanum=min(length(cutfea),expfeanum-actexpfeanum);
    RankedFea=[RankedFea,cutfea(1:splfeanum)];
end
end


%Description: Read matrices containing feature and class information of samples, output the ranked features and corresponding MIC value.
%fearray      - Input/m*n feature array, in which m is the number of samples and n is the total number of candidate features.
%classflag    - Input/m*1 vector, each element in the vector is the flag number of the class it belongs to.
%varargin     - Input/The number of ranked features expected (The parameter is optional, all the ranked features will be output without this parameter).
%FeaNum       - Output/Ranked features. Each feature is indicated with its serial number.
%MICValue     - Output/MIC value corresponding to each ranked feature.
%Example1:
%RankedFea = MIC(fearray,classflag);
%Example2:
%RankedFea = MIC(fearray,classflag,psfeanum);

function [FeaNum,MICValue] = MIC(fearray,classflag,varargin)
%Check the input parameters
[inrow,incol]=size(fearray);
if (nargin>2)
    expfeanum=varargin{1};
elseif (nargin==2)
    expfeanum=incol;
end 
if (incol<expfeanum)
    disp('Warning: the number of required features is bigger than input candidate features.');
    return; 
end
cd MINE/MIC;
%Calculate MIC score
for i=1:incol
    fearray=round(fearray*10^4)/10^4;
    minestats= mine(fearray(:,i)',classflag);
    MICValue(i)=minestats.mic;    
end
%Sort
[MICValue,FeaNum]=sort(MICValue,'descend');
[outrow,outcol]=size(MICValue);
%Check the result
if (incol~=outcol)
    for i=1:incol
        if(isempty(find(FeaNum==i)))
            FeaNum=[FeaNum,i];
            MICValue=[MICValue,0];
        end
    end    
end
if (incol>expfeanum)
    FeaNum=FeaNum(1:expfeanum);
    MICValue=MICValue(1:expfeanum);
end
cd ../..;
end