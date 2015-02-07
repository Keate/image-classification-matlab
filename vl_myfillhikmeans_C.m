function mytree = vl_myfillhikmeans_C(Ensemble,EnsembleWeight,K, L)

if isempty(Ensemble)
    load Ensemble
end

MaxEnsembleSize=1e+8;
randindex=randperm(size(Ensemble,1));
randindex=randindex(1:min(size(Ensemble,1),MaxEnsembleSize));
Ensemble=Ensemble(randindex,:);
EnsembleWeight=EnsembleWeight(randindex);

if size(Ensemble,1)<10*K,
    %Ensemble=[Ensemble;(repmat(mean(Ensemble),K+1-size(Ensemble,1),1)+repmat(std(single(Ensemble),0),K+1-size(Ensemble,1),1) .* randn(K+1-size(Ensemble,1),size(Ensemble,2)))]; %append
    Ensemble=[Ensemble;repmat(Ensemble,min(K,ceil((5*K+1)/size(Ensemble,1))),1)]; %append
    EnsembleWeight=[EnsembleWeight;repmat(EnsembleWeight,min(K,ceil((5*K+1)/size(EnsembleWeight,1))),1)]; 
end


%[C0,Index]=mykmeans_c(single(Ensemble),K); Index=Index+1; C0=single(C0'); % my kmeans with souce code

%[C0,Index]=vl_kmeans(single((bsxfun(@times,Ensemble(:,1:end-1),[ones(1,128),gps_sift_ratio*ones(1,2)]))'), K, 'distance', 'l2', 'algorithm', 'elkan') ;  % fast kmeans without source code
EnsembleWeight=EnsembleWeight/mean(EnsembleWeight);
%EnsembleWeight=ones(size(EnsembleWeight));
[C0,Index]=mykmeans_c(single(Ensemble),K,single(EnsembleWeight)); Index=Index+1; C0=single(C0'); % my weighted kmeans with souce code


C=C0;meanC=mean(C);[sC,iC]=sort(meanC);C=C(:,iC);
Index0=zeros(1,size(Ensemble,1));
for i=1:size(Ensemble,1),Index0(i)=find(Index(i)==iC);end
Index=Index0;





if 1==L,
    mytree.centers=int32(C*2^0);

else
    mytree.centers=int32(C*2^0);

    for k=1:K,
        index=(Index==k);
        if sum(index)>=K,
            ensemble=Ensemble(index,:);
            ensembleweight=EnsembleWeight(index,:);
        else
            tempensemble=Ensemble(index,:);
            if isempty(tempensemble), tempensemble1=Ensemble(1,:) ;else tempensemble1=tempensemble(1,:);end
            tempdist=bsxfun(@minus,Ensemble,tempensemble1);
            tempdist=sum(tempdist.^2,2);
            [tempsorted,tempid]=sort(tempdist);
            ensemble=Ensemble(tempid(1:K),:);
            ensembleweight=EnsembleWeight(tempid(1:K),:);
        end
        if L>2, display(strcat('level(',num2str(L-1),'),cluster(',num2str(k),') clustered.'));end
        mytree.sub(k) = vl_myfillhikmeans_C(ensemble,ensembleweight,K, L-1);
    end
end

