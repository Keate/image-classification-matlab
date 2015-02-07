% An example to show how approximate  kernel kmeans  is invoked
% Generates two-dimensional Gaussian mixture Ensemble of size 500 and finds the 2
% clusters using 50 samples
clear,clc,close all

max_iter = 100;
k =10000;
N = 2e+4;
m = 1;
lambda = 1/100;

Ensemble = single(rand(N,128));
% for i=1:k
%     data=single(randn(ceil(N/k),128) + 5*(i-1));
%     Ensemble((i-1)*ceil(N/k)+1:i*ceil(N/k),:) = data;
% end

plot(Ensemble(:,1),Ensemble(:,2),'.b');

%Sample m Ensemble points
perm = single(randperm(N));
indices = perm(1:m);

[centers,mincenter,mindist,q2,quality] = fastkmeans(Ensemble,k,1);
%Compute m x N RBF kernel
Krect = single(exp(-lambda*(ones(m,1)*sum(Ensemble.^2,2)' - 2 * Ensemble(indices,:)*Ensemble' + sum(Ensemble(indices,:).^2,2)*ones(1,N))));
%Krect=rand(size(Krect));
labels=approx_kkmeans(Krect,k,max_iter,indices);
figure
plot(Ensemble(labels==1,1),Ensemble(labels==1,2),'.r',Ensemble(labels==2,1),Ensemble(labels==2,2),'.b',Ensemble(labels==3,1),Ensemble(labels==3,2),'.k',Ensemble(labels==4,1),Ensemble(labels==4,2),'.g')