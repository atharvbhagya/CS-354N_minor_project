[numVec, distMat] = filereading_writing();
n = size(numVec, 1);
color = ['w', 'k', 'b', 'r', 'c', 'g', 'm', 'y'];
clustersNo = 13;
rng('default');
[idx,C] = kmedoids(numVec, n, 'Distance', 'hamming', 'replicates', 4);
nDistMat = pdist(C, 'hamming');

clusteringUPGMA = phylogenetictree(distMat, clustersNo);
[clusteringNJ1 , NJdistMat] = phylogenetictreeNJ(distMat, clustersNo);

[CIndex, ~] = Run(nDistMat, clustersNo);
nCIndex = zeros(n, 1);
for i = 1:n
    to = idx(i);
    nCIndex(i) = CIndex(to);
end

NJsilhoutteVal1 = mean(silhouette([], clusteringNJ1, NJdistMat));
SCSilhoutte = mean(silhouette([], nCIndex, distMat));
UPGMAsilhoutteVal = mean(silhouette([], clusteringUPGMA, distMat));

[idx,C] = kmedoids(numVec, round(n/2), 'Distance', 'hamming', 'replicates', 4);
nDistMat = pdist(C, 'hamming');
[CIndex, ~] = Run(nDistMat, clustersNo);

nCIndex = zeros(n, 1);
for i = 1:n
    to = idx(i);
    nCIndex(i) = CIndex(to);
end
SCVQSilhoutte = mean(silhouette([], nCIndex, distMat));