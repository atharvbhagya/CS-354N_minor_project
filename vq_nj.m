[numVec, distMat] = filereading_writing();
[~,seqs] = fastaread('snphylo.output.fasta');
n = size(numVec, 1);

%distMat1 = pdist(B, 'hamming');
distMat1 = seqpdist(seqs, 'Method', 'p-distance');

clustersNo = 16;
rng('default');
%alignment-score
%p-distance
%Jukes-Cantor
[idx,C] = kmedoids(numVec, round(n/2), 'Distance', 'hamming', 'replicates', 4);
s(:,1)=seqs(:,13);
 s(:,2)=seqs(:,20);
 s(:,3)=seqs(:,7);
 s(:,4)=seqs(:,10);
 s(:,5)=seqs(:,17);
 s(:,6)=seqs(:,11);
 s(:,7)=seqs(:,9);
 s(:,8)=seqs(:,15);
 s(:,9)=seqs(:,31);
 s(:,10)=seqs(:,8);
 s(:,11)=seqs(:,6);
 s(:,12)=seqs(:,16);
 s(:,13)=seqs(:,25);
 s(:,14)=seqs(:,12);
 s(:,15)=seqs(:,21);
 s(:,16)=seqs(:,5);
%nDistMat = pdist(C, 'hamming');
nDistMat = seqpdist(s, 'Method', 'Jukes-Cantor');
[CIndex , xyz] = phylogenetictreeNJ(nDistMat, clustersNo);
Silhoutte = mean(silhouette([], CIndex, xyz));
nCIndex = zeros(n, 1);

%[clusteringNJ1 , NJdistMat] = phylogenetictreeNJ(nDistMat, clustersNo);
%nDistMat = pdist(C, 'hamming');
% %-- Distance matrix for all 31 seq--
% DistMat = pdist(numVec, 'hamming');
% NJtree = seqneighjoin(DistMat);
%     distmat = pdist(NJtree)';

for i = 1:n
    to = idx(i);
    nCIndex(i) = CIndex(to);
end

NJsilhoutteVal = mean(silhouette([], nCIndex, distMat1));
%NJsilhoutteVal = mean(silhouette([], nCIndex, distMat1));