[numVec, distMat] = filereading_writing();
[~,seqs] = fastaread('snphylo.output.fasta');
distMat1 = seqpdist(seqs, 'Method', 'p-distance');
n = size(numVec, 1);
clustersNo = 14;
rng('default');



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
nDistMat = seqpdist(s, 'Method', 'p-distance');

CIndex = phylogenetictree(nDistMat, clustersNo);
%[useless]NJCIndex = phylogenetictreeNJ(nDistMat, clustersNo);
Silhouette = mean(silhouette([], CIndex, nDistMat));
%[useless]NJSilhouette = mean(silhouette([], NJCIndex, nDistMat));
nCIndex = zeros(n, 1);

for i = 1:n
    to = idx(i);
    nCIndex(i) = CIndex(to);
end

UPGMAsilhoutteVal = mean(silhouette([], nCIndex, distMat1));
%[useless]NJsilhoutteVal = mean(silhouette([], nCIndex, distMat1));