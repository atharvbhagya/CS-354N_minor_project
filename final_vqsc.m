[~,seqs] = fastaread('snphylo.output.fasta');
distMat1 = seqpdist(seqs, 'Method', 'alignment-score');
  rng('default');  
n = size(seqs, 2);
m = size(seqs{1}, 2);
numVec = zeros(n, m);
    for i = 1:n
        cur = seqs{i};
        for j = 1:m
            if cur(j) == 'A'
                numVec(i, j) = 0;
            elseif cur(j) == 'T'
                numVec(i, j) = 1;
            elseif cur(j) == 'G'
                numVec(i, j) = 2;
            elseif cur(j) == 'C'
                numVec(i, j) = 3;
            elseif cur(j) == 'N'
                numVec(i, j) = 4;
            elseif cur(j) == 'M'
                numVec(i, j) = 5;
            elseif cur(j) == 'Y'
                numVec(i, j) = 6;
            elseif cur(j) == 'R'
                numVec(i, j) = 7;
            elseif cur(j) == 'S'
                numVec(i, j) = 8;
            elseif cur(j) == 'W'
                numVec(i, j) = 9;
            end
        end
    end
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
 
 distMat = seqpdist(s, 'Method', 'alignment-score');
clustersNo = 13;
rng('default');
Silhoutte = -2;
M = squareform(distMat);
clusteringType = 2;
dmat = sort(distMat);
    for i = 1:n
        sigma = dmat(i);
        simGraph = exp(-M.^2 ./ (2*sigma^2));
        nCIndex = SpectralClustering(simGraph, clustersNo, clusteringType);
        nSilhoutte = mean(silhouette([], nCIndex, distMat));
        if (nSilhoutte > Silhoutte) 
            Silhoutte = nSilhoutte;
            CIndex = nCIndex;
        end
    end
    
    nCIndex = zeros(n, 1);
for i = 1:n
    to = idx(i);
    nCIndex(i) = CIndex(to);
end
SCVQSilhoutte = mean(silhouette([], nCIndex, distMat1));

            