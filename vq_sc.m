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

x = ismember(numVec,C,'rows');


j=1;
for i=1:31
    if x(i)==1
    sequence(j,:)=seqs(:,i);
    j=j+1;
    end
end
    
sequence = sequence';

distMat = seqpdist(sequence, 'Method', 'alignment-score');
clustersNo = 11;
rng('default');
Silhoutte = -2;
M = squareform(distMat);
clusteringType = 3;
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
