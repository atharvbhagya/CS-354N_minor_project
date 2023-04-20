no_cluster = 13;
[~,seqs] = fastaread('snphylo.output.fasta');
distances = seqpdist(seqs,'Method','Jukes-Cantor');

upgma_tree = seqlinkage(distances,'UPGMA',seqs);
distMat_upgma = pdist(upgma_tree)';
Z = linkage(distMat_upgma);
index = cluster(Z, 'MaxClust', no_cluster);

nj_tree = seqneighjoin(distances,'equivar',seqs);
distMat_nj = pdist(nj_tree)';
z=linkage(distMat_nj);
index1=cluster(z, 'MaxClust', no_cluster);

nj_Silhoutte = mean(silhouette([], index1, distMat_nj));
upgma_Silhouette = mean(silhouette([], index1, distMat_upgma));