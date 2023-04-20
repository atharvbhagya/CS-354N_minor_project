function clusteringWPGMA = phylogenetictree(distVect, clustersNo)
    Z = linkage(distVect, 'UPGMA');
    clusteringWPGMA = cluster(Z, 'MaxClust', clustersNo);
end