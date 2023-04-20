function [clusteringNJ , distMat] = phylogenetictreeNJ( distVect, clustersNo )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    NJtree = seqneighjoin(distVect);
    distMat = pdist(NJtree)';
    Z = linkage(distMat);
    clusteringNJ = cluster(Z, 'MaxClust', clustersNo);
end

