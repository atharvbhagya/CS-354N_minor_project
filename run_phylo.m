function [CIndex, Silhoutte] = run_phylo(distMat, clustersNo)
    n = size(distMat, 2);
    M = squareform(distMat);
    Silhoutte = -2;
    dmat = sort(distMat);
    for i = 1:n
        
        
        clusteringUPGMA = phylogenetictree(distMat, clustersNo);
        nSilhoutte = mean(silhouette([], nCIndex, distMat));
        if (nSilhoutte > Silhoutte) 
            Silhoutte = nSilhoutte;
            CIndex = nCIndex;
        end
    end
end