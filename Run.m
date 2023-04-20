function [CIndex, Silhoutte] = Run(distMat, clustersNo)
    n = size(distMat, 2);
    M = squareform(distMat);
    clusteringType = 2;
    Silhoutte = -2;
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
end