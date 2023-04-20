function accuracy = checkaccuracy(n, clustering, clusteringWPGMA)
    matched = 0;
    for i = 1:n
        h = i + 1;
        for j = h:n
            tx = 0;
            ty = 0;
            if clustering(i) == clustering(j)
                tx = 1;
            end
            if clusteringWPGMA(i) == clusteringWPGMA(j)
                ty = 1;
            end
            if tx == ty
                matched = matched + 1;
            end
        end
    end
    
    accuracy = (matched * 2) / (n * (n - 1));
    
end

