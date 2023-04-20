function [numVec, distMat] = filereading_writing()
    [~,seqs] = fastaread('snphylo.output.fasta');
    distMat = seqpdist(seqs, 'Method', 'p-distance');
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
end