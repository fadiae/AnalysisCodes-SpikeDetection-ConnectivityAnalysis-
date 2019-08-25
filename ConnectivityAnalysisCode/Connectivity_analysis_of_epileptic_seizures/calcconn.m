function [ conn ] = calcconn( data, dmWin )
%CALCCONN Calculates the connectivity matrix for ``data'' using a time-lag
% window of ``dmWin'', based on the Euclidean distance between the 
% stationary distributions of the kernels.
% Input:    data  -    time X features (e.g. ROIs).
%           dmWin -    time lag window to use.

N  = size(data,2);
sD = cell(1,N);

for s = 1:N
    
    % Computing the lag-map of each feature (for dmWin>1):
    if dmWin>1
        lagdata = lagmatrix(data(:,s),0:dmWin);
        lagdata = lagdata((dmWin+1):end,:);
    else
        lagdata = data(:,s);
    end
    
    % Computing the distance matrix:
    dist = squareform(pdist(lagdata));
    
    % Computing the kernel scale:
    if median(dist(:)) ~= 0
        eps = median(dist(:));
    else
        eps = mean(dist(:));
    end
    
    % Computing the normalized affinity kernel (Markovian matrix):
    K = exp(-dist.^2/(eps^2));
    K = K * diag(1./sum(K,1));
        
    % Computing the stationary distributions of each matrix (1st left
    % eigenvector):
    [stDist, ~] = eigs(K,1);
    stDist      = sign(stDist(1))*stDist;
    sD{s}       = stDist/sum(stDist);

end

% Calculating the connectivity matrix based on the distances between the
% stationary distributions:
conn    = zeros(N,N);
cmetric = @(p1,p2) sqrt(sum((p2 - p1).^2));
for row = 1:N
    for col = 1:(row-1)
        conn(row,col) = cmetric(sD{row},sD{col});
    end
end
conn = conn + conn.';

% Creating an affinity matrix from the distances:
conn = exp(-conn/(4*median(conn(:))));

end