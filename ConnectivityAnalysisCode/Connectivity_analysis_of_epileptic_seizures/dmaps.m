function [ psi ] = dmaps( data, neigs, ep_sc )
%DMAPS Calculates the diffusion maps eigenvectors (for dimensionality
% reduction.
% Input:    data    -   data matrix of size: (samples X features)
%           neigs   -   number of eigenvectors to return (dimension of the
%                       new representation
%           ep_sc   -   kernel scale parameter - multiplying the mean
% Output:   psi     -   new representation vectors (samples X neigs)

dis = squareform(pdist(data));
ep  = ep_sc*median(dis(:));
W   = exp(-dis.^2./(ep.^2));
K   = diag(1./sum(W,1)) * W;

[v, e]  = eigs(K,neigs+1);
[~, ie] = sort(diag(e),'descend');
psi     = v(:,ie(2:end));

end

