function [ feat ] = calc_rieman_feat( mats, ep, maxiter )
%CALC_RIEMAN_FEAT Converts SPD matrices to features in a Euclidean space, by: 
% (1) finding the Riemannian mean of the matrices ``mats'', denoted by
%     ``Cref''
% (2) projecting all the matrices to the tangent space (of the Riemannian 
%     manifold of SPD matrices) at ``Cref''
% Input:    mats    -    SPD matrices array of size (N X N X number of matrices)
%           ep      -    threshold on the difference between iterations in
%                        the calculation of the mean.
%           maxiter -    maximum number of iterations to perform in the
%                        calculation of the mean.
% Output:   feat    -    features on the tangent space (Euclidean space),
%                        reprenting the matrices in ``mats''.

N       = size(mats,1);
fN      = N * ( N +1 ) /2;
featNrm = sqrt(2)*ones(N) + (1 - sqrt(2)) * eye(N);
feat    = zeros(size(mats,3), fN);

% Calculating the Riemannian mean of the matrices:
cref    = calc_rieman_mean( mats, ep, maxiter );
crefsqt = sqrtm(cref)^(-1);

% Projecting the matrices onto the tangent space at cref and creating the
% Euclidean features:
for kk = 1:size(mats,3)
    tmpMat     = logm(crefsqt * mats(:,:,kk) * crefsqt) .* featNrm;
    feat(kk,:) = tmpMat(triu(true(size(tmpMat))));
end

end

