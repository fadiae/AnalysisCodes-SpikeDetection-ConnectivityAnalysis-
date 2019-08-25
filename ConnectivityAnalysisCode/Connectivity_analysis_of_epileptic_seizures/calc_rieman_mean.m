function [ cref ] = calc_rieman_mean( mats, ep, maxiter )
%CALC_RIEMAN_MEAN Calculates the Riemannian mean of the matrices ``mats''
% iteratively.
% Input:    mats    -    SPD matrices array of size (N X N X number of matrices)
%           ep      -    threshold on the difference between iterations in
%                        the calculation of the mean.
%           maxiter -    maximum number of iterations to perform in the
%                        calculation of the mean.
% Output:   cref    -    A matrix which is the calculated Riemannian mean
%                        of ``mats''.

% Initialization:
cref = mean(mats,3);

% Calculating the mean iteratively:
for ite = 1:maxiter
    csqrt       = sqrtm(cref);
    csqrtn      = csqrt^(-1);
    logCp       = zeros(size(mats,1));
    for kk = 1:size(mats,3)
        logCp = logCp + csqrt * logm(csqrtn * mats(:,:,kk) * csqrtn) * csqrt;
    end
    sbar = logCp / size(mats,3);
    cref = csqrt * expm(csqrtn * sbar * csqrtn) * csqrt;
    
    if norm(sbar,'fro') < ep
        break;
    end
end




end

