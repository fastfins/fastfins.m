function [U, s, V, Jv] = svd_explicit_WJ(model, obs, prior, sol, tol, nmax)
%SVD_EXPLICIT_WJ
%
% SVD for factorizing the forward model
% requires the sol by runing 
% [~,~,~, sol] = minus_log_post(model, prior, v);
%
% Tiangang Cui, 17/Jan/2014

if ~isstruct(sol)
    [~,~,~,~,sol] = minus_log_post(model, obs, prior, sol);
end

Ju = explicit_J(model, sol);
Ju = Ju./obs.std(:);
Jv = matvec_prior_Lt(prior, Ju')';

[U,S,V] = svd(Jv,'econ');

[dS,ind] = sort(diag(S), 'descend');
r = min(nmax, sum(dS>=tol));
U = U(:,ind(1:r));
V = V(:,ind(1:r));
s = dS(1:r);


end