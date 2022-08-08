function [V, d, sol] = eigen_PPFisher(model, obs, prior, sol, tol, nmax)
%EIGEN_PPGNH
%
% Eigendecomposition of the prior-preconditioned Hessian
% requires the HI by runing
% [~,~,~, HI] = minus_log_post(model, prior, v);
%
% tol:       is the truncation tolerance
% max_eigen: is the maximum number of modes solved by eigs
%
% Tiangang Cui, 19/Mar/2014

if ~isstruct(sol)
    [~,~,~,~,sol] = minus_log_post(model, obs, prior, sol);
end
    
if obs.jacobian
    Ju = explicit_jacobian(model, sol);
    Jv = matvec_prior_Lt(prior, Ju')';
    H  = Jv'*sol.I*Jv;
    
    [V,D]   = eig(H);
    [d,ind] = sort(diag(D), 'descend');
    
    r = min(nmax, sum(d>=tol));
    V = V(:,ind(1:r));
    d = d(1:r);
else
    opts.issym  = 1;
    opts.isreal = 1;
    
    [V, D]  = eigs(@(dv)  matvec_PPFisher(model, prior, sol, dv), prior.dof,  nmax, 'LA', opts);
    d = diag(D);
    i = d>tol;
    V = V(:,i);
    d = d(i);
end

end