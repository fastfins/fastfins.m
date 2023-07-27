function [V, d, sol] = eigen_PPHess(model, obs, prior, sol, tol, nmax)
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

%opts.issym  = 1;
%opts.isreal = 1;
%[V, D]  = eigs(@(dv)  matvec_PPHess(model, prior, sol, dv), prior.dof,  nmax,'LA', opts);

[Vl,Dl] = eigs(@(dv)matvec_PPHess(model,prior,sol,dv),prior.dof,nmax,...
    'largestreal','IsFunctionSymmetric',true);

dl = diag(Dl);
i  = abs(dl)>tol;
Vl = Vl(:,i);
dl = dl(i);

%[Vs,Ds] = eigs(@(dv)matvec_PPHess(model,prior,sol,dv),prior.dof,nmax,...
%    'smallestreal','IsFunctionSymmetric',true);
[Vs,Ds] = eigs(@(dv) matvec_PPHess(model,prior,sol,dv),prior.dof,nmax,...
    'smallestreal','IsFunctionSymmetric',true);

ds = diag(Ds);
i  = abs(ds)>tol;
Vs = Vs(:,i);
ds = ds(i);

V = [Vl,Vs];
d = [dl;ds];

end