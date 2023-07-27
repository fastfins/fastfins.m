function [V, d] = eigen_PPHess_FD(model, obs, prior, v, fd_tol, tol, nmax)
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

[Vl,Dl] = eigs(@(dv)matvec_PPHess_FD(model,obs,prior,v,dv,fd_tol),prior.dof,nmax,...
    'largestreal','IsFunctionSymmetric',true);

dl = diag(Dl);
i  = abs(dl)>tol;
Vl = Vl(:,i);
dl = dl(i);

[Vs,Ds] = eigs(@(dv) matvec_PPHess_FD(model,obs,prior,v,dv,fd_tol),prior.dof,nmax,...
    'smallestreal','IsFunctionSymmetric',true);

ds = diag(Ds);
i  = abs(ds)>tol;
Vs = Vs(:,i);
ds = ds(i);

V = [Vl,Vs];
d = [dl;ds];

end

function w = matvec_PPHess_FD(model, obs, prior, v, dv, tol)
%MATVEC_PPH
%
% compute the matvec with PPH
%
% Tiangang Cui, 19/Mar/2014

if nargin == 5
    tol = 1E-5;
end

[~,~,~,gmllkd_p] = minus_log_post(model, obs, prior, v+tol*dv);
[~,~,~,gmllkd_m] = minus_log_post(model, obs, prior, v-tol*dv);

w = (gmllkd_p-gmllkd_m)/(2*tol);

end
