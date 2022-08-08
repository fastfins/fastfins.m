function logrho = rto_density(model, obs, prior, rto_lin, u)
%GET_MAP_MATLAB
%
% Evaluate RTO density for a given v
%
% Tiangang Cui, August, 1, 2017


[x, dxdu]   = u2x( prior, u );
HI          = forward_solve(model, x);
HI.dxdu     = dxdu;

misfit  = (HI.d - obs.data)./obs.std; 
llkd    = - 0.5*sum(misfit(:).^2);
vr      = rto_lin.Phi'*matvec_prior_invL (prior, u - prior.mean_u);

tmp     = scale_rows(rto_lin.Psi'*misfit, rto_lin.s);
res     = scale_rows(tmp + vr, rto_lin.d);

logd    = rto_func_d(model, obs, prior, rto_lin, HI);

logrho  = 0.5*sum(res(:).^2) - 0.5*sum(vr(:).^2) + llkd - logd + rto_lin.logc;

end