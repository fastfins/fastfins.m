function curr_p = rto_update_density(model, obs, prior, rto_lin, curr_p)

misfit          = (curr_p.d - obs.data)./obs.std; 
curr_p.llkd     = - 0.5*sum(misfit(:).^2);
curr_p.logd     = rto_func_d(model, obs, prior, rto_lin, curr_p);
curr_p.v        = matvec_prior_invL (prior, curr_p.u - prior.mean_u);
curr_p.lp       = - 0.5*sum(curr_p.v.^2, 1);

vr      = rto_lin.Phi'*curr_p.v;
tmp     = scale_rows(rto_lin.Psi'*misfit, rto_lin.s);
res     = scale_rows(tmp + vr, rto_lin.d);

curr_p.logrho   = 0.5* ( sum(res(:).^2) - sum(vr(:).^2) ) + curr_p.llkd - curr_p.logd + rto_lin.logc;

%0.5*sum(res(:).^2)
%0.5*sum(vr(:).^2) 

end