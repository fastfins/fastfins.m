function rto_lin = rto_local_lin(model, obs, prior, v, svd_tol, svd_nmax)
% compute the local linearization
%
% Tiangang Cui, August, 1, 2017

[~,~,~,~,HI] = minus_log_post(model, obs, prior, v);

if model.explicit_ja
    [rto_lin.Psi, rto_lin.s, rto_lin.Phi] ...
        = svd_explicit_WJ(model, obs, prior, v, svd_tol, svd_nmax);
else
    [rto_lin.Psi, rto_lin.s, rto_lin.Phi] ...
        = svd_rand_WJ(model, obs, prior, v, svd_tol, svd_nmax);
end

% linearisation
rto_lin.Fref    = HI.d./obs.std;
rto_lin.xref    = rto_lin.Phi'*v;
rto_lin.linref  = rto_lin.Psi.*rto_lin.s(:)';

% normalise data
rto_lin.data    = obs.data./obs.std;

rto_lin.d       = (rto_lin.s.^2 + 1).^(-0.5);
rto_lin.LPhi    = matvec_prior_L(prior, rto_lin.Phi);
rto_lin.type    = 'low_rank';
rto_lin.nrank   = length(rto_lin.s);

rto_lin.v       = v;
rto_lin.u       = matvec_prior_L (prior, v) + prior.mean_u;
rto_lin.mlp     = 0.5*sum(v.^2);
%{
rto_lin.logd    = sum(log(rto_lin.s.^2 + 1));
if length(obs.std) > 1
    rto_lin.logc    = 0.5*sum(log(rto_lin.s.^2 + 1)) - 0.5*sum(log(obs.std));
else
    rto_lin.logc    = 0.5*sum(log(rto_lin.s.^2 + 1)) - log(obs.std)*obs.n_data;
end
%}
rto_lin.logc    = -0.5*sum(log(rto_lin.s.^2 + 1));

rto_lin.dof     = prior.dof;

end