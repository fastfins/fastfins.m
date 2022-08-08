function [res, logrhos] = rto_density(model, obs, prior, rto_lin, us)
%GET_MAP_MATLAB
%
% Runs optimization algorithms to get the MAP estitimate
%
% Tiangang Cui, August, 1, 2017

m = size(us, 2);

res  = zeros(obs.n_data, m);
logrhos = zeros(m, 1);

for i = 1:m
    v           = matvec_prior_invL(prior, us(:,i) - prior.mean_u);
    vr          = rto_lin.Phi'*v;
    
    HI          = forward_solve(model, us(:,i));
    res(:,i)    = HI.d - obs.data;
    misfit      = res(:,i)./obs.std;   % G(v) - d
    llkd        = - 0.5*sum(misfit(:).^2);
    lp          = - 0.5*sum(vr.^2, 1);
    
    % RTO density 
    tmp         = rto_lin.linref'*misfit(:);
    Tvr         = rto_lin.d(:).*(tmp + vr);
    lrto        = - 0.5*sum(Tvr(:).^2);
    logd        = rto_func_d(model, obs, rto_lin, HI, false);
    
    logrhos(i) = lp + llkd - logd - lrto - rto_lin.logc;
end

end

