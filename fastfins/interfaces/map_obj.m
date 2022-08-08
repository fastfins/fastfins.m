function [sol,g] = map_obj(model, obs, prior, v)
%MODEL_SOLVE
%
% solve the forward model
% compute the minus log-likelihood, and minus log-posterior
% solve the adjoint gradient, optional
% assemble the information for evaluating the matvec with PPH, optional
%
% Tiangang Cui, 19/August/2019

% prior, assuming the input v is wsoltten parameters
grad_p  = v;
mlp = 0.5*(v'*grad_p);
% map to correlated Gaussian
u   = matvec_prior_L(prior, v) + prior.mean_u;
% forward solve
sol = forward_solve(model, u);
%
if nargout == 1
    mllkd   = minus_log_like(obs, sol.d);
    sol.f   = mllkd + mlp;           % minus log posterior
else
    [mllkd,gd,sol.I] = minus_log_like(obs, sol.d);
    sol.f   = mllkd + mlp;           % minus log posterior
    gu      = matvec_Jty(model, sol, gd);
    gmllkd  = matvec_prior_Lt(prior, gu);
    g       = gmllkd + grad_p;
end

end