function [mlpt, mllkd, gmlpt, gmllkd, sol] = minus_log_post(model, obs, prior, v, mask)
%MODEL_SOLVE
%
% solve the forward model
% compute the minus log-likelihood, and minus log-posterior
% solve the adjoint gradient, optional
% assemble the information for evaluating the matvec with PPH, optional
%
% Tiangang Cui, 19/August/2019

if nargin < 5
    mask = [];
end

% prior, assuming the input v is wsoltten parameters
grad_p  = v;
mlp = 0.5*sum(v.^2, 1);
% map to correlated Gaussian
u   = matvec_prior_L(prior, v) + prior.mean_u;
% forward solve
sol = forward_solve(model, u);
%
if nargout < 3
    mllkd   = minus_log_like(obs, sol.d, mask);
    mlpt    = mllkd + mlp;           % minus log posterior
else
    [mllkd,gd,sol.I] = minus_log_like(obs, sol.d, mask);
    mlpt    = mllkd + mlp;           % minus log posterior
    if isempty(mask)
        gu  = matvec_Jty(model, sol, gd);
    else
        gd_full = zeros(obs.n_data, 1);
        gd_full(mask) = gd;
        gu  = matvec_Jty(model, sol, gd_full);
    end
    gmllkd  = matvec_prior_Lt(prior, gu);
    gmlpt   = gmllkd + grad_p;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
