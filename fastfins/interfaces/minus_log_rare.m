function [mlr, gmlr, sol] = minus_log_rare(model, prior, c, sigma, v)
%MODEL_SOLVE
%
% solve the forward model
% compute the minus log-likelihood, and minus log-posterior
% solve the adjoint gradient, optional
% assemble the information for evaluating the matvec with PPH, optional
%
% Tiangang Cui, 19/August/2019

% map to correlated Gaussian
u   = matvec_prior_L(prior, v) + prior.mean_u;
% forward solve
sol = forward_solve(model, u);
%
% f = 1/(1 + exp((c - x)/sigma));
% log(f)   = - log( 1 + exp((c - x)/sigma) )
% - log(f) = log( 1 + exp((c - x)/sigma) )

tmp = exp((c - sol.qoi)/sigma);
mlr = log(1 + tmp);

if nargout > 1
    gmlr  = - tmp/((1+tmp)*sigma);
    gmlr  = gmlr * matvec_Qty(model, sol);
    gmlr  = matvec_prior_Lt(prior, gmlr);
end

end
