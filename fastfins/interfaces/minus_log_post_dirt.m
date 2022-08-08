function [mllkd, mlp, gmllkd] = minus_log_post_dirt(model, obs, prior, v)
%MODEL_SOLVE
%
% solve the forward model
% compute the minus log-likelihood, and minus log-posterior
% solve the adjoint gradient, optional
% assemble the information for evaluating the matvec with PPH, optional
%
% Tiangang Cui, 19/August/2019

mlp = 0.5*sum(v.^2, 1);
% map to correlated Gaussian
u  = matvec_prior_L(prior, v) + prior.mean_u;

n = size(u,2);

if nargout < 3
    mllkd = zeros(1,n);    
    for i = 1:n
        sol = forward_solve(model, u(:,i));
        mllkd(i) = minus_log_like(obs, sol.d);
    end
else
    mllkd = zeros(1,n);
    gmllkd = zeros(prior.dof,n);
    for i = 1:n
        sol = forward_solve(model, u(:,i));
        [mllkd(i),gd] = minus_log_like(obs, sol.d);
        gu = matvec_Jty(model, sol, gd);
        gmllkd(:,i) = matvec_prior_Lt(prior, gu);
    end
end

end
