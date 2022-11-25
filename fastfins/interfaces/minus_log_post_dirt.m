function [mllkd, mlp, gmllkd] = minus_log_post_dirt(model, obs, prior, v, masks)
%MODEL_SOLVE
%
% solve the forward model
% compute the minus log-likelihood, and minus log-posterior
% solve the adjoint gradient, optional
% assemble the information for evaluating the matvec with PPH, optional
%
% Tiangang Cui, 19/August/2019

if nargin < 5
    masks = [];
end

mlp = 0.5*sum(v.^2, 1);
% map to correlated Gaussian
u  = matvec_prior_L(prior, v) + prior.mean_u;

n = size(u,2);

if iscell(masks)
    if nargout < 3
        mllkd = cell(size(masks));
        for mm = 1:length(masks)
            mllkd{mm} = zeros(1,n);
            for i = 1:n
                sol = forward_solve(model, u(:,i));
                mllkd{mm}(i) = minus_log_like(obs, sol.d, masks{mm});
            end
        end
    else
        mllkd = cell(size(masks));
        gmllkd = cell(size(masks));
        for mm = 1:length(masks)
            mllkd{mm} = zeros(1,n);
            gmllkd{mm} = zeros(prior.dof,n);
            for i = 1:n
                sol = forward_solve(model, u(:,i));
                [mllkd{mm}(i),gd] = minus_log_like(obs, sol.d, masks{mm});
                gd_full = zeros(obs.n_data, 1);
                gd_full(masks{mm}) = gd;
                gu  = matvec_Jty(model, sol, gd_full);
                gmllkd{mm}(:,i) = matvec_prior_Lt(prior, gu);
            end
        end
    end
else
    if nargout < 3
        mllkd = zeros(1,n);
        for i = 1:n
            sol = forward_solve(model, u(:,i));
            mllkd(i) = minus_log_like(obs, sol.d, masks);
        end
    else
        mllkd = zeros(1,n);
        gmllkd = zeros(prior.dof,n);
        for i = 1:n
            sol = forward_solve(model, u(:,i));
            [mllkd(i),gd] = minus_log_like(obs, sol.d, masks);
            if isempty(masks)
                gu  = matvec_Jty(model, sol, gd);
            else
                gd_full = zeros(obs.n_data, 1);
                gd_full(masks) = gd;
                gu  = matvec_Jty(model, sol, gd_full);
            end
            gmllkd(:,i) = matvec_prior_Lt(prior, gu);
        end
    end
end

end
