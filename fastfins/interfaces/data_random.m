function y = data_random(model, obs, prior, v, ind, beta)
%generate data
%
% Tiangang Cui, 19/August/2019

% map to correlated Gaussian
u = matvec_prior_L(prior, v) + prior.mean_u;
y = zeros(numel(ind), size(u,2));
% forward solve
for i = 1:size(u,2)
    sol = forward_solve(model, u(:,i));
    switch obs.like
        case {'normal'}
            tmp = sol.d + randn(size(obs.data)).*(obs.std/sqrt(beta));
            y(:,i) = tmp(ind);
        otherwise
            error('likelihood not implemented')
    end
end

end
