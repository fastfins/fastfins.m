function [mlpts, mllkds] = minus_log_post_vec(model, obs, prior, vs)
%MODEL_SOLVE
%
% solve the forward model
% compute the minus log-likelihood, and minus log-posterior
% solve the adjoint gradient, optional
% assemble the information for evaluating the matvec with PPH, optional
%
% Tiangang Cui, 19/August/2019

% map to correlated Gaussian
us = matvec_prior_L(prior, vs) + prior.mean_u;
% forward solve
ds = forward_solve_vec(model, us);
%

switch obs.like
    case {'normal'}    
        misfit = (ds - obs.data(:))./obs.std(:);
        mllkds  = 0.5*sum(misfit(:).^2, 1); % minus log-likelihood
    case {'poisson'}
        mllkds = - obs.data(:)'*log(ds) + sum(ds,1); %- sum(log(factorial(obs.data(:))));
    otherwise
        error('likelihood not implemented')
end

mlpts  = mllkds + 0.5*sum(vs.^2, 1);   % minus log posterior

end