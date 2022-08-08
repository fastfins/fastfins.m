function [mlpt, mllkd, gmlpt, gmllkd, sol] = minus_log_post_rare(model, obs, prior, c, sigma, v)
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
mlp = 0.5*sum(v.^2, 1);
% map to correlated Gaussian
u   = matvec_prior_L(prior, v) + prior.mean_u;
% forward solve
sol = forward_solve(model, u);
%
% for the rare event part
tmp = exp((c - sol.qoi)/sigma);
mlr = log(1 + tmp);
%
if nargout < 3
    mllkd   = minus_log_like(obs, sol.d) + mlr;
    mlpt    = mllkd + mlp;           % minus log posterior
else
    gmlr  = - tmp/((1+tmp)*sigma);
    gmlr  = gmlr * matvec_Qty(model, sol);
    %
    [mllkd,gd,sol.I] = minus_log_like(obs, sol.d);
    %
    mllkd   = mllkd + mlr;
    %
    mlpt    = mllkd + mlp;           % minus log posterior
    gu      = matvec_Jty(model, sol, gd);
    gmllkd  = matvec_prior_Lt(prior, gu+gmlr);
    gmlpt   = gmllkd + grad_p;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [mllkd, gmllkd, I] = minus_log_like(obs, d)
%
% minus log likelihood for a given observable model output 
% gradient w.r.t. observable model output
% Fisher information matrix at observable model output

switch obs.like
    case {'normal'}    
        misfit = (d - obs.data)./obs.std;
        mllkd  = 0.5*sum(reshape(misfit, obs.n_data, []).^2, 1); % minus log-likelihood
        gmllkd = misfit./obs.std;
        if isscalar(obs.std)
            I = speye(obs.n_data, obs.n_data)./obs.std^2;
        else
            I = spdiags(1./obs.std.^2, 0, obs.n_data, obs.n_data);
        end
    case {'poisson'}
        d = reshape(d, obs.n_data, []);
        mllkd = - obs.data(:)'*log(d) + sum(d, 1); %- sum(log(factorial(obs.data(:))));
        gmllkd = 1 - obs.data(:)./d;
        d(d<eps) = eps;
        I = spdiags(1./d(:), 0, obs.n_data*size(d,2), obs.n_data*size(d,2));
    otherwise
        error('likelihood not implemented')
end

end