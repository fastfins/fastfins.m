function [mlpt, mllkd, gmlpt, gmllkd, sol] = minus_log_post_temp(model, obs, prior, v, b1, b2)
%MODEL_SOLVE
%
% solve the forward model
% compute the minus log-likelihood, and minus log-posterior
% solve the adjoint gradient, optional
% assemble the information for evaluating the matvec with PPH, optional
%
% Tiangang Cui, 19/August/2019

pp = 0.01; % Prior temperature power

% prior, assuming the input v is wsoltten parameters
grad_p  = v;
mlp = 0.5*sum(v.^2, 1);

if (nargout==1)
    % Parallelise this stuff
    % map to correlated Gaussian
    u   = matvec_prior_L(prior, v) + prior.mean_u;
    mllkd = mlp;
    parfor i=1:size(v,2)
        % forward solve
        soli = forward_solve(model, u(:,i));
        mllkd(i)   = minus_log_like(obs, soli.d);
    end
    mlpt   = mllkd*(b2 - b1) + mlp*(b2^pp - b1^pp);           % minus log posterior
    
    return;
end


% map to correlated Gaussian
u   = matvec_prior_L(prior, v) + prior.mean_u;
% forward solve
sol = forward_solve(model, u);
%
if nargout < 3
    mllkd   = minus_log_like(obs, sol.d);
    mlpt    = mllkd*(b2 - b1) + mlp*(b2^pp - b1^pp);           % minus log posterior
else
    [mllkd,gd,sol.I] = minus_log_like(obs, sol.d);
    mlpt    = mllkd*(b2 - b1) + mlp*(b2^pp - b1^pp);           % minus log posterior
    gu      = matvec_Jty(model, sol, gd);
    gmllkd  = matvec_prior_Lt(prior, gu);
    gmlpt   = gmllkd*(b2 - b1) + grad_p*(b2^pp - b1^pp);
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