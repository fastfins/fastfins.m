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