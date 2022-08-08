function [gsvd, P, S] = redu_init(hessian, options, v)
% Initialize the reduced subspace
% Tiangang Cui, 01/Oct/2013

% initialize

switch options.hess_type
    case {'Eig'}
        [gsvd.V,d_loc] = hessian(v, options.local_trunc_tol, options.local_max_dim);
        gsvd.S  = d_loc.^(0.5);
    case {'SVD'}
        [~,gsvd.S,gsvd.V] = hessian(v, options.local_trunc_tol, options.local_max_dim);
end

gsvd.dof = length(gsvd.S);
gsvd.n   = 1;

r2 = sum(gsvd.S>=options.lis_trunc_tol);
P  = gsvd.V(:,1:r2);
S  = gsvd.S(1:r2);

end