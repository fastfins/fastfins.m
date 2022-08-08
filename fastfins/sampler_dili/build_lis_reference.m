function [P,S,gsvd] = build_lis_reference(method, hessian, vmap, options)

np  = length(vmap);

switch method
    case{'Laplace'}
        lap = get_laplace(@(x) hessian(x, options.local_trunc_tol, options.local_max_dim), vmap);
        v   = vmap;
    case{'Prior'}
        v = randn(np,1);
end

ni = 10;
ii = 0;
num_basis   = zeros(options.max_hess, 1);
param_dist  = zeros(options.max_hess,1);
[gsvd,P,S]  = redu_init(hessian, options, v);

for i = 2:options.max_hess
    % random number
    switch method
        case{'Laplace'}
            r = randn(np,1);
            v = vmap + r + lap.V_ihalf*(lap.V'*r);
        case{'Prior'}
            v = randn(np,1);
    end

    P1 = P;
    S1 = S;
    % update lis, re initilaize the MCMC
    [gsvd,P,S] = redu_enrich(hessian, options, gsvd, v);
    
    % check for convergence of the global basis
    ii              = ii + 1;
    d               = dist_fm(P1, S1, P, S);
    param_dist(ii)  = d;
    num_basis(ii)   = size(P,2);
    
    if ii < ni
        md = sum(param_dist(1:ii))/ii;
    else
        md = sum(param_dist(ii-ni+1:ii))/ni;
    end
    
    % recompute the tolerance
    if ii == ni
        options.global_conv_tol = md * options.global_conv_tol;
        fprintf('Convergence tol: %E\n\n', options.global_conv_tol);
    end
    fprintf('%5d%5d%5d\t%E\t%E\n', [gsvd.n, gsvd.dof, length(S), d, md]);
    
    if md < options.global_conv_tol && gsvd.n > options.min_hess
        break;
    end
end

end