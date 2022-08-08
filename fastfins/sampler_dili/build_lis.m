function [P,S,gsvd,out] = build_lis(ml_target, hessian, vmap, options)
%BUILD_LIS
%
% builds likelihood informed subspace 
%
% Tiangang Cui, 17/Jan/2019

switch options.method
    case{'Prior','Laplace'}
        [P,S,gsvd] = build_lis_reference(options.method, hessian, vmap, options);
        out = [];
    case{'DILI'} 
        options = mcmc_options(options, 'proposal', 'MALA', 'using_gibbs', true);
        options.nstep = options.max_hess*options.nbatch;
        [P,S,gsvd,out] = build_lis_dili(ml_target, hessian, vmap, options);
    case{'MALA'}
        error('not implemented')
        %[P,S,gsvd,out] = build_lis_mala(ml_target, hessian, vmap, options);
end

end
