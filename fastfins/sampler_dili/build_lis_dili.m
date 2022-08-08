function [P,S,gsvd,out] = build_lis_dili(ml_target, hessian, vmap, options)
%BUILD_LIS_DILI
%
% Build the LIS using DILI
% Tiangang Cui, 01/August/2019
%

curr.v  = vmap;
np      = length(vmap);
sigma   = options.sigma;
dt      = options.dt;
out     = mcmc_outputs(np, options);

% initialise lis, mean and covariance
[gsvd,P,S] = redu_init(hessian, options, vmap);
tmp     = (S.^2 + 1).^(-1);
stat.n  = gsvd.dof;
stat.sum   = (P'*vmap)*stat.n;
stat.cross = diag(tmp)*stat.n + stat.sum*stat.sum'/stat.n;
stat.ref   = P'*vmap;
% end of initialisation

[kernel,stat] = build_kernel(options.proposal, stat, sigma);

switch options.proposal
    case {'MALA'}
        [curr.mlpt, curr.mllkd, curr.grad] = ml_target(curr.v);
        curr.grad_v_sub = P'*curr.grad;
    otherwise
        [curr.mlpt, curr.mllkd] = ml_target(curr.v);
end
curr.v_sub  = P'*curr.v;
curr.v_null = curr.v - P*curr.v_sub;

% MCMC
acc_s = 0;
acc_n = 0;
acc   = 0;
batch = 0;

ni = 10;
ii = 0;
num_basis   = zeros(options.max_hess, 1);
param_dist  = zeros(options.max_hess,1);

nstep = options.max_hess*options.nbatch;
stop_iter = nstep;
for i = 1:nstep
    % random number
    r = randn(np,1);
    
    if  options.using_gibbs % gibbs update
        % project random number
        r_sub  = P'*r;
        r_null = r - P*r_sub;
        
        % propose in the subspace, and evaluate the acceptance rate
        [alpha, next] = propose_sub(options.proposal, ml_target, P, kernel, curr, r_sub);
        if  log(rand) < alpha
            acc   = acc + 1;
            acc_s = acc_s + 1;
            curr  = next; % update
        end
        stat.cross = stat.cross + curr.v_sub(:)*curr.v_sub(:)';
        stat.sum   = stat.sum + curr.v_sub;
        stat.n     = stat.n + 1;
        
        % propose in the null, and evaluate the acceptance rate
        [alpha, next] = propose_null(options.proposal, ml_target, P, dt, curr, r_null);
        if  log(rand) < alpha
            acc_n = acc_n + 1;
            curr  = next; % update
        end
    else % all in one update
        % propose in the null, and evaluate the acceptance rate
        [alpha, next] = propose(options.proposal, ml_target, P, kernel, dt, curr, r);
        if  log(rand) < alpha
            acc   = acc + 1;
            acc_s = acc_s + 1;
            curr  = next; % update
        end
        stat.cross = stat.cross + curr.v_sub(:)*curr.v_sub(:)';
        stat.sum   = stat.sum + curr.v_sub;
        stat.n     = stat.n + 1;
    end % end of update
      
    batch = batch + 1;
    if  batch == options.nbatch
        if  options.adapt
            delta = min(0.1,sqrt(options.nbatch/i));
            if (acc/options.nbatch) < options.rate
                sigma = sigma - delta;
            else
                sigma = sigma + delta;
            end
        end
        acc   = 0;
        batch = 0;
        
        P1 = P;
        S1 = S;
        % update lis, re initilaize the MCMC
        [gsvd,P,S] = redu_enrich(hessian, options, gsvd, curr.v);
        if  strcmp(options.proposal, 'MALA')
            curr.grad_v_sub = P'*curr.grad;
        end
        curr.v_sub  = P'*curr.v;
        curr.v_null = curr.v - P*curr.v_sub;
        
        if  options.sbatch > 1
            tmp = P'*out.samples(:,1:out.j);
            stat.n = out.j;
        else
            tmp = P'*out.samples(:,i);
            stat.n = i;
        end
        stat.sum    = sum(tmp, 2);
        stat.cross  = tmp*tmp';
        stat.ref    = P'*vmap;
        
        % rebuild kernel
        [kernel,stat] = build_kernel(options.proposal, stat, sigma);
        
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
            stop_iter = i;
            break;
        end
    end
    
    if  options.sbatch > 1
        if  mod(i, options.sbatch) == 0
            out.j = out.j + 1;
            out.samples(:,out.j) = curr.v;
            out.mlpt(out.j)  = curr.mlpt;
            out.mllkd(out.j) = curr.mllkd;
            out.mh(out.j)    = alpha;
            out.sigma(out.j) = sigma;
        end
    else
        out.samples(:,i) = curr.v;
        out.mlpt(i)  = curr.mlpt;
        out.mllkd(i) = curr.mllkd;
        out.mh(i)    = alpha;
        out.sigma(i) = sigma;
    end
end

out.stat = stat;
out.acc_rate = [acc_s, acc_n]/min(nstep, stop_iter);

end
