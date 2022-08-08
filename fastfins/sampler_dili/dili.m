function [out, stat] = dili(ml_target, P, stat, init, options)
%DILI_MCMC
%
% Dimension invariant, likelihood informed (DILI) MCMC sampler
% Tiangang Cui, 01/August/2019
%
%
% LIS:          the likelihood induced subspace
% CS:           the complement of the likelihood induced subspace

curr.v  = init;
np      = length(init);
sigma   = options.sigma;
dt      = options.dt;
out     = mcmc_outputs(np, options);
out.sub_samples = zeros(size(P,2),size(out.samples,2));

% initialise the mean and covariance
if ~isfield(stat,'M') || ~isfield(stat,'C') || ~isfield(stat,'n')
    stat.n     = 1;
    stat.sum   = P'*init(:);
    stat.cross = stat.sum(:)*stat.sum(:)';
else
    stat.sum   = stat.M*stat.n;
    stat.cross = stat.C*stat.n + stat.sum(:)*stat.sum(:)'/stat.n;
end
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

for i = 1:options.nstep
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
    
    if  options.adapt    
        batch = batch + 1;
        if  batch == options.nbatch
            delta = min(0.1,sqrt(options.nbatch/i));
            if (acc/options.nbatch) < options.rate
                sigma = sigma - delta;
            else
                sigma = sigma + delta;
            end
            batch = 0;
            acc   = 0;
        end
        [kernel,stat] = build_kernel(options.proposal, stat, sigma);
    end
    
    if  options.sbatch > 1
        if  mod(i, options.sbatch) == 0
            out.j = out.j + 1;
            out.samples(:,out.j)     = curr.v;
            out.sub_samples(:,out.j) = curr.v_sub;
            out.mlpt(out.j)  = curr.mlpt;
            out.mllkd(out.j) = curr.mllkd;
            out.mh(out.j)    = alpha;
            out.sigma(out.j) = sigma;
        end
    else
        out.samples(:,i)     = curr.v;
        out.sub_samples(:,i) = curr.v_sub;
        out.mlpt(i)  = curr.mlpt;
        out.mllkd(i) = curr.mllkd;
        out.mh(i)    = alpha;
        out.sigma(i) = sigma;
    end
end

out.acc_rate = [acc_s, acc_n]/options.nstep;

end

