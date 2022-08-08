function [out, stat] = amcmc(ml_target, stat, init, options)
%AMCMC
%
% Adaptive MCMC
%
% Tiangang Cui, 07/July/2019

v_curr  = init;
np      = length(init);
sigma   = options.sigma;
stat    = init_cov(options.proposal, stat, sigma, init);
out     = mcmc_outputs(np, options);

switch options.proposal
    case {'MALA'}
        [mlpt_curr, mllkd_curr, mg_curr] = ml_target(v_curr);
    case {'RW'}
        [mlpt_curr, mllkd_curr] = ml_target(v_curr);
        mg_curr = 0;
    case {'OW_Prior'}
        [mlpt_curr, mllkd_curr] = ml_target(v_curr);
        mg_curr = 0;
end

% start MCMC
acc_t = 0;
acc   = 0;
batch = 0;

for i = 1:options.nstep
    switch options.proposal
        case {'MALA'}
            drift_curr  = -(stat.C*mg_curr)*(stat.scale^2/2);
            r_curr  = randn(np,1);
            v_next  = v_curr + drift_curr + stat.scale*(stat.L*r_curr);
            
            [mlpt_next, mllkd_next, mg_next] = ml_target(v_next);
            
            % drift_next is not re-used, because the adaptation
            drift_next  = -(stat.C*mg_next)*(stat.scale^2/2);
            r_next  = -( (stat.L\(drift_curr + drift_next))/stat.scale + r_curr );
            log_n2c = - 0.5 * r_next(:)'*r_next(:);
            log_c2n = - 0.5 * r_curr(:)'*r_curr(:);
            alpha   = (mlpt_curr - mlpt_next) + (log_n2c - log_c2n);
            
        case {'RW'}
            v_next  = v_curr + stat.scale*(stat.L*randn(np,1));
            [mlpt_next, mllkd_next] = ml_target(v_next);
            mg_next = 0;
            alpha   = mlpt_curr - mlpt_next;
            
        case {'OW_Prior'}
            r       = randn(np,1);
            v_next  = stat.A*v_curr + stat.B*r;
            % eval density
            [mlpt_next, mllkd_next] = ml_target(v_next);
            mg_next = 0;
            alpha   = mllkd_curr - mllkd_next;
    end
    
    if  log(rand) < alpha
        v_curr      = v_next;
        mg_curr     = mg_next;
        mlpt_curr   = mlpt_next;
        mllkd_curr  = mllkd_next;
        acc         = acc+1;
        acc_t       = acc_t+1;
    end
    
    if  options.adapt
        batch = batch + 1;
        % update covariance
        stat.cross = stat.cross + v_curr(:)*v_curr(:)';
        stat.sum   = stat.sum + v_curr(:);
        stat.n     = stat.n + 1;
        
        if  batch == options.nbatch
            delta = min(0.1,sqrt(options.nbatch/i));
            if (acc/options.nbatch) < options.rate
                sigma = sigma - delta;
            else
                sigma = sigma + delta;
            end
            batch = 0;
            acc   = 0;
            stat  = update_cov(options.proposal, stat, sigma);
        end
    end
    
    if  options.sbatch > 1
        if mod(i, options.sbatch) == 0
            out.j = out.j + 1;
            out.samples(:,out.j) = v_curr;
            out.mlpt(out.j)  = mlpt_curr;
            out.mllkd(out.j) = mllkd_curr;
            out.mh(out.j)    = alpha;
            out.sigma(out.j) = sigma;
        end
    else
        out.samples(:,i) = v_curr;
        out.mlpt(i)  = mlpt_curr;
        out.mllkd(i) = mllkd_curr;
        out.mh(i)    = alpha;
        out.sigma(i) = sigma;
    end
end

out.acc_rate = acc_t/options.nstep;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stat = init_cov(proposal, stat, sigma, init)

if ~isfield(stat,'M') || ~isfield(stat,'C') || ~isfield(stat,'n')
    stat.n     = 1;
    stat.sum   = init(:);
    stat.cross = init(:)*init(:)';
else
    stat.sum   = stat.M*stat.n;
    stat.cross = stat.C*stat.n + stat.sum(:)*stat.sum(:)'/stat.n;
end
stat = update_cov(proposal, stat, sigma);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stat = update_cov(proposal, stat, sigma)
% Compute the operators for proposals duirng the LIS update
% Tiangang Cui, 25/Mar/2013
%

np     = size(stat.cross, 1);
stat.M = stat.sum/stat.n;
stat.C = stat.cross/stat.n - stat.M(:)*stat.M(:)' + eye(np)*1E-3;

switch proposal
    case {'MALA', 'RW'}
        stat.scale = exp(sigma)/sqrt(max(diag(stat.C)));
        stat.L     = chol(stat.C)';
    case {'OW_Prior'}
        [V,D]   = eig(stat.C);
        [d,ind] = sort(abs(diag(D)), 'descend');
        V       = V(:,ind);
        
        dt      = exp(sigma)/sqrt(sum(d));
        tmp     = ( (0.5*dt)*d + 1 ).^(-1); % eigen values of the B operator
        db      = sqrt(2*dt)*sqrt(d).*tmp;
        da      = (1 - 0.5*dt*d).*tmp;
        
        stat.B  = V*(db(:).*V');
        stat.A  = V*(da(:).*V');
end

end


