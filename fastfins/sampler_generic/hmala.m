function out = hmala(ml_target, hessian, init, options)
%HMALA
%
% Hessian preconditioned Metropolis adjusted Langevin MCMC sampler
%
% Tiangang Cui, 17/Jan/2014
%
% Reference:
% Martin, J., Wilcox, L. C., Burstedde, C., & Ghattas, O. (2012). A stocha-
% stic Newton MCMC method for large-scale statistical inverse problems with
% application to seismic inversion. SIAM Journal on Scientific Computing,
% 34(3), A1460-A1487.
%
% Tiangang Cui, 07/July/2019

v_curr  = init;
np      = length(init);
sigma   = options.sigma;
scale   = exp(sigma);
out     = mcmc_outputs(np, options);

[mlpt_curr,mllkd_curr,mg_curr,~] = ml_target(v_curr);
la_curr = get_laplace(hessian, v_curr);
la_next = la_curr;
Hg_curr = mg_curr + la_curr.V_iH*(la_curr.V'*mg_curr);

if  options.lis_flag
    t1          = mg_curr - v_curr;
    out.cross_g = t1*t1';
    out.cross_v = v_curr*v_curr';
    out.sum_v   = v_curr;
    out.nsample = 1;
end

% start MCMC
acc_t = 0;
acc   = 0;
batch = 0;

for i = 1:options.nstep
    r_curr  = randn(np,1);
    Hr_curr = r_curr + la_curr.V_ihalf*(la_curr.V'*r_curr);
    v_next  = v_curr - (0.5*scale^2)*Hg_curr + scale*Hr_curr;
    
    if  strcmp(options.proposal, 'Newton')
        [mlpt_next,mllkd_next,mg_next,~,HI] = ml_target(v_next);
        la_next = get_laplace(hessian, HI);
    else
        [mlpt_next,mllkd_next,mg_next] = ml_target(v_next);
    end
    
    Hg_next = mg_next + la_next.V_iH*(la_next.V'*mg_next);
    Hr_next = (0.5*scale)*(Hg_curr+Hg_next) - Hr_curr;
    tmp     = la_next.V_half' * Hr_next;
    log_c2n = - 0.5 * (r_curr'*r_curr) + la_curr.ldet_half;
    log_n2c = - 0.5 * ( norm(Hr_next)^2 + norm(tmp)^2 ) + la_next.ldet_half;
    
    alpha   = (mlpt_curr - mlpt_next) + (log_n2c - log_c2n); % log acceptance prob.
    
    if  log(rand) < alpha
        v_curr  = v_next;
        mlpt_curr   = mlpt_next;
        mllkd_curr  = mllkd_next;
        %
        mg_curr = mg_next;
        Hg_curr = Hg_next;
        %
        if strcmp(options.proposal, 'Newton')
            la_curr = la_next;
        end
        acc   = acc+1;
        acc_t = acc_t+1;
    end
    
    if  options.adapt
        batch = batch + 1;
        if  batch == options.nbatch
            delta = min(0.1,sqrt(options.nbatch/i));
            if (acc/options.nbatch) < options.rate
                sigma = sigma - delta;
            else
                sigma = sigma + delta;
            end
            scale = exp(sigma);
            batch = 0;
            acc   = 0;
        end
    end
    
    if  options.lis_flag
        t1          = mg_curr - v_curr;
        out.cross_g = out.cross_g + t1*t1';
        out.cross_v = out.cross_v + v_curr*v_curr';
        out.sum_v   = out.sum_v + v_curr;
        out.nsample = out.nsample+1;
    end
    
    % save    
    if  options.sbatch > 1
        if mod(i, options.sbatch) == 0
            out.j = out.j + 1;
            out.samples(:,out.j) = v_curr;
            out.mlpt(out.j)      = mlpt_curr;
            out.mllkd(out.j)     = mllkd_curr;
            out.mh(out.j)        = alpha;
            out.sigma(out.j)     = sigma;
        end
    else
        out.samples(:,i) = v_curr;
        out.mlpt(i)      = mlpt_curr;
        out.mllkd(i)     = mllkd_curr;
        out.mh(i)        = alpha;
        out.sigma(i)     = sigma;
    end
    
end

out.acc_rate = acc_t/options.nstep;

end

