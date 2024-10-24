function [data_sets, betas, log_zs, liss, kernels] = smc_lis_pm(opts)

% initial samples set
betas   = zeros(opts.max_iter, 1);
liss    = cell(opts.max_iter, 1);
kernels = cell(opts.max_iter, 1);
log_zs  = zeros(opts.max_iter, 1);
%
% 1st col: samples, 2nd col: mllkds, 3rd col: grad
data_sets = cell(opts.max_iter, 3);

% initial iteration, use i.i.d. prior
samples = randn(opts.np, opts.nsamples);
mllkds  = zeros(1, opts.nsamples);
gmllkds = zeros(opts.np, opts.nsamples);
for i = 1:opts.nsamples
    [~,mllkds(i),~,gmllkds(:,i)] = opts.mlpost(samples(:,i));
end

beta = 0;
sigma = opts.mc_sigma;

for iter = 1:opts.max_iter
    beta_p = beta;
    %
    disp(['iteration ' num2str(iter)])
    %
    % adapt on the temperature, increase temperature to 1
    beta = max(beta_p,opts.min_temp);
    % compute ess over sample size
    ess = ess_ratio(beta_p, beta, mllkds);
    while ess > opts.ess_tol
        beta = beta*opts.beta_fac;
        ess = ess_ratio(beta_p, beta, mllkds);
    end
    beta = min(1, beta);
    ess = ess_ratio(beta_p, beta, mllkds);
    
    betas(iter) = beta;
    % log of normalising const between two adjacent layers
    log_weight = (beta_p - beta)*mllkds;
    ref_weight = max(log_weight);
    log_zs(iter) = log(mean(exp(log_weight-ref_weight))) + ref_weight;
    %
    disp(['    temperature ' num2str(beta) ', ess ' num2str(ess) ', log const ' num2str(log_zs(iter))])
    
    % compute subspace for the log-likelihood
    weights = exp(log_weight-ref_weight)/sum(exp(log_weight-ref_weight));
    [V,S,~] = svd( gmllkds.*sqrt(weights(:)')*beta, 'econ' );
    s  = diag(S);
    cs = cumsum(s.^2);
    %
    err = 0.5*sqrt(cs(end)-cs);
    r = max(sum(err>opts.tru_tol), opts.min_rank);
    liss{iter}.V = V(:,1:r);
    liss{iter}.s = s;
    liss{iter}.r = r;
    
    disp(['    LIS dimension ' num2str(r)])
    
    % resampling
    ind     = datasample(1:opts.nsamples, opts.nsamples, 'weights', weights);
    samples = samples(:,ind);
    mllkds  = mllkds(ind);
    gmllkds = gmllkds(:,ind);
    v_sub_curr = liss{iter}.V'*samples;
    
    % build MCMC kernel
    kernel  = build_kernel(opts.mc_prop, v_sub_curr, sigma);
    
    switch opts.mc_prop
        case {'MALA'}
            [mllkd_pm_curr, ~, ~, gmllkd_pm_curr, ~] = ...
                minus_sub_pm(opts.mlpost, liss{iter}.V, beta, v_sub_curr, opts.mc_npm);
        otherwise
            mllkd_pm_curr = ...
                minus_sub_pm(opts.mlpost, liss{iter}.V, beta, v_sub_curr, opts.mc_npm);
    end
    
    % iterate on all samples using the MCMC kernel
    for k = 1:opts.mc_iter
        switch opts.mc_prop
            case {'MALA'}
                gmlpt_curr  = gmllkd_pm_curr + v_sub_curr;
                drift_curr  = -(kernel.C*gmlpt_curr)*(kernel.scale^2/2);
                r_curr      = randn(size(v_sub_curr));
                v_sub_next  = v_sub_curr + drift_curr + kernel.scale*(kernel.L*r_curr);
                
                [mllkd_pm_next, mllkd_next, v_next, gmllkd_pm_next, gmllkd_next] = ...
                    minus_sub_pm(opts.mlpost, liss{iter}.V, beta, v_sub_next, opts.mc_npm);
                
                gmlpt_next  = gmllkd_pm_next + v_sub_next;
                drift_next  = -(kernel.C*gmlpt_next)*(kernel.scale^2/2);
                r_next      = -( (kernel.L\(drift_curr + drift_next))/kernel.scale + r_curr );
                log_n2c     = - 0.5 * sum(r_next.^2, 1);
                log_c2n     = - 0.5 * sum(r_curr.^2, 1);
                % row vector
                mlpt_curr   = mllkd_pm_curr + 0.5*sum(v_sub_curr.^2, 1);
                mlpt_next   = mllkd_pm_next + 0.5*sum(v_sub_next.^2, 1);
                
                alpha       = (mlpt_curr - mlpt_next) + (log_n2c - log_c2n);
            case {'RW'}
                v_sub_next  = v_sub_curr + kernel.scale*(kernel.L*randn(size(v_sub_curr)));
                [mllkd_pm_next, mllkd_next, v_next] = ...
                    minus_sub_pm(opts.mlpost, liss{iter}.V, beta, v_sub_next, opts.mc_npm);
                %
                mlpt_curr   = mllkd_pm_curr + 0.5*sum(v_sub_curr.^2, 1);
                mlpt_next   = mllkd_pm_next + 0.5*sum(v_sub_next.^2, 1);
                alpha       = mlpt_curr - mlpt_next;
            case {'OW'}
                v_sub_next  = kernel.A*v_sub_curr + kernel.B*randn(size(v_sub_curr));
                [mllkd_pm_next, mllkd_next, v_next] = ...
                    minus_sub_pm(opts.mlpost, liss{iter}.V, beta, v_sub_next, opts.mc_npm);
                alpha       = mllkd_pm_curr - mllkd_pm_next;
        end
        
        acc = log(rand(1, opts.nsamples)) < alpha;
        
        disp(['    MC iter ' num2str(k) ...
            ', acceptance rate ' num2str(sum(acc)/opts.nsamples) ...
            ', sigma ' num2str(sigma)])
        % update
        v_sub_curr(:,acc) = v_sub_next(:,acc);
        mllkd_pm_curr(acc) = mllkd_pm_next(acc);
        %
        samples(:,acc) = v_next(:,acc);
        mllkds(acc) = mllkd_next(acc);
        %
        switch opts.mc_prop
            case {'MALA'}
                gmllkd_pm_curr(:,acc) = gmllkd_pm_next(:,acc);
                gmllkds(:,acc) = gmllkd_next(:,acc);
        end
        %
        if opts.adapt
            delta   = min(0.05,1/sqrt(opts.mc_iter*iter+k));
            if sum(acc)/opts.nsamples < opts.mc_rate
                sigma = sigma - delta;
            else
                sigma = sigma + delta;
            end
            kernel  = build_kernel(opts.mc_prop, v_sub_curr, sigma);
        end
    end
    % end of MCMC
    switch opts.mc_prop
        case {'OW', 'RW'}
            for i = 1:opts.nsamples
                [~,~,~,gmllkds(:,i)] = opts.mlpost(samples(:,i));
            end
    end
    
    % end of iteration
    data_sets{iter,1} = samples;
    data_sets{iter,2} = mllkds;
    data_sets{iter,3} = gmllkds;
    kernels{iter} = kernel;
    %
    [V,S,~] = svd( gmllkds/sqrt(opts.nsamples)*beta, 'econ' );
    s = diag(S);
    liss{iter}.post_V = V(:,1:100);
    liss{iter}.post_s = s;
    
    if abs(beta - 1) < 1E-10
        break;
    end
end

ind = find(abs(betas-1)<1E-5);
betas = betas(1:ind);
liss = liss(1:ind);
kernels = kernels(1:ind);
log_zs = log_zs(1:ind);
data_sets = data_sets(1:ind,:);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kernel = build_kernel(proposal, samples, sigma)
% Compute the operators for proposals duirng the LIS update
% Tiangang Cui, 03/August/2020
%

kernel.C = cov(samples') + eye(size(samples,1))*1E-5;

switch proposal
    case {'MALA', 'RW'}
        kernel.scale = exp(sigma)/sqrt(max(diag(kernel.C)));
        kernel.L = chol(kernel.C)';
    case {'OW'}
        [V,D] = eig(kernel.C);
        [d,ind] = sort(abs(diag(D)), 'descend');
        V   = V(:,ind);
        
        %dt      = exp(sigma)/sqrt(sum(d));
        dt = exp(sigma);
        %disp([dt, sum(d)])
        tmp = ( (0.5*dt)*d + 1 ).^(-1); % eigen values of the B operator
        DB  = sqrt(2*dt)*sqrt(d).*tmp;
        DA  = (1 - 0.5*dt*d).*tmp;
        
        kernel.B = (V.*DB(:)')*V';
        kernel.A = (V.*DA(:)')*V';
end

end