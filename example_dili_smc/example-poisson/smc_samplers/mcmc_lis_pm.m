function [data, mllkds] = mcmc_lis_pm(opts, lis, vinit, kernel_init)

beta = 1;
v_curr = vinit;
sigma = opts.mc_sigma;

nsamples = size(v_curr, 2);
data = zeros(size(v_curr,1), nsamples*opts.mc_iter);
mllkds = zeros(1, nsamples*opts.mc_iter);

% initial iteration
v_sub_curr = lis.V'*v_curr;
switch opts.mc_prop
    case {'MALA'}
        [mllkd_pm_curr, mllkd_curr, v_curr, gmllkd_pm_curr] = ...
            minus_sub_pm(opts.mlpost, lis.V, beta, v_sub_curr, opts.mc_npm);
    otherwise
        mllkd_pm_curr = ...
            minus_sub_pm(opts.mlpost, lis.V, beta, v_sub_curr, opts.mc_npm);
end

% build MCMC kernel
kernel = build_kernel(opts.mc_prop, kernel_init, sigma);

acc_num = 0;
batch = 0;

for iter = 1:opts.mc_iter

    switch opts.mc_prop
        case {'MALA'}
            gmlpt_curr  = gmllkd_pm_curr + v_sub_curr;
            drift_curr  = -(kernel.C*gmlpt_curr)*(kernel.scale^2/2);
            r_curr      = randn(size(v_sub_curr));
            v_sub_next  = v_sub_curr + drift_curr + kernel.scale*(kernel.L*r_curr);

            [mllkd_pm_next, mllkd_next, v_next, gmllkd_pm_next] = ...
                minus_sub_pm(opts.mlpost, lis.V, beta, v_sub_next, opts.mc_npm);

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
                minus_sub_pm(opts.mlpost, lis.V, beta, v_sub_next, opts.mc_npm);
            %
            mlpt_curr   = mllkd_pm_curr + 0.5*sum(v_sub_curr.^2, 1);
            mlpt_next   = mllkd_pm_next + 0.5*sum(v_sub_next.^2, 1);
            alpha       = mlpt_curr - mlpt_next;
        case {'OW'}
            v_sub_next  = kernel.A*v_sub_curr + kernel.B*randn(size(v_sub_curr));
            [mllkd_pm_next, mllkd_next, v_next] = ...
                minus_sub_pm(opts.mlpost, lis.V, beta, v_sub_next, opts.mc_npm);
            alpha       = mllkd_pm_curr - mllkd_pm_next;
    end

    acc = log(rand(1, nsamples)) < alpha;

    % update
    v_sub_curr(:,acc) = v_sub_next(:,acc);
    mllkd_pm_curr(acc) = mllkd_pm_next(acc);
    %
    v_curr(:,acc) = v_next(:,acc);
    mllkd_curr(acc) = mllkd_next(acc);
    %
    ind = (1:nsamples) + (iter-1)*nsamples;
    data(:,ind) = v_curr;
    mllkds(ind) = mllkd_curr;
    %
    acc_num = acc_num + sum(acc);
    %
    switch opts.mc_prop
        case {'MALA'}
            gmllkd_pm_curr(:,acc) = gmllkd_pm_next(:,acc);
    end
    %

    if opts.adapt
        batch = batch + 1;
        % update covariance
        kernel.cross = kernel.cross + v_sub_curr*v_sub_curr';
        kernel.sum = kernel.sum + sum(v_sub_curr,2);
        kernel.num = kernel.num + nsamples;
        
        if  batch == opts.nbatch
            delta = min(0.05,sqrt(opts.nbatch/iter));
            if (acc_num/(opts.nbatch*nsamples)) < opts.mc_rate
                sigma = sigma - delta;
            else
                sigma = sigma + delta;
            end
            kernel = build_kernel(opts.mc_prop, kernel, sigma);

            disp(['    MC iter ' num2str(iter) ...
                ', acceptance rate ' num2str(sum(acc_num)/(nsamples*opts.nbatch)) ...
                ', sigma ' num2str(sigma)])
            batch = 0;
            acc_num = 0;
            
        end
    end

end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kernel = build_kernel(proposal, kernel, sigma)
% Compute the operators for proposals duirng the LIS update
% Tiangang Cui, 25/Mar/2013
%

M = kernel.sum/kernel.num;
kernel.C = kernel.cross/kernel.num - M(:)*M(:)' + eye(length(kernel.sum))*1E-3;

switch proposal
    case {'MALA', 'RW'}

        kernel.scale    = exp(sigma)/sqrt(max(diag(kernel.C)));
        kernel.L        = chol(kernel.C)';

    case {'OW'}

        [V,D]   = eig(kernel.C);
        [d,ind] = sort(abs(diag(D)), 'descend');
        V       = V(:,ind);

        %dt      = exp(sigma)/sqrt(sum(d));
        dt      = exp(sigma);
        %disp([dt, sum(d)])
        tmp     = ( (0.5*dt)*d + 1 ).^(-1); % eigen values of the B operator
        DB      = sqrt(2*dt)*sqrt(d).*tmp;
        DA      = (1 - 0.5*dt*d).*tmp;

        kernel.B    = scale_cols(V, DB)*V';
        kernel.A    = scale_cols(V, DA)*V';

end

end


