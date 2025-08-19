
%%
opts.mlpost   = @(v) minus_log_post(model, obs, prior, v);   % the log-likelihood function
opts.np       = prior.dof;
%
opts.mc_iter  = 5E3;
opts.mc_prop  = 'MALA';
opts.mc_sigma = -0.5;
opts.mc_rate  = 0.57;
opts.adapt    = true;
%opts.mc_prop  = 'OW';
%opts.mc_rate  = 0.23;
%opts.mc_sigma = -1;
opts.nbatch = 50; % adaptation step interval
%
opts.mc_npm   = 5; % number of samples for marginlisation

%%
% either build LIS using example_dili, or using the data free one
n_MC = 200;
prior_samples = randn(prior.dof, n_MC);
H = zeros(opts.np);
for i = 1:n_MC
    [~,~,~,~,sol] = minus_log_post(model, obs, prior, prior_samples(:,i));
    Ju = explicit_jacobian(model, sol);
    Jv = matvec_prior_Lt(prior, Ju')';
    H  = H + Jv'*sol.I*Jv;
end
H = H / n_MC;
[V, S, ~] = svd(H);
s = diag(S).^0.5;

% select s > 0.1
lis_df.r = sum(s>0.1);
lis_df.s = s(1:lis_df.r);
lis_df.V = V(:,1:lis_df.r);

%%
% build an initial kernel for MCMC, sample sum and sum of outter product

lis = lis_dili;

m = 3; % number of chains
v_sub = lis.V'*vmap;
kernel_init.num = 100;
kernel_init.sum = (v_sub)*kernel_init.num;
kernel_init.cross = (diag(1./(lis.s.^2+1)) + v_sub*v_sub')*kernel_init.num;

% running three chains in parallel
[data, mllkds]  = mcmc_lis_pm(opts, lis, vmap*ones(1,m), kernel_init);

%%
% plot the first element of the first chain
figure
plot(data(1,1:m:end))
title('the first element of full space MCMC samples, 1st chain out of 3')


figure
plot(mllkds(1:m:end))
title('minus log likelihhod of full space MCMC samples, 1st chain out of 3')