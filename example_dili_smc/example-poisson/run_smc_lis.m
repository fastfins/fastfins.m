
%% options for set the SMC

opts.max_iter = 50; % max number of stages (temperatures)
opts.min_temp = 1E-3;
opts.nsamples = 1E3;  % number of samples used in SMC
opts.ess_tol  = 0.5;  % tol on ESS for finding the next temperature
opts.beta_fac = 1+1E-3;
opts.lis_sq   = false;
opts.tru_tol  = 0.1; % truncation based on error bound of DH, not DH^2
opts.min_rank = 5;
opts.lis_err  = 'Hell';
opts.mlpost   = @(v) minus_log_post(model, obs, prior, v);   % the log-likelihood function
opts.np       = prior.dof;
%
opts.mc_iter  = 5; % number of MCMC iteration per temperature
opts.mc_prop  = 'MALA'; % use a MALA proposal, preconditioned by posterior covairance
opts.mc_sigma = -0.2;
opts.mc_rate  = 0.57; % taget acceptance rate
opts.adapt    = true;
%opts.mc_prop  = 'OW';
%opts.mc_rate  = 0.23;
opts.mc_sigma = -1; % scaling of the proposal
%
opts.mc_npm   = 5; % number of samples for marginlisation

%% run SMC, this is very slow on a serial machine, only for demonstration
[data_sets, betas, log_zs, liss, kernels]  = smc_lis_pm(opts);

%%
% plot the SMC states for the last temperature
figure
plot(data_sets{end,1}(1,:))
title('the first element of full space SMC samples')


figure
plot(data_sets{end,2})
title('minus log likelihood of SMC samples')

