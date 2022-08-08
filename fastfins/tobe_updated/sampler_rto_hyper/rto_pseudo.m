function out = rto_pseudo(model, obs, prior, hyper_setup, rto_setup, nstep, nmar, a_init, b_init)
%RTO_MCMC
%
% RTO sampler
% Tiangang Cui, August, 1, 2017

% rto samples
out.alphas  = zeros(size(a_init(:),2), nstep+1);
out.betas   = zeros(size(b_init(:),2), nstep+1);
out.lps_h   = zeros(1, nstep+1);
out.llkds_m = zeros(1, nstep+1);

curr        = update_model(obs, prior, hyper_setup, a_init, b_init);
[vmap,HI]   = rto_cond_map(model, curr.obs, curr.prior, zeros(prior.DoF, 1));
curr.rto_lin    = rto_local_lin(model, curr.obs, curr.prior, rto_setup, HI, vmap);
curr.llkd_m     = rto_marginal_lkd(model, curr.obs, curr.prior, curr.rto_lin, nmar);

out.alphas(:,1) = a_init;
out.betas(:,1)  = b_init;
out.lps_h(1)    = curr.lp;
out.llkds_m(1)  = curr.llkd_m; 


adapt           = hyper_proposal_init(hyper_setup);
if ~adapt.isgibbs
    out.sigma   = zeros(floor(nstep/adapt.nbatch),1);
end

k = 0;
batch = 0;
for i = 2:(nstep+1)
    % given lambda, delta, v and a current linearization
    % propose
    [next,lr,inc]   = hyper_propose(obs, prior, hyper_setup, adapt, curr);
    [vmap,HI]       = rto_cond_map    (model, next.obs, next.prior, zeros(prior.DoF, 1));
    next.rto_lin    = rto_local_lin   (model, next.obs, next.prior, rto_setup, HI, vmap);
    [next.llkd_m,s] = rto_marginal_lkd(model, next.obs, next.prior, next.rto_lin, nmar);
    
    disp(s)
    disp([next.llkd_m + next.lp, curr.llkd_m + curr.lp])
    if log(rand) < (next.llkd_m + next.lp - curr.llkd_m - curr.lp + lr)
        curr = next;
        adapt.acc   = adapt.acc + inc;
        disp('ACC')
    end
    
    batch = batch + inc;
    
    if batch == adapt.nbatch
        samples = [out.alphas(:,1:i-1); out.betas(:,1:i-1)];
        adapt   = hyper_proposal_adapt(adapt, i, samples);
        
        batch           = 0;
        k               = k+1;
        out.sigma(k,:)  = adapt.sigma;
    end
    
    out.alphas(:,i) = curr.alpha;
    out.betas(:,i)  = curr.beta;
    out.lps_h(i)    = curr.lp;
    out.llkds_m(i)  = curr.llkd_m;
    
    disp(i)
end

end
