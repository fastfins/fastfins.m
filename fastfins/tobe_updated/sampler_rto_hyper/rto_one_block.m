function out = rto_one_block(model, obs, prior, hyper_setup, rto_setup, nstep, a_init, b_init)
%RTO_MCMC
%
% RTO sampler
% Tiangang Cui, August, 1, 2017

rto_setup   = rto_init(rto_setup);

% rto samples
out.alphas  = zeros(size(a_init(:),2), nstep+1);
out.betas   = zeros(size(b_init(:),2), nstep+1);
out.lps_h   = zeros(1, nstep+1);
out.us      = zeros(prior.DoF, nstep+1);
out.llkds   = zeros(1, nstep+1);
out.lps     = zeros(1, nstep+1);

curr_h      = update_model(obs, prior, hyper_setup, a_init, b_init);
% find map and initialise RTO
[vmap,HI]   = rto_cond_map(model, curr_h.obs, curr_h.prior, zeros(prior.DoF, 1));
curr_h.rto_lin  = rto_local_lin(model, curr_h.obs, curr_h.prior, rto_setup, HI, vmap);
curr_p      = rto_one_sample(model, curr_h.obs, curr_h.prior, curr_h.rto_lin);

i = 1;
out.alphas(:,i) = curr_h.alpha;
out.betas(:,i)  = curr_h.beta;
out.lps_h(i)    = curr_h.lp;
out.us(:,i)     = curr_p.u;
out.llkds(i)    = curr_p.llkd;
out.lps(i)      = curr_p.lp;

adapt_H         = hyper_proposal_init(hyper_setup);

k = 0;
la = 0;
lb = 0;
batch = 0;

for i = 2:(nstep+1)
    % given lambda, delta, v and a current linearization
    % propose
    [next_h,lr,inc] = hyper_propose(obs, prior, hyper_setup, adapt_H, curr_h);
    
    [vmap,HI]       = rto_cond_map  (model, next_h.obs, next_h.prior, zeros(prior.DoF, 1));
    next_h.rto_lin  = rto_local_lin (model, next_h.obs, next_h.prior, rto_setup, HI, vmap);
    next_p          = rto_one_sample(model, next_h.obs, next_h.prior, next_h.rto_lin);
    
    %disp([next_p.logrho next_h.lp curr_p.logrho curr_h.lp])
    %disp(next_p.logrho + next_h.lp - curr_p.logrho - curr_h.lp)
    
    disp([curr_h.alpha, curr_h.beta, next_h.alpha, next_h.beta])
    disp([curr_p.logrho + curr_h.lp, next_p.logrho + next_h.lp])
    %disp(lr)
    if log(rand) < (next_p.logrho + next_h.lp - curr_p.logrho - curr_h.lp + lr)
        curr_h      = next_h;
        curr_p      = next_p;
        adapt_H.acc = adapt_H.acc + inc;
        disp('>>>> ACC')
        
        la = la + 1;
    end

    batch = batch + inc;
    
    if batch == adapt_H.nbatch
        samples = [out.alphas(:,1:i-1); out.betas(:,1:i-1)];
        adapt_H   = hyper_proposal_adapt(adapt_H, i, samples);
        
        batch           = 0;
        k               = k+1;
        out.sigma(k,:)  = adapt_H.sigma;
    end
    
    out.alphas(:,i) = curr_h.alpha;
    out.betas(:,i)  = curr_h.beta;
    out.lps_h(i)    = curr_h.lp;
    out.us(:,i)     = curr_p.u;
    out.llkds(i)    = curr_p.llkd;
    out.lps(i)      = curr_p.lp;
    
    disp(i)
end

out.la          = la;
out.lb          = lb;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
