function out = rto_centered(model, obs, prior, hyper_setup, rto_setup, nstep, a_init, b_init)
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
rto_lin     = rto_local_lin(model, curr_h.obs, curr_h.prior, rto_setup, HI, vmap);
curr_p      = rto_one_sample(model, curr_h.obs, curr_h.prior, rto_lin);

i = 1;
out.alphas(:,i) = curr_h.alpha;
out.betas(:,i)  = curr_h.beta;
out.lps_h(i)    = curr_h.lp;
out.us(:,i)     = curr_p.u;
out.llkds(i)    = curr_p.llkd;
out.lps(i)      = curr_p.lp;
out.acc         = 0;

adapt_H         = hyper_proposal_init(hyper_setup);
if ~adapt_H.isgibbs
    out.sigma   = zeros(floor(nstep/adapt_H.nbatch),1);
end

k = 0;
l = 0;
batch = 0;
for i = 2:(nstep+1)
    % hyper parameters
    
    if adapt_H.isgibbs
        % sample alpha and beta
        %{
        tmp1 = (curr_p.d - obs.data)./obs.std;
        tmp2 = matvec_prior_invL (prior, curr_p.u - prior.mean_u);
        norm(full(prior.RQ(:) - curr_h.prior.RQ(:)/sqrt(curr_h.beta)))
        disp([0.5*norm(tmp1)^2, -curr_p.llkd/curr_h.alpha])
        disp([0.5*norm(tmp2)^2, -curr_p.lp/curr_h.beta])
        %}
        alpha       = gamrnd(hyper_setup.a0 + obs.Ndata*0.5, 1/(hyper_setup.t0 - curr_p.llkd/curr_h.alpha) );
        beta        = gamrnd(hyper_setup.a0 + prior.DoF*0.5, 1/(hyper_setup.t0 - curr_p.lp/curr_h.beta) );
        curr_h      = update_model(obs, prior, hyper_setup, alpha, beta);
        [vmap,HI]   = rto_cond_map(model, curr_h.obs, curr_h.prior, zeros(prior.DoF, 1));
        rto_lin     = rto_local_lin(model, curr_h.obs, curr_h.prior, rto_setup, HI, vmap);
    else
        [next_h,lr,inc] = hyper_propose(obs, prior, hyper_setup, adapt_H, curr_h);
        cl_cond_post    = log_cond_post(curr_h.obs, curr_h.prior, curr_p.d, curr_p.u);
        nl_cond_post    = log_cond_post(next_h.obs, next_h.prior, curr_p.d, curr_p.u);
        
        if  log(rand) < (nl_cond_post + next_h.lp - cl_cond_post - curr_h.lp + lr)
            curr_h      = next_h;
            [vmap,HI]   = rto_cond_map(model, curr_h.obs, curr_h.prior, zeros(prior.DoF, 1));
            rto_lin     = rto_local_lin(model, curr_h.obs, curr_h.prior, rto_setup, HI, vmap);
            adapt_H.acc = adapt_H.acc + inc;
            disp('ACC_H')
        end
        batch = batch + inc;
    end
    
    disp([curr_h.alpha, curr_h.beta])
    
    curr_p  = rto_update_density(model, curr_h.obs, curr_h.prior, rto_lin, curr_p);
    % RTO
    for j = 1:rto_setup.substep
        next_p      = rto_one_sample(model, curr_h.obs, curr_h.prior, rto_lin);
        
        %disp([next_p.logrho curr_p.logrho])
        
        if log(rand) < (next_p.logrho - curr_p.logrho)
            curr_p  = next_p;
            disp('ACC_P')
            l = l+1;
        end
    end
    
    if ~adapt_H.isgibbs && batch == adapt_H.nbatch
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

out.l           = l;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = log_cond_post(obs, prior, md, u)

misfit  = (md - obs.data)./obs.std;
mllkd   = 0.5*sum(misfit(:).^2);
v       = matvec_prior_invL (prior, u - prior.mean_u);
mlp     = 0.5*sum(v(:).^2);

f       = - mllkd - mlp;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
