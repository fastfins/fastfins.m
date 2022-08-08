function [LISs, outs] = rto_lis(model, prior, obs, rto_setup)
% building LIS using RTO, and run RTO within LIS
%
% Tiangang Cui, Nov 11, 2015
%{
% need to define the following
rto_setup.neig
rto_setup.lis_thres
rto_setup.debug
rto_setup.Nlis
rto_setup.Nstep

rto_setup.large_scale_flag
rto_setup.svd_tol
rto_setup.nmax
rto_setup.max_iter
%}

vmap            = get_map_matlab(model, obs, prior, randn(prior.DoF, 1)); 

LISs            = cell(rto_setup.max_iter+1, 1);
outs            = cell(rto_setup.max_iter+1, 1);

% Build LIS using Laplace at MAP
LISs{1}         = build_lis_laplace(model, prior, obs, vmap, rto_setup);

for i = 2:(rto_setup.max_iter+1)
    outs{i-1}   = rto_mcmc(model, obs, LISs{i-1}, rto_setup.Nlis*2, rto_setup);
    ind         = randperm(rto_setup.Nlis*2, rto_setup.Nlis);
    sub_samples = outs{i-1}.v_samples(ind,:)';
    lpt_jr      = outs{i-1}.lpt(ind);
    
    LISs{i}     = build_lis_IS(model, prior, obs, ...
        LISs{i-1}, lpt_jr, sub_samples, rto_setup);
end

outs{rto_setup.max_iter+1} ...
                = rto_mcmc(model, obs, LISs{romset.max_iter+1}, rto_setup.Nstep, rto_setup);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LIS = build_lis_laplace(model, prior, obs, vmap, rto_setup)

n           = rto_setup.Nlis;

[V, d]      = eigen_PPGNH(model, obs, prior, vmap, 1E-2, rto_setup.neig);
s           = (d+1).^(-0.5);
r           = randn(prior.DoF, n);
css         = r - V*(V'*r);
subs        = V*scale_rows(V'*r, s);
samples     = css + subs + repmat(vmap, 1, n);

crossV      = zeros(prior.DoF);
tic;
for i = 1:n
    [V, d]  = eigen_PPGNH(model, obs, prior, samples(:,i), 1E-2, rto_setup.neig);
    crossV  = crossV + V*scale_cols(V, d)';
end

cross       = crossV/n;

% dim redu
[pU,pD]     = eig(cross);
[pd,ind]    = sort(real(diag(pD)), 'descend');
pU          = pU(:,ind);
jnd         = pd>=rto_setup.lis_thres;
LIS         = basis_LIS(prior, pU(:,jnd)); % new LIS

LIS.ess_c   = n;
LIS.ess_t   = 0;

LIS.cross   = cross;
    
if rto_setup.debug
    LIS.crossV  = crossV;
end

LIS.DoF

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LIS = build_lis_IS(model, prior, obs, LISp, lpt_jr, sub_samples, rto_setup)

% romset.deim_factor
% romset.lognorm_flag
% romset.lis_thres
% romset.Nlis

n               = size(sub_samples, 2);
r               = randn(prior.DoF, n);
nul_samples     = r - LISp.P*(LISp.P'*r);
samples         = LISp.P*sub_samples + nul_samples;

lpr_nu          = 0.5*sum(nul_samples.^2, 1)';

lpt             = zeros(n, 1);
crossV          = zeros(prior.DoF);

tic;
for i = 1:n
    lpt(i)      = minus_log_post(model, obs, prior, samples(:,i));
end

tmp     = lpt_jr + lpr_nu - lpt;
tmp1    = sort(tmp, 'descend');
tmp2    = tmp - mean(tmp1(1:10));

weights         = exp(tmp2);
 
%{
LIS.lpt         = lpt;
LIS.lpr_nu        = lpr_nu2;
LIS.lpt_jr        = lpt_jr2;
LIS.weights   = weights;
return
%}

%

for i = 1:n
    [V, d]      = eigen_PPGNH(model, obs, prior, samples(:,i), 1E-2, rto_setup.neig);
    crossV      = crossV + V*scale_cols(V, d)'*weights(i);
end

ess_c           = n*(mean(weights))^2/mean(weights.^2);
cross           = ( ( crossV/sum(weights) )*ess_c + LISp.cross*LISp.ess_t ) / (ess_c + LISp.ess_t);

% dim redu
[pU,pD]         = eig(cross);
[pd,ind]        = sort(real(diag(pD)), 'descend');
pU              = pU(:,ind);
jnd             = pd>=rto_setup.lis_thres;
LIS             = basis_LIS(prior, pU(:,jnd)); % new LIS

LIS.ess_c       = ess_c;
LIS.ess_t       = ess_c + LISp.ess_t;
LIS.weights     = weights;
LIS.cross       = cross;

if rto_setup.debug
    LIS.crossV  = crossV;
end

LIS.DoF

end

