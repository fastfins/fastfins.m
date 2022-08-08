function adapt = hyper_proposal_init(hyper_setup)

if ~isfield(hyper_setup, 'proposal_type')
    adapt.proposal_type = 'AM';
    disp('ERR: no proposal type!')
else
    adapt.proposal_type = hyper_setup.proposal_type;
end

adapt.isgibbs = strcmp(adapt.proposal_type, 'gibbs');

if adapt.isgibbs
    return;
end

if ~isfield(hyper_setup, 'nbatch')
    adapt.nbatch    = 50;
else
    adapt.nbatch    = hyper_setup.nbatch;
end
% default acceptance rate
if ~isfield(hyper_setup, 'rate')
    adapt.rate      = 0.23;
else
    adapt.rate      = hyper_setup.rate;
end
% default jump size
if ~isfield(hyper_setup, 'sigma')
    adapt.sigma     = -1;
else
    adapt.sigma     = hyper_setup.sigma;
end
%{
if ~isfield(hyper_setup, 'init_weight')
    adapt.init_weight = hyper_setup.n_hyper;
else
    adapt.init_weight = hyper_setup.init_weight;
end
%}
%
if ~isfield(hyper_setup, 'init_samples')
    adapt.init_samples  = exp(randn(hyper_setup.n_hyper));
    disp('IID proposal')
else
    adapt.init_samples  = hyper_setup.init_samples;
end
%{
if ~isfield(hyper_setup, 'init_mean')
    adapt.init_mean = mean(adapt.init_samples, 2);
else
    adapt.init_mean = hyper_setup.init_mean;
end
%
if ~isfield(hyper_setup, 'init_2nd')
    adapt.init_2nd  = cov(adapt.init_samples') + adapt.init_mean*adapt.init_mean';
else
    adapt.init_2nd  = hyper_setup.init_2nd;
end

% factorize new covariance
C       = adapt.init_2nd - adapt.init_mean(:)*adapt.init_mean(:)';
%}

%{
C           = cov(log(adapt.init_samples)');
scale       = exp(adapt.sigma)/sqrt(max(diag(C)));
adapt.L     = chol(C)'*scale;
adapt.acc   = 0;
%}

C           = cov(log(adapt.init_samples)');
adapt.L     = chol(C)'*exp(adapt.sigma);
adapt.acc   = 0;

if strcmp(hyper_setup.proposal_type, 'GM')
    if ~isfield(hyper_setup, 'GM_ngauss')
        adapt.GM_ngauss = 5;
    else
        adapt.GM_ngauss = hyper_setup.GM_ngauss;
    end
    %
    if ~isfield(hyper_setup, 'GM_weight')
        adapt.AM_weight = 0.35;
    else
        adapt.AM_weight = 1 - hyper_setup.GM_weight;
    end
    %
    if ~isfield(hyper_setup, 'GM_batch')
        adapt.GM_batch  = 1E3;
    else
        adapt.GM_batch  = hyper_setup.GM_batch;
    end
    %
    adapt.GM_options    = statset('MaxIter',1000);
    adapt.GM            = fitgmdist(log(adapt.init_samples)', adapt.GM_ngauss, 'Options', adapt.GM_options);
else
    adapt.AM_weight     = 1;
end

end