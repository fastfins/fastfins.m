function adapt = hyper_proposal_adapt(adapt, i, samples)

tmp     = log([adapt.init_samples, samples]);
% factorize new covariance
C       = cov(tmp');
C       = max(diag(C));

%{
delta   = min(0.1,sqrt(adapt.nbatch/i));
if (adapt.acc/adapt.nbatch) < adapt.rate
    adapt.sigma = adapt.sigma - delta;
else
    adapt.sigma = adapt.sigma + delta;
end

%cond(C)
scale   = exp(adapt.sigma)/sqrt(max(diag(C)));
adapt.L = chol(C)'*scale;
%}
adapt.acc   = 0;


disp(adapt.sigma)
adapt.L = chol(C)'*exp(adapt.sigma);


if adapt.AM_weight < 1 && mod(i, adapt.GM_batch) == 0
    adapt.GM    = fitgmdist(tmp', adapt.GM_ngauss, 'Options', adapt.GM_options);
end

end