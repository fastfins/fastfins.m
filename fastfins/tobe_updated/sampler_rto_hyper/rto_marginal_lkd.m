function [llkd_m, stdrho] = rto_marginal_lkd(model, obs, prior, rto_lin, nsample)
%Compute marginal likelihood
%
% Tiangang Cui, August, 1, 2017

out     = rto_samples(model, obs, prior, rto_lin, nsample);

mlw     = max(out.logrhos);
llkd_m  = log( sum(exp(out.logrhos - mlw))/length(out.logrhos) ) + mlw;
stdrho  = std(out.logrhos)/mean(out.logrhos);

end