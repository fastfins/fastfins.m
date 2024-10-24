function ess2n = ess_ratio(beta_p, beta, mllkds)

log_weight = (beta_p - beta)*mllkds;
log_weight = log_weight - max(log_weight);
ess2n = ( sum(exp(log_weight)).^2/sum(exp(2*log_weight)) ) / length(mllkds(:));

end