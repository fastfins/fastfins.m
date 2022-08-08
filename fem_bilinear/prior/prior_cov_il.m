function v = prior_cov_il(prior, u)
%COV_INVLU 
%
% whitening transformation
%
% Tiangang Cui, 17/Jan/2014

switch prior.cov_type
    case {'MRF'}
        v = prior.RQ *u;
    case {'GP'}
        v = prior.RC'\u;
end

end