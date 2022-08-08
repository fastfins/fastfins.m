function u = prior_cov_lt(prior, v)
%COV_LTV
%
% Tiangang Cui, 17/Jan/2014

switch prior.cov_type
    case {'MRF'}
        tmp = prior.R\(prior.R'\v(prior.per,:));
        u(prior.per,:) = tmp.*prior.sca(:);
    case {'GP'}
        u = prior.RC*v;
end

end