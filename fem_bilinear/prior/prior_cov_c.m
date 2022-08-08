function u = prior_cov_c(prior, v)
%COV_CV
%
% Tiangang Cui, 17/Jan/2014

switch prior.cov_type
    case {'MRF'}
        tmp1 = prior.R\(prior.R'\v(prior.per,:));
        tmp2(prior.per,:) = tmp1.*prior.sca(:);
        u    = zeros(size(v));
        tmp3 = tmp2(prior.per,:).*prior.sca(:);
        u(prior.per,:) = prior.R\(prior.R'\tmp3);
    case {'GP'}
        u = prior.C*v;
end

end