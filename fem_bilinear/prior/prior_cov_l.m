function u = prior_cov_l(prior, v)
%COV_LV
%
% v ~ N(0, I)
% u ~ N(0, C)
%
% Tiangang Cui, 17/Jan/2014

switch prior.cov_type
    case {'MRF'}
        u   = zeros(size(v));
        tmp = v(prior.per,:).*prior.sca(:);
        u(prior.per,:) = prior.R\(prior.R'\tmp);
    case {'GP'}
        u   = prior.RC'*v;
end

end
