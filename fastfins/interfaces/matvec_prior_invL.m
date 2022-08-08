function v = matvec_prior_invL(prior, u)
%MATVEC_PRIOR_L
% 
% L*v
%
% v ~ N(0, I)
% u ~ N(0, C)
%
% Tiangang Cui, 17/Jan/2014

switch prior.type
    case {'Field'}
        v = prior_cov_il(prior, u);
    case {'Basis'}
        v = prior.basis_w'*u;
end

end
