function lis = basis_LIS(prior, P)
%BASIS_LIS
%
% P:     new basis
% S:     singular values associated with P
% v_0:   reference variable, without prior mean
%
% Tiangang Cui, 20/Oct/2012

lis.P       = P;
lis.basis   = prior_cov_l  (prior, P); 
lis.basis_w = prior_cov_ilt(prior, P); 
lis.dof     = size(P,2);
lis.mean_u  = prior.mean_u;

lis.type    = 'Basis';
lis.note    = 'LIS';

end

