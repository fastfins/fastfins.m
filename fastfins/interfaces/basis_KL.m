function prior = basis_KL(prior, thres)
%BASIS_KL
% 
% Tiangang Cui, 17/Jan/2014

[V, d]          = prior_cov_eig(prior);

if thres < 1
    jnd         = cumsum(d)/sum(d) <= thres; % truncation 
    prior.dof   = sum(jnd);
else
    jnd         = 1:floor(thres);
    prior.dof   = floor(thres);
end
prior.P         = V(:, jnd);
S               = d(jnd).^(0.5);
prior.basis     = prior.P.*S(:)';
prior.basis_w   = prior.P.*(S(:)'.^(-1));
prior.mean_v    = zeros(prior.dof, 1);
prior.d         = d;
prior.chol2w    = prior_cov_lt(prior, prior.basis_w);
prior.w2chol    = prior_cov_il(prior, prior.basis);

prior.type      = 'Basis';
prior.note      = 'KL';

%figure
%semilogy(cumsum(d)/sum(d))
%title('Prior cummulative energy')

end