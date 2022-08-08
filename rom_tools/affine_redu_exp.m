function redu_x = affine_redu_exp(prior_redu, sub_samples, opts)

xs  = exp( prior_redu.basis*sub_samples + prior_redu.mean_u(:));

[U,S,~] = svd(xs, 'econ');
d   = diag(S).^2;
ii  = d >= d(2)*opts.prior_deim_tol;
nU  = sum(ii);

[redu_x.e, redu_x.B, redu_x.K] = deim(U, nU, opts.deim_reg_factor);

redu_x.LU   = prior_redu.basis(redu_x.e, :);
redu_x.umu  = prior_redu.mean_u(redu_x.e);

redu_x.d    = d;
redu_x.U    = U(:,1:nU);

disp(sprintf('DEIM dim = %5i', sum(ii)));

end