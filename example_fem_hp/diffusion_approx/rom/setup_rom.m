function rom = setup_rom(model, states, prior_redu, sub_vs, weights, rom_opts)
%set up the reduced order model for the 2nd order PDE
%

%
[U,S] = svd((model.mass_L'*states).*sqrt(weights(:)'), 'econ');
s  = diag(S);
nU = truncate_energy(s(2:end).^2, rom_opts.pod_state_tol) + 1; % take out the first basis from the energy
U  = model.mass_L'\U(:,1:nU);
rom.dof = nU;
%
nd = model.local_elem.num_dim;

% build reduced basis for the params
n = length(prior_redu.bases);

% DEIM for mu_a
mus = exp(prior_redu.bases*sub_vs + prior_redu.mean_us(:));
[V,S,~] = svd((model.mass_L'*mus).*sqrt(weights(:)'), 'econ');
V  = model.mass_L'\V;
d  = diag(S).^2;
ii = d >= d(2)*rom_opts.prior_deim_tol;
nP = sum(ii);
V_mu = V(:,1:nP);
%
[rom.e_mu,~,rom.K_mu] = deim(V, nP, rom_opts.deim_reg_factor);
rom.redu_u_mu = prior_redu.bases(rom.e_mu,:);
rom.redu_mean_mu = prior_redu.mean_us(rom.e_mu);

%
% DEIM for kappa
ks = 1./(2*(mus + 0.2*model.mu_s));
[V,S,~] = svd((model.mass_L'*ks).*sqrt(weights(:)'), 'econ');
V  = model.mass_L'\V;
d  = diag(S).^2;
ii = d >= d(2)*rom_opts.prior_deim_tol;
nP = sum(ii);
V_k = V(:,1:nP);
%
[rom.e_k,~,rom.K_k] = deim(V, nP, rom_opts.deim_reg_factor);
rom.redu_u_k = prior_redu.bases(rom.e_k,:);
rom.redu_mean_k = prior_redu.mean_us(rom.e_k);

[rom.Ms,rom.Ks] = build_rom_matrices(model, U, V_mu, V_k);

% other matrices, rhs, and observation operators
rom.Abnd = U'*( model.M_bnd*U )*(2/pi);
rom.b = U'*model.b;
%
rom.states = U;
rom.param_redu = prior_redu.bases;
rom.mu_s = model.mu_s;
end