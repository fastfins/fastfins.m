function rom = setup_p_poisson_rom(model, prior_redu, sub_vs, states, weights, rom_opts)
%set up the reduced order model for the 2nd order PDE
%
%inputs:

ne = model.mesh.num_elem;
nd = model.local_elem.num_dim;
nb = model.local_elem.num_node;
nq = model.local_elem.num_quad;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
[U,S] = svd((model.mass_L'*states).*sqrt(weights(:)'), 'econ');
s  = diag(S);
nU = truncate_energy(s(2:end).^2, rom_opts.pod_state_tol) + 1; % take out the first basis from the energy
rom.states = model.mass_L'\U(:,1:nU);
rom.dof = nU;
%
grad_s = cell(nd, 1);
for di = 1:nd
    grad_s{di} = zeros(nq*ne, nU);
end
%
for i = 1:nU
    gs = calc_grad(model.mesh, model.local_elem, rom.states(:,i));
    for di = 1:nd
        grad_s{di}(:,i) = gs{di}(:);
    end
end
               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% build function basis for all states (not on the POD basis)
ns = size(states, 2);
sqrt_etas = zeros(nq*ne, ns);
bs = zeros(nq*ne*nd, ns);
% build DEIM bases on quadrature points
for i = 1:ns
    [fstate, gs] = eval_pre_data(model, states(:,i));
    % sqrt of eta
    sqrt_etas(:,i) = fstate(:).^(model.p_rate/4-1/2);
    % vector b
    tmpb = fstate(:).^(model.p_rate/4-1);
    for di = 1:nd
        ind = (1:nq*ne) + (di-1)*nq*ne;
        bs(ind,i) = tmpb.*gs{di}(:);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEIM for eta 
[Veta,S,~] = svd(sqrt(model.WdetJ(:)).*sqrt_etas, 'econ');
d  = diag(S).^2;
ii = d >= d(2)*rom_opts.deim_state_tol;
neta = sum(ii);
Veta = Veta(:,1:ceil(neta*rom_opts.deim_reg_factor))./sqrt(model.WdetJ(:));
[rom.eta_e,~,rom.eta_K] = deim(Veta, neta, rom_opts.deim_reg_factor);
Veta = Veta(:,1:neta);
%
% build interpolation basis for sqrt_eta
%rom.eta_ei = ceil(rom.eta_e/nq);
%rom.eta_qi = rom.eta_e - (rom.eta_ei-1)*nq;
rom.eta_ne = length(rom.eta_e);
rom.eta_U  = zeros(rom.eta_ne*nd, nU);
for di = 1:nd
    ind = (1:rom.eta_ne) + (di-1)*rom.eta_ne;
    rom.eta_U(ind,:) = grad_s{di}(rom.eta_e,:);
end

% DEIM for b
[Vb,S,~] = svd(repmat(sqrt(model.WdetJ(:)),nd,1).*bs, 'econ');
d  = diag(S).^2;
ii = d >= d(2)*rom_opts.deim_state_tol;
nb = sum(ii);
Vb = Vb(:,1:ceil(nb*rom_opts.deim_reg_factor))./repmat(sqrt(model.WdetJ(:)),nd,1);
[rom.b_e,~,rom.b_K] = deim(Vb, nb, rom_opts.deim_reg_factor);
Vb = Vb(:,1:nb);
%
rom.b_U  = cell(nd, 1);
rom.b_ne = zeros(nd, 1);
%rom.b_ei = cell(nd, 1);
%rom.b_qi = cell(nd, 1);
for di = 1:nd
    % build the interpolation basis for b
    ind = ( rom.b_e > (nq*ne*(di-1)) ) & ( rom.b_e <= (nq*ne*di) );
    be_di = rom.b_e(ind) - nq*ne*(di-1);
    %rom.b_ei{di} = ceil(be_di/nq);
    %rom.b_qi{di} = be_di - (rom.b_ei{di}-1)*nq;
    rom.b_ne(di) = length(be_di);
    %
    rom.b_U{di} = zeros(rom.b_ne(di)*nd, nU);
    for dj = 1:nd
        ind = (1:rom.b_ne(di)) + (dj-1)*rom.b_ne(di);
        rom.b_U{di}(ind,:) = grad_s{dj}(be_di,:);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DEIM for parameter
xs = sqrt(model.beta_weight).*exp(0.5*(prior_redu.basis*sub_vs + prior_redu.mean_u(:)));
[Vx,S,~] = svd(xs.*sqrt(weights(:)'), 'econ');
d  = diag(S).^2;
ii = d >= d(2)*rom_opts.deim_state_tol;
nx = sum(ii);
Vx = Vx(:,1:ceil(nx*rom_opts.deim_reg_factor));
[rom.x_e,~,rom.x_K] = deim(Vx, nx, rom_opts.deim_reg_factor);
rom.beta_weight = sqrt(model.beta_weight(rom.x_e,:));
rom.redu_u = prior_redu.basis(rom.x_e,:);
rom.redu_mean = prior_redu.mean_u(rom.x_e);
Vx = Vx(:,1:nx);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rom.As, rom.inds] = build_poisson_rom_matrices(model, grad_s, Veta, Vb, rom.states, Vx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rhs and observation operators
rom.b = rom.states'*model.b;
rom.obs_operator = model.obs_operator*rom.states;
%
rom.epsilon = model.epsilon;
rom.p_rate = model.p_rate;
rom.nd = model.local_elem.num_dim;
%
rom.res_tol = model.res_tol;
rom.iter_disp = model.iter_disp;
rom.max_newton_iter = model.max_newton_iter;

end