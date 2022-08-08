
%
p_rate = 3;
epsilon = 1E-6;
[state_func, eta_func, force_func, beta_func] = ref_soln2d(p_rate, epsilon, true);
%
tol = 1E-4;
bnd_funcs = {@(p) find(abs(p(:,2))<tol); @(p) find(abs(p(:,1))<tol)};
%
obs_locs    = [linspace(0.1,0.9,30)', 0.1*ones(30,1)];
%
% test regular mesh
reg_flag = true;    
%
ps = 1:4;
hs = 1./[10, 20, 40, 80];
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
model_opts = hp_options('h', hs(end)/2, 'poly_order', ps(end), 'quad_order', ps(end)*4, ...
    'xyratio', 1, 'obs_locs', obs_locs, 'bnd_funcs', bnd_funcs, 'force_func', force_func, ...
    'nl_solver', 'line_search', 'res_tol', 1E-10, 'exp_param', false);
%
model = build_model(model_opts, p_rate, epsilon, reg_flag);

% models and priors
models = cell(length(hs), length(ps));
for hi = 1:length(hs)
    for pi = 1:length(ps)
        
        model_opts = hp_options('h', hs(hi), 'poly_order', ps(pi), 'quad_order', ps(pi)*4, ...
            'xyratio', 1, 'obs_locs', obs_locs, 'bnd_funcs', bnd_funcs, 'force_func', force_func, ...
            'nl_solver', 'line_search', 'res_tol', 1E-10, 'exp_param', false);
        %
        models{hi,pi} = build_model(model_opts, p_rate, epsilon, reg_flag);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% full solution
k = beta_func(model.mesh.nodes(model.pa_bnd_mesh.node_indices,1));
sol = forward_solve(model, k);
ref = state_func(model.mesh.nodes(:,1), model.mesh.nodes(:,2));
disp((ref-sol.state)'*model.mass*(ref-sol.state))
figure
trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), sol.state, 'edgecolor', 'none')
figure
trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), ref, 'edgecolor', 'none')

% hp solution
sols = cell(length(hs), length(ps));
L2_error = zeros(length(hs), length(ps));
L2_obs   = zeros(length(hs), length(ps));
L2_qoi   = zeros(length(hs), length(ps));
dofs = zeros(length(hs), length(ps));
for hi = 1:length(hs)
    for pi = 1:length(ps)   
        k = beta_func(models{hi,pi}.mesh.nodes(models{hi,pi}.pa_bnd_mesh.node_indices,1));
        sols{hi,pi} = forward_solve(models{hi,pi}, k);
        ref = state_func(models{hi,pi}.mesh.nodes(:,1), models{hi,pi}.mesh.nodes(:,2));
        %
        L2_error(hi,pi) = sqrt(((ref-sols{hi,pi}.state)'*models{hi,pi}.mass*(ref-sols{hi,pi}.state)));
        L2_obs  (hi,pi) = norm(sol.d-sols{hi,pi}.d);
        L2_qoi  (hi,pi) = abs(sol.qoi-sols{hi,pi}.qoi);
        dofs(hi,pi) = models{hi,pi}.mesh.dof;
    end
end

figure
%loglog(sqrt(dofs), L2_error)
loglog(1./hs, L2_error)
title('norm(u)')
figure
%loglog(sqrt(dofs), L2_obs)
loglog(1./hs, L2_obs)
title('norm(d)')
figure
loglog(1./hs, L2_qoi)
title('norm(d)')

disp('rate l2')
disp( (log2(L2_error(end,:))-log2( L2_error(1,:)))/(log2(hs(end)) - log2(hs(1))) )
disp('rate obs')
disp( (log2(L2_obs(end,:))-log2( L2_obs(1,:)) )/(log2(hs(end)) - log2(hs(1))) )
disp('rate qoi')
disp( (log2(L2_qoi(end,:))-log2( L2_qoi(1,:)) )/(log2(hs(end)) - log2(hs(1))) )


