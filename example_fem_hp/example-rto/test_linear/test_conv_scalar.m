% test squared parameter
sq_flag = false;

% test regular mesh
reg_flag = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test scalar kappa field
[state_f, kappa_f, force_f, bottom_flux, bnd_funs, bc_types, bc_funs] = ref_soln2d_scalar(true);

% flux function for qoi
flux_fx = @(x) uniform_fun (x, 0, 1); 
flux_fy = @(y) triangle_fun(y, 0, 1, 'right');
flux = @(x1, x2) flux_fx(x1).*flux_fy(x2);

% obs locations
gx  = linspace(1/4, 3/4, 4);
gy  = linspace(1/4, 3/4, 4);
[xx0,yy0] = meshgrid(gx, gy);
xx0 = xx0(:);
yy0 = yy0(:);

gx  = linspace(0.75, 0.95, 2);
gy  = linspace(0.05, 0.25, 2);
[xx1,yy1] = meshgrid(gx, gy);
xx1 = xx1(:);
yy1 = yy1(:);

obs_locs  = unique([xx0, yy0; xx1 yy1], 'rows');     

%
ps = 1:4;
hs = 1./[10, 20, 40, 80];
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% test general mesh

model_opts = hp_options('exp_param', true, 'h', hs(end), 'poly_order', ps(end)+1, 'quad_order', 20, ...
    'xyratio', 1, 'force_func', force_f, 'qoi_func', flux, 'obs_locs', obs_locs, 'sq_param', sq_flag, ...
    'bnd_funcs', bnd_funs,'bc_types', bc_types, 'bc_funcs', bc_funs);

% high resolution
model = build_model(model_opts, reg_flag);


% models and priors
models = cell(length(hs), length(ps));
for hi = 1:length(hs)
    for pi = 1:length(ps)
        model_opts = hp_options(model_opts, 'h', hs(hi), 'poly_order', ps(pi), 'quad_order', ps(pi)*3+1);
        models{hi,pi} = build_model(model_opts, reg_flag);
    end
end

% full solution
k = kappa_f(model.mesh.nodes(:,1), model.mesh.nodes(:,2));
if sq_flag
    sol = forward_solve(model, 0.5*log(k));
else
    sol = forward_solve(model, log(k));
end
ref = state_f(model.mesh.nodes(:,1), model.mesh.nodes(:,2));
disp((ref-sol.state)'*model.mass*(ref-sol.state))
figure
trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), sol.state, 'edgecolor', 'none')

% hp solution
sols = cell(length(hs), length(ps));
L2_error = zeros(length(hs), length(ps));
L2_qoi   = zeros(length(hs), length(ps));
L2_obs   = zeros(length(hs), length(ps));
dofs = zeros(length(hs), length(ps));
for hi = 1:length(hs)
    for pi = 1:length(ps)
        k = kappa_f(models{hi,pi}.mesh.nodes(:,1), models{hi,pi}.mesh.nodes(:,2));
        if sq_flag
            sols{hi,pi} = forward_solve(models{hi,pi}, 0.5*log(k));
        else
            sols{hi,pi} = forward_solve(models{hi,pi}, log(k));
        end
        ref = state_f(models{hi,pi}.mesh.nodes(:,1), models{hi,pi}.mesh.nodes(:,2));
        %
        L2_error(hi,pi) = sqrt(((ref-sols{hi,pi}.state)'*models{hi,pi}.mass*(ref-sols{hi,pi}.state)));
        L2_qoi  (hi,pi) = abs(sol.qoi-sols{hi,pi}.qoi);
        L2_obs  (hi,pi) = norm(sol.d-sols{hi,pi}.d);
        dofs(hi,pi) = models{hi,pi}.mesh.dof;
    end
end

figure
%loglog(sqrt(dofs), L2_error)
loglog(1./hs, L2_error)
title('norm(u)')
figure
%loglog(sqrt(dofs), L2_qoi)
loglog(1./hs, L2_qoi)
title('norm(qoi)')
figure
%loglog(sqrt(dofs), L2_obs)
loglog(1./hs, L2_obs)
title('norm(d)')

disp('rate l2')
disp( (log2(L2_error(end,:))-log2( L2_error(1,:)))/(log2(hs(end)) - log2(hs(1))) )
disp('rate obs')
disp( (log2(L2_obs(end,:))-log2( L2_obs(1,:)) )/(log2(hs(end)) - log2(hs(1))) )
disp('rate qoi')
disp( (log2(L2_qoi(end,:))-log2( L2_qoi(1,:)) )/(log2(hs(end)) - log2(hs(1))) )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


