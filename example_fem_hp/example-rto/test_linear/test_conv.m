% test squared parameter
sq_flag = true;

% test regular mesh
reg_flag = true;

% test case
test_case = 'scalar';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

renew_model = false;

% test scalar kappa field
switch test_case
    case {'scalar'}
        [state_f, kappa_f, force_f, bottom_flux, bnd_funs, bc_types, bc_funs] = ref_soln2d_scalar(true);
    case {'vector'}
        [state_f, kappa_f, force_f, bottom_flux, bnd_funs, bc_types, bc_funs] = ref_soln2d_vector(true);
    case {'tensor'}
        [state_f, kappa_f, force_f, bottom_flux, bnd_funs, bc_types, bc_funs] = ref_soln2d_tensor(true);
end

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
hs = 1./[8, 16, 32, 64];
%


% test general mesh

model_opts = hp_options('exp_param', true, 'h', hs(end)/2, 'poly_order', ps(end), 'quad_order', 20, ...
    'xyratio', 1, 'force_func', force_f, 'qoi_func', flux, 'obs_locs', obs_locs, 'sq_param', sq_flag, ...
    'bnd_funcs', bnd_funs,'bc_types', bc_types, 'bc_funcs', bc_funs);

% high resolution
model = build_model(model_opts, reg_flag);


% models and priors
models = cell(length(hs), length(ps));
for hi = 1:length(hs)
    for pi = 1:length(ps)
        model_opts = hp_options(model_opts, 'h', hs(hi), 'poly_order', ps(pi), 'quad_order', ps(pi)*4+1);
        models{hi,pi} = build_model(model_opts, reg_flag);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% full solution
switch test_case
    case {'scalar'}
        u = kappa_f(model.mesh.nodes(:,1), model.mesh.nodes(:,2));
        if sq_flag
            u = 0.5*log(u);
        else
            u = log(u);
        end
    case {'vector'}
        nd = model.local_elem.num_dim;
        u = zeros(model.mesh.dof, nd);
        for di = 1:nd
            u(:,di) = kappa_f{di}(model.mesh.nodes(:,1), model.mesh.nodes(:,2));
        end
        if sq_flag
            u = 0.5*log(u);
        else
            u = log(u);
        end
    case {'tensor'}
        nd = model.local_elem.num_dim;
        u = zeros(model.mesh.dof, nd+1);
        u(:,1) = kappa_f{1}(model.mesh.nodes(:,1), model.mesh.nodes(:,2));
        if sq_flag
            u(:,1) = 0.5*log(u(:,1));
        else
            u(:,1) = log(u(:,1));
        end
        for di = 1:nd
            u(:,di+1) = kappa_f{di+1}(model.mesh.nodes(:,1), model.mesh.nodes(:,2));
        end
end

%
sol = forward_solve(model, u);
ref = state_f(model.mesh.nodes(:,1), model.mesh.nodes(:,2));
disp((ref-sol.state)'*model.mass*(ref-sol.state))
figure
trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), sol.state, 'edgecolor', 'none')
figure
trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2),ref, 'edgecolor', 'none')
L2 = sqrt(sol.state'*model.mass*sol.state);
H1 = sqrt(sol.state'*model.stiff*sol.state + sol.state'*model.mass*sol.state);

% hp solution
sols = cell(length(hs), length(ps));
L2_error = zeros(length(hs), length(ps));
L2_qoi   = zeros(length(hs), length(ps));
L2_obs   = zeros(length(hs), length(ps));
dofs = zeros(length(hs), length(ps));
for hi = 1:length(hs)
    for pi = 1:length(ps)
        switch test_case
            case {'scalar'}
                u = kappa_f(models{hi,pi}.mesh.nodes(:,1), models{hi,pi}.mesh.nodes(:,2));
                if sq_flag
                    u = 0.5*log(u);
                else
                    u = log(u);
                end
            case {'vector'}
                nd = models{hi,pi}.local_elem.num_dim;
                u = zeros(models{hi,pi}.mesh.dof, nd);
                for di = 1:nd
                    u(:,di) = kappa_f{di}(models{hi,pi}.mesh.nodes(:,1), models{hi,pi}.mesh.nodes(:,2));
                end
                if sq_flag
                    u = 0.5*log(u);
                else
                    u = log(u);
                end
            case {'tensor'}
                nd = models{hi,pi}.local_elem.num_dim;
                u = zeros(models{hi,pi}.mesh.dof, nd+1);
                u(:,1) = kappa_f{1}(models{hi,pi}.mesh.nodes(:,1), models{hi,pi}.mesh.nodes(:,2));
                if sq_flag
                    u(:,1) = 0.5*log(u(:,1));
                else
                    u(:,1) = log(u(:,1));
                end
                for di = 1:nd
                    u(:,di+1) = kappa_f{di+1}(models{hi,pi}.mesh.nodes(:,1), models{hi,pi}.mesh.nodes(:,2));
                end
        end
        %
        sols{hi,pi} = forward_solve(models{hi,pi}, u);
        ref = state_f(models{hi,pi}.mesh.nodes(:,1), models{hi,pi}.mesh.nodes(:,2));
        %
        L2_error(hi,pi) = sqrt(((ref-sols{hi,pi}.state)'*models{hi,pi}.mass*(ref-sols{hi,pi}.state)));
        L2_diff(hi,pi)  = abs(sqrt(sols{hi,pi}.state'*models{hi,pi}.mass*sols{hi,pi}.state) - L2);
        %
        H1_error(hi,pi) = sqrt( (ref-sols{hi,pi}.state)'*models{hi,pi}.stiff*(ref-sols{hi,pi}.state) + (ref-sols{hi,pi}.state)'*models{hi,pi}.mass*(ref-sols{hi,pi}.state) );
        H1_diff(hi,pi)  = abs(sqrt(sols{hi,pi}.state'*models{hi,pi}.mass*sols{hi,pi}.state + sols{hi,pi}.state'*models{hi,pi}.stiff*sols{hi,pi}.state) - H1);
        %
        L2_qoi  (hi,pi) = abs(sol.qoi-sols{hi,pi}.qoi);
        L2_obs  (hi,pi) = norm(sol.d-sols{hi,pi}.d);
        %
        L2_qoi2 (hi,pi) = abs(mean(sol.d)-mean(sols{hi,pi}.d));
        %
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


disp( (log2(L2_diff(end,:))-log2( L2_diff(1,:)) )/(log2(hs(end)) - log2(hs(1))) )
disp( (log2(L2_qoi2(end,:))-log2( L2_qoi2(1,:)) )/(log2(hs(end)) - log2(hs(1))) )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


