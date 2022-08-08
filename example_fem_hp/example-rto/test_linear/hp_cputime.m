

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
hp_list = [16 1;
32 1;
64 1;
128 1;
14 2;
22 2;
35 2;
56 2;
11 3;
15 3;
19 3;
25 3;
10 4;
13 4;
17 4;
23 4];
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model_opts = hp_options('exp_kappa', false, 'h', 1/20, 'poly_order', 5, 'quad_order', 20, ...
    'xyratio', 1, 'force_func', force_f, 'qoi_func', flux, 'obs_locs', obs_locs, ...
    'bnd_funcs', bnd_funs,'bc_types', bc_types, 'bc_funcs', bc_funs);

% models and priors
models = cell(size(hp_list,1), 1);
As = cell(size(hp_list,1), 1);
for i = 1:size(hp_list,1)
    model_opts = hp_options(model_opts, 'h', 1/hp_list(i,1), 'poly_order', hp_list(i,2), 'quad_order', hp_list(i,2)*3+1);
    models{i} = build_model_reg(model_opts);
    k = kappa_f(models{i}.mesh.nodes(:,1), models{i}.mesh.nodes(:,2));
    As{i} = build_stiff_matrix(models{i}, log(k));
end

for i = 1:size(hp_list,1)
    dofs(i) = size(As{i},1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_test = 1E2;
%time_lu = zeros(size(hp_list,1), 1);
time_ch = zeros(size(hp_list,1), 1);
time_lup = zeros(size(hp_list,1), 1);
for i = 1:size(hp_list,1)
    disp(i)
    bs = randn(models{i}.mesh.dof,n_test) + models{i}.b;
    us = zeros(models{i}.mesh.dof,n_test);
    
    tic;
    for j = 1:n_test
        [sol.L,~,sol.p] = chol(As{i},'lower', 'vector');
        % solve
        us(sol.p,j)  = sol.L'\(sol.L\bs(sol.p,j));
    end
    time_ch(i) = toc;
    
    tic;
    for j = 1:n_test
        [sol.L,sol.U] = lu(As{i}(sol.p,sol.p));
        us(sol.p,j)  = sol.U\(sol.L\bs(sol.p,j));
    end
    time_lup(i) = toc;
    %{
    tic;
    for j = 1:n_test
        [sol.L,sol.U] = lu(As{i});
        us(:,j)  = sol.U\(sol.L\bs(:,j));
    end
    time_lu(i) = toc;
    %}
end
