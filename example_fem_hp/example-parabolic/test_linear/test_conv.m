%load_dir;

test_case = 'tensor';
sq_flag = true;

init_func = @(x1,x2) (cos(x1.*pi.*2)-1).*sin(x2*pi*2)*2;

% flux function for qoi
flux_fx = @(x) uniform_fun (x, 0, 1);
flux_fy = @(y) triangle_fun(y, 0, 1, 'right');
flux = @(x1, x2) flux_fx(x1).*flux_fy(x2);

% obs locations
gx  = linspace(1/4, 3/4, 2);
gy  = linspace(1/4, 3/4, 2);
[xx0,yy0] = meshgrid(gx, gy);
xx0 = xx0(:);
yy0 = yy0(:);

obs_locs  = unique([xx0, yy0; 0.8, 0.2], 'rows');

dh = 1/20;
dt = 1/100;

hs = 2.^(0:3);
for i = 1:length(hs)
    t_fac = 1/round(1/(1/hs(i))^(3/2));
    model_opts = hp_options('h', dh/hs(i), 'poly_order', 2, 'quad_order', 15, ...
        'xyratio', 1, 'qoi_func', flux, 'obs_locs', obs_locs, 'sq_param', sq_flag, ...
        'time_dep', true, 'dt', dt*t_fac, 'init_func', init_func, ...
        'obs_tstart', 0.05, 'obs_tfinal', 0.15, 'obs_ntime', 6, 'pred_t', 0.2);
    
    % setup the process convolution prior
    centers_x = 0.1:0.2:0.9;
    centers_y = 0.1:0.2:0.9;
    [xx,yy] = meshgrid(centers_x, centers_y);
    centers = [xx(:), yy(:)];
    radii   = ones(size(centers, 1), 1)*0.1;
    weights = ones(size(centers, 1), 1)*1.5;
    k_func  = @(x) exp( - 0.5*x.^2);
    
    inv_opts = inverse_options('s2n', 10, 'cov_type', 'Conv', 'mean', 0, ...
        'centers', centers, 'radii', radii, 'weights', weights, ...
        'kernel_func', k_func, 'test_type', 'Prior');
    
    [models{i}, ~, priors{i}] = setup_parabolic(model_opts, inv_opts, test_case, []);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 5;
v = randn(priors{1}.dof,n);
qois = [];
sols = {};
for hi = 1:length(hs)
    tic;
    for j = 1:n
        %
        u = matvec_prior_L(priors{hi}, v(:,j)) + priors{hi}.mean_u;
        sol = forward_solve(models{hi}, u);
        sols{hi}(:,j) = sol.d(:);
        qois(hi,j) = sol.qoi;
    end
    toc
end
%
for hi = 1:length(hs)
    L2_qoi  (hi) = mean(abs(qois(length(hs),:)-qois(hi,:)));
    L2_obs  (hi) = mean(sqrt(sum((sols{length(hs)}-sols{hi}).^2, 1)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
%loglog(sqrt(dofs), L2_qoi)
loglog(hs, L2_qoi)
title('norm(qoi)')
figure
%loglog(sqrt(dofs), L2_obs)
loglog(hs, L2_obs)
title('norm(d)')

disp('rate obs')
disp( (log2(L2_obs(end-1))-log2( L2_obs(1)) )/(log2(hs(end-1)) - log2(hs(1))) )
disp('rate qoi')
disp( (log2(L2_qoi(end-1))-log2( L2_qoi(1)) )/(log2(hs(end-1)) - log2(hs(1))) )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


