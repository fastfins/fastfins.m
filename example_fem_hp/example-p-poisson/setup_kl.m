load_dir;

renew_setup = false;

if isfile('setup.mat') && ~renew_setup
    load('setup.mat');
else
    p_rate = 3;
    epsilon = 1E-6;
    tol = 1E-4;
    bnd_funcs = {@(p) find(abs(p(:,2))<tol); @(p) find(abs(p(:,1))<tol)};
    %
    obs_locs    = [linspace(0.1,0.9,20)', 0.05*ones(20,1)];
    force_func  = @(x1, x2) 1;
    %
    model_opts = hp_options('h', 1/50, 'poly_order', 2, 'quad_order', 10, ...
        'xyratio', 0.2, 'obs_locs', obs_locs, 'bnd_funcs', bnd_funcs, 'force_func', force_func, ...
        'nl_solver', 'line_search', 'res_tol', 1E-10, 'exp_param', true, 'nl_disp', 1, ...
        'num_newton', 100); 
    % set 'num_newton' for the number of Newton iterations
    % set 'nl_disp' for display options: 0 -- nothing except warning 
    %                                    1 -- display at final newton iteration
    %                                    2 -- display at each newton iteration
    %
%     inv_opts = inverse_options('s2n', 50, 'cov_type', 'GP', 'mean', 0, ...
%         'power', 2, 'scale', 200, 'sigma', 1);
    inv_opts = inverse_options('s2n', 10, 'cov_type', 'GP', 'mean', 0, ...
        'power', 2, 'scale', 200, 'sigma', 1);    
    %
    [model, obs, prior] = setup_p_poisson(model_opts, inv_opts, p_rate, epsilon, []);
    
    % build kl basis
    prior_kl = basis_KL(prior, 1-1E-3);
    
    solt = forward_solve(model, obs.true_u);
    figure
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), solt.state, 'edgecolor', 'none')
    figure
    plot(model.mesh.nodes(model.pa_bnd_mesh.node_indices, 1), obs.true_u)
    
    
    mlpost  = @(v) minus_log_post(model, obs, prior_kl, v);
    mvhess  = @(sol, dv) matvec_PPFisher(model, prior_kl, sol, dv);
    vini = randn(prior_kl.dof, 1);
    vmap = get_map_2020a(mlpost, mvhess, vini);
    umap = matvec_prior_L(prior_kl, vmap) + prior_kl.mean_u;
    solm = forward_solve(model, umap);
    
    figure
    trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), solm.state, 'edgecolor', 'none')
    figure
    plot(model.mesh.nodes(model.pa_bnd_mesh.node_indices, 1), umap)
    
    figure
    plot(obs.data, 'o')
    hold on
    plot(solm.d, 'x')
    plot(solt.d, '-')
    legend('obs', 'map sol', 'true sol');
   
    save('setup.mat');
end

% test the qoi, random samples
tic;
for i = 1:100
    u = matvec_prior_L(prior_kl, randn(prior_kl.dof,1)) + prior.mean_u;
    sol = forward_solve(model, u);
    qois(i) = sol.qoi;
end
toc

% test the qoi, random samples on the corners of the hypercube [-4, 4]
tic;
for i = 1:100
    v = 4 * ((randn(prior_kl.dof,1) > 0)*2-1 );
    u = matvec_prior_L(prior_kl, v) + prior.mean_u;
    sol = forward_solve(model, u);
    qois(i) = sol.qoi;
end
toc

xstate = repmat({(-4+8/16:8/16:4)'}, prior_kl.dof, 1);
% xdata  = repmat({linspace(0.9, 1.2, 15)'}, numel(obs.data), 1);

% lpfun = @(x,b1,b2)-minus_log_post_temp(model, ...
%                         setfield(obs, 'data', reshape(x(1:numel(obs.data)), [], 1)), ...
%                         prior_kl, ...
%                         reshape(x(numel(obs.data)+1:end), prior_kl.dof, 1), ...
%                         b1, b2);

lpfun = @(x,b1,b2)-minus_log_post_temp(model, ...
                        obs, ...
                        prior_kl, ...
                        x', ...
                        b1, b2)';


IRT = tt_dirt_approx(xstate, lpfun, 10.^(-4:1/2:0), ...
                     'nswp', 1, 'kickrank', 0, 'y0', 10, ...
                     'reference', 'n4', 'interpolation', 'f', 'boundary', true);

