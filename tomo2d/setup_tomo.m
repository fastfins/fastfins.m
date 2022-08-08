function [model, obs, prior] = setup_tomo(tomo_opts, inv_opts, output)
%SETUP_TOMO
%
% Tiangang Cui, 11/August/2022

model = build_model(tomo_opts.mesh_size, tomo_opts.angle, ...
    tomo_opts.n_sourc, tomo_opts.n_detec, tomo_opts.d_width);
model.Is    = tomo_opts.Is;
% transformation
model.func.type = tomo_opts.f_type;
model.func.log_thres = tomo_opts.log_thres;
model.func.erf_scale = tomo_opts.erf_scale;
model.func.erf_shift = tomo_opts.erf_shift;

%%%%

dh      = 0.5/tomo_opts.mesh_size;
gx      = linspace(dh, 1-dh, tomo_opts.mesh_size);
gy      = linspace(dh, 1-dh, tomo_opts.mesh_size);
mesh    = rectmesh2d(gx, gy); % mesh generator

% prior
switch  inv_opts.cov_type
    case {'GP'}
        prior = make_prior_gp  (mesh, inv_opts.scale, inv_opts.power, inv_opts.sigma);
    case {'MRF'}
        prior = make_prior_mrf (mesh, inv_opts.k, inv_opts.cond, inv_opts.sigma);
    case {'Conv'}
        prior = make_prior_conv(mesh, inv_opts.centers, inv_opts.radii, inv_opts.weights, inv_opts.kernel_func);
end
prior.NP  = prior.mesh.Nnode;

if length(inv_opts.mean) == 1
    prior.mean_u = inv_opts.mean*ones(prior.NP,1);
else
    prior.mean_u = inv_opts.mean;
end

%%%%%%%%%%%%%%%%%%% loading data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obs.like = inv_opts.like;
obs.jacobian = inv_opts.jacobian;

if  ~isempty(output) && exist(output,'file')
    tmp         = load(output);
    obs.data    = tmp.data;
    prior.truex = tmp.truex;
else
    trueu = load_test_image(prior, inv_opts.test_type, inv_opts.test_base, inv_opts.test_range);
    %truex = inv_opts.test_range - truex + inv_opts.test_base*2;
    sol = forward_solve(model, trueu);
    
    data  = poissrnd(sol.d);
    obs.data  = data;
    prior.trueu = trueu;
end

if  ~isempty(output) && ~exist(output,'file')
    save(output, 'data', 'trueu');
end


obs.n_sensors    = length(obs.data);
obs.n_datasets   = 1;
obs.n_data       = obs.n_sensors*obs.n_datasets;
obs.like         = 'poisson';


%%%%%%%%%%%%%%%%%%% end of loading data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end