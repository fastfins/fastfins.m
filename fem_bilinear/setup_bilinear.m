function [model, obs, prior, soln] = setup_bilinear(model_opt, prior_opt, output)
%SETUP_PDE
%
% Tiangang Cui, 11/May/2014

gx              = linspace(0, model_opt.xyratio, model_opt.xyratio*model_opt.mesh_size+1);
gy              = linspace(0, 1, model_opt.mesh_size+1);
[fmesh,pmesh]   = rectmesh2d(gx,gy);  % mesh maker
fmesh.xyratio   = model_opt.xyratio;
pmesh.xyratio   = model_opt.xyratio;

% model
switch model_opt.test_case
    case {'EIT'} % EIT
        model = EIT_make_model(fmesh, model_opt);
    case {'Laplace'} % test case for DILI
        model = laplace_make_model(fmesh, model_opt);
    case{'Heat'} % test case 1
        model = heat_make_model(fmesh, model_opt);
    case{'RD'}   % test case 0
        model = RD_make_model(fmesh);
end
model.GMRES_flag  = model_opt.gmres;
model.problem     = model_opt.test_case;
model.explicitJ   = false;
model.beta        = model_opt.beta;

% prior
switch  prior_opt.cov_type
    case {'GP'}
        prior = make_prior_gp  (pmesh, prior_opt.scale, prior_opt.power, prior_opt.sigma);
    case {'MRF'}
        prior = make_prior_mrf (pmesh, prior_opt.k, prior_opt.cond, prior_opt.sigma);
    case {'Conv'}
        prior = make_prior_conv(pmesh, prior_opt.centers, prior_opt.radii, prior_opt.weights, prior_opt.kernel_func);
end
prior.NP  = prior.mesh.Nnode;

if length(prior_opt.mean) == 1
    prior.mean_u = prior_opt.mean*ones(prior.NP,1);
else
    prior.mean_u = prior_opt.mean;
end

% transformation
model.func.type = model_opt.f_type;
model.func.log_thres = model_opt.log_thres;
model.func.erf_scale = model_opt.erf_scale;
model.func.erf_shift = model_opt.erf_shift;

%%%%%%%%%%%%%%%%%%% loading data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if  ~isempty(output) && exist(output,'file')
    tmp         = load(output);
    obs.data    = tmp.data;
    obs.std     = tmp.std;
    prior.truex = tmp.truex;
else
    truex = load_test_image(prior, model.func, prior_opt.test_type, prior_opt.test_base, prior_opt.test_range);
    soln  = forward_solve(model, truex);  % reference solution
    
    % the s.t.d. is calculated from the signal to noise ratio
    if model_opt.s2n > 0
        std = mean(abs(soln.d(:)))/model_opt.s2n;
    else
        std = model_opt.std;
    end
    
    % generate data
    data    = soln.d + randn(model.Nsensors, model.Ndatasets)*std;
    
    obs.data    = data;
    obs.std     = std;
    prior.truex = truex;
end

if  ~isempty(output) && ~exist(output,'file')
    save(output, 'data', 'std', 'truex');
end

obs.Nsensors    = model.Nsensors;
obs.Ndatasets   = model.Ndatasets;
obs.Ndata       = obs.Nsensors*obs.Ndatasets;

%%%%%%%%%%%%%%%%%%% end of loading data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end