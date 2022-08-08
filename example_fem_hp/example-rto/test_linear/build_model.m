function model = build_model(model_opts, reg_flag)

if reg_flag
    [p,t]   = simplemesh2d(model_opts.h, 1, model_opts.xyratio);
else
    func_mesh = @(x) drectangle0(x, 0, 1, 0, 1);
    [p,t]   = distmesh2d(func_mesh, @huniform, model_opts.h, [0, 0; 1, 1], []);
end
model   = setup_2nd_order(p, t, model_opts);

%
model.gmres_flag  = model_opts.gmres;
model.res_tol = model_opts.res_tol;
model.explicit_ja = false;

% assemble forcing term
force   = model_opts.force_func(model.mesh.nodes(:,1), model.mesh.nodes(:,2));
model.f = model.mass*force(:); % weak form of the force, on nodes

% assemble boundary condition
model.bnd_b = zeros(model.mesh.dof, 1);
model.bmass = spalloc(model.mesh.dof,model.mesh.dof,0);
for bi = 1:length(model_opts.bnd_funcs)
    b_ind   = model_opts.bnd_funcs{bi}(model.mesh.nodes);
    b_force = model_opts.bc_funcs{bi}(model.mesh.nodes(:,1),model.mesh.nodes(:,2));
    if length(b_force) == 1 
        b_force = ones(size(b_ind(:)))*b_force;
    else
        b_force = b_force(b_ind);
    end
    bf_weak = model.mass_bnds{bi}*sparse(b_ind(:),ones(size(b_ind(:))),b_force(:),model.mesh.dof, 1);
    switch model_opts.bc_types{bi}
        case {'flux'}
            disp('double check non-zero flux b.c.')
            model.bnd_b = model.bnd_b + bf_weak;
        case {'essential'}
            model.bmass = model.bmass + model.mass_bnds{bi}*model.penalty;
            model.bnd_b = model.bnd_b + bf_weak*model.penalty;
        case {'mixed'}
            disp('Mixed b.c. is not implemented')
    end
end

model.b = model.f + model.bnd_b;

% apply observation operator
model.obs_operator  = mesh_interpolate(model.mesh, model.local_elem, model_opts.obs_locs, model.h, false);
model.n_sensors     = size(model_opts.obs_locs,1);
model.n_datasets    = 1;

% apply qoi function
tol = 1E-10;
disp('flux QoI')
model.qoi_flag = false;
if ~isempty(model_opts.qoi_func)
    phi = reshape(model_opts.qoi_func(model.mesh.nodes(:,1),model.mesh.nodes(:,2)), 1, []);
    phi(abs(phi)<tol) = 0;
    model.phi = sparse(phi);
    model.qoi_flag  = true;
end

model.exp_param = model_opts.exp_param;
model.exp_thres = model_opts.exp_thres;
model.sq_param  = model_opts.sq_param;

end
