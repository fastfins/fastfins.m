function FEM = RD_make_model(mesh, options)
%HEAT_MAKE_OUTPUTS
%
% generates measurements points and indices w.r.t. the FEM mesh
%
% Tiangang Cui, 03/May/2014

tol = 1E-10;
FEM = RD_make_FEM(mesh); % setup the basic structure
f   = options.force_func(mesh.node);
FEM.fs = FEM.M*f;

FEM.fast_adj_flag   = options.fast_adj;
FEM.ML_flag         = options.use_ML;
FEM.a               = options.param_a;
FEM.kappa           = options.param_k;
FEM.dt_max          = options.dt_max;
FEM.dt_init         = options.dt_init;
FEM.dt_multiplier   = options.dt_multiplier;
FEM.Nmaxsteps       = 500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch options.obs_type
    case{'Nodes'}
        N     = size(options.obs_locs,1);
        FEM.C = [];
        for k = 1:N
            ind = mesh.centers(1,:)>(options.obs_locs(k,1)+tol) & mesh.centers(1,:)<(options.obs_locs(k,2)-tol) & ...
                mesh.centers(2,:)>(options.obs_locs(k,3)+tol) & mesh.centers(2,:)<(options.obs_locs(k,4)-tol);
            tmp = 1:mesh.Nel;
            elems = tmp(ind);
            C   = sparse([],[],[],1,mesh.Nnode,ceil((sqrt(sum(ind))+1)^2));
            for i = 1:sum(ind)
                nodes = mesh.node_map(:,elems(i));
                dx    = mesh.node(:,nodes(3)) - mesh.node(:,nodes(1));
                detJ  = prod(abs(dx)); % iJ = diag(1./dx);
                locs  = 0.25*detJ*ones(1,4);
                C(1,nodes)  = C(1,nodes) + locs;
            end
            FEM.C = [FEM.C; C];
        end
        FEM.Nsensors = size(FEM.C, 1);
        
    case{'Interp'}
        FEM.Nsensors = size(options.obs_locs,1);
        FEM.C        = sparse([],[],[], FEM.Nsensors, mesh.Nnode, FEM.Nsensors*4); % make the observation matrix
        
        bl = mesh.node(:, mesh.node_map(1,:));
        %br = mesh.node(:, mesh.node_map(2,:));
        tr = mesh.node(:, mesh.node_map(3,:));
        %tl = mesh.node(:, mesh.node_map(4,:));
        
        FEM.sensor_ind = zeros(FEM.Nsensors,4);
        FEM.sensor_wei = zeros(FEM.Nsensors,4);
        for i = 1:FEM.Nsensors
            tmpx = options.obs_locs(i,1) + tol;
            tmpy = options.obs_locs(i,2) + tol;
            ind  = bl(1,:) <= tmpx & bl(2,:) <= tmpy & tr(1,:) > tmpx & tr(2,:) > tmpy;
            if sum(ind) ~= 1
                disp([i tmpx tmpy sum(ind)])
                error('cannot locate sensor')
            end
            FEM.sensor_ind(i,:) = mesh.node_map(:,ind);
            FEM.C(i,FEM.sensor_ind(i,:)) = bilinear_weights_rec(bl(:,ind), tr(:,ind), [tmpx;tmpy]);
        end
end

FEM.Tstart      = options.Tstart;
FEM.Tfinal      = options.Tfinal;
FEM.Ndatasets   = options.Ntime;

end