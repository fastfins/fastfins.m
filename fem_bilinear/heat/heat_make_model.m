function FEM = heat_make_model(mesh, options)
%HEAT_MAKE_OUTPUTS
%
% generates measurements points and indices w.r.t. the FEM mesh
%
% Tiangang Cui, 03/May/2014

tol = 1e-10;

FEM = heat_make_FEM(mesh); % setup the basic structure

f   = options.force_func(mesh.node);
FEM.force   = FEM.M*f;
FEM.fs      = FEM.force;

% set the enssential boundary condition
FEM.c = spalloc(FEM.DoF, FEM.DoF, FEM.DoF);
if sum(options.ess_bc_flag) > 0
    for i = 1:4
        FEM.c = FEM.c + FEM.bnd_mat{i}*options.ess_bc_flag(i)*FEM.penalty;
        FEM.fs(FEM.bnd_ind{i}) = FEM.fs(FEM.bnd_ind{i}) + options.ess_bc_vals(i)*FEM.penalty;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
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
%}
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
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FEM.Tstart      = options.obs_tstart;
FEM.Tfinal      = options.obs_tfinal;

FEM.Ndatasets   = options.obs_ntime;

Tsteps          = [logspace(log10(options.tinit), log10(options.tnormal), options.ntransit), ...
    options.tnormal*ones(1, options.nnormal)]; % setup the time

FEM.Tsteps      = Tsteps/sum(Tsteps)*FEM.Tfinal;
FEM.Tsteps(end) = FEM.Tsteps(end)+FEM.Tsteps(1)*1E-3;
FEM.Nsteps      = length(FEM.Tsteps);

% process observations
obs_Tsteps      = linspace(0, FEM.Tfinal, FEM.Ndatasets);
FEM.obs_int     = sparse([],[],[], FEM.Nsteps, FEM.Ndatasets, FEM.Ndatasets*2);
t_start         = 0;
for           i = 1:FEM.Nsteps
    t_end       = t_start + FEM.Tsteps(i);
    t_ind       = find(obs_Tsteps>t_start & obs_Tsteps<=t_end);
    temp1       = (t_end - obs_Tsteps(t_ind))/FEM.Tsteps(i);   % weighting at t_start
    temp2       = (obs_Tsteps(t_ind) - t_start)/FEM.Tsteps(i); % weighting at t_end
    t_start     = t_end; % increment the time
    
    FEM.obs_int(i,  t_ind)  = temp1;
    FEM.obs_int(i+1,t_ind)  = temp2;
end

end