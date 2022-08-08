function FEM = laplace_make_model(mesh, options)
%LAPLACE_MAKE_OUTPUTS
%
% generates measurements points and indices w.r.t. the FEM mesh
%
% Tiangang Cui, 03/May/2014

tol = 1E-10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FEM = laplace_make_FEM(mesh);
f   = options.force_func(mesh.node);
FEM.force   = FEM.M*f;
FEM.fs      = FEM.force;

% set the enssential boundary condition
if sum(options.ess_bc_flag) == 0
    FEM.c = FEM.c_mean;
else
    FEM.c = 0;
    for i = 1:4
        FEM.c = FEM.c + FEM.bnd_mat{i}*options.ess_bc_flag(i)*FEM.penalty;
        FEM.fs(FEM.bnd_ind{i}) = FEM.fs(FEM.bnd_ind{i}) + options.ess_bc_vals(i)*FEM.penalty;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% observation operator
%
% switch options.obs_type
%     case{'Nodes'}
%         N     = size(options.obs_locs,1);
%         FEM.C = [];
%         for k = 1:N
%             ind = mesh.centers(1,:)>(options.obs_locs(k,1)+tol) & mesh.centers(1,:)<(options.obs_locs(k,2)-tol) & ...
%                 mesh.centers(2,:)>(options.obs_locs(k,3)+tol) & mesh.centers(2,:)<(options.obs_locs(k,4)-tol);
%             tmp = 1:mesh.Nel;
%             elems = tmp(ind);
%             C   = sparse([],[],[],1,mesh.Nnode,ceil((sqrt(sum(ind))+1)^2));
%             for i = 1:sum(ind)
%                 nodes = mesh.node_map(:,elems(i));
%                 dx    = mesh.node(:,nodes(3)) - mesh.node(:,nodes(1));
%                 detJ  = prod(abs(dx)); % iJ = diag(1./dx);
%                 locs  = 0.25*detJ*ones(1,4);
%                 C(1,nodes)  = C(1,nodes) + locs;
%             end
%             FEM.C = [FEM.C; C];
%         end
%         FEM.Nsensors = size(FEM.C, 1);
%
%     case{'Interp'}
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

FEM.Ndata       = FEM.Nsensors;
FEM.Ndatasets   = 1;
FEM.Nsteps      = 1;
FEM.c           = sparse(FEM.c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make qoi
FEM.flux_flag = false;
if ~isempty(options.qoi_flux)
    %{
    nq  = length(options.qoi_flux);
    phi = zeros(nq, mesh.Nnode);
    for i = 1:nq
        phi(i,:) = options.qoi_flux{i}(mesh.node);
    end
    %}
    phi = reshape(options.qoi_flux(mesh.node), 1, []);
    phi(abs(phi)<tol) = 0;
    FEM.phi = sparse(phi);
    FEM.flux_flag   = true;
end

end