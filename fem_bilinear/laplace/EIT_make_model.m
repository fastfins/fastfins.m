function FEM = EIT_make_model(mesh, options)
%EIT_MAKE_OUTPUTS   
%
% generates measurements points and indices w.r.t. the FEM mesh
%
% Tiangang Cui, 03/May/2012

tol = 1E-10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FEM = laplace_make_FEM(mesh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FEM.sensors     = zeros(4*length(options.obs_locs),1);
count           = 0;
% y = 0
for i   = 1:length(options.obs_locs)
    j       = find(abs(mesh.node(2,:))<tol & abs(mesh.node(1,:)-options.obs_locs(i))<tol);
    count   = count + 1;
    FEM.sensors(count)  = j;
end

% x = xyratio
for i   = 1:length(options.obs_locs)
    j       = find(abs(mesh.node(1,:)-mesh.xyratio)<tol & abs(mesh.node(2,:)-options.obs_locs(i))<tol);
    count   = count + 1;
    FEM.sensors(count)  = j;
end

% y = 1
for i   = length(options.obs_locs):-1:1
    j       = find(abs(mesh.node(2,:)-1)<tol & abs(mesh.node(1,:)-options.obs_locs(i))<tol);
    count   = count + 1;
    FEM.sensors(count)  = j;
end

% x = 0
for i   = length(options.obs_locs):-1:1
    j       = find(abs(mesh.node(1,:))<tol & abs(mesh.node(2,:)-options.obs_locs(i))<tol);
    count   = count + 1;
    FEM.sensors(count)  = j;
end


FEM.Nsensors    = 4*length(options.obs_locs);
FEM.fs  = sparse([],[],[], mesh.Nnode, FEM.Nsensors, mesh.Nnode*FEM.Nsensors); 

for i = 1:FEM.Nsensors
    % implement boundary reference, i.e. current removed evenly round
    % boundary
    % FEM.fs(mesh.node_map_bnd(1,:),count) = -(1/mesh.N_bnd_f)*ones(mesh.N_bnd_f,1); 
    % fsp(FEM.sensors(count),count)        = FEM.fs(FEM.sensors(count),count) + 1;
    
    fsm = zeros(mesh.Nnode, 1);
    
    fsm(mesh.node_map_bnd(1,:)) = -(1/4)*ones(mesh.Nbndf,1); 
    fsm = FEM.Mb*fsm;
    
    xn  = mesh.node(:,FEM.sensors(i));
    dis = sum( (mesh.node - repmat(xn, 1, mesh.Nnode)).^2 )';
    fsp = FEM.M*exp(-0.5*dis/(2*0.02^2))/(pi*0.02^2);
    
    %sum(fsp)
    %sum(fsm)
    
    FEM.fs(:,i)  = fsp + fsm;
end


FEM.Ndatasets   = FEM.Nsensors;
FEM.Ndata       = FEM.Nsensors^2;
FEM.Nsteps      = 1;

FEM.C           = sparse([],[],[], FEM.Nsensors, mesh.Nnode, FEM.Nsensors); % make the observation matrix
FEM.C(:,FEM.sensors)= speye(FEM.Nsensors);

FEM.c = FEM.c_mean;
FEM.flux_flag = false;
end
