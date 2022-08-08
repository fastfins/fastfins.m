function FEM = laplace_make_FEM(mesh)
%LAPLACE_MAKE_FEM    
%
% makes the basic  FEM structure
%
% Tiangang Cui, 03/May/2014

FEM.W1              = sparse([],[],[],mesh.Nnode,mesh.Nel,4*mesh.Nel);
FEM.W2              = sparse([],[],[],mesh.Nnode,mesh.Nel,4*mesh.Nel);
FEM.W3              = sparse([],[],[],mesh.Nnode,mesh.Nel,4*mesh.Nel);
FEM.M               = sparse([],[],[],mesh.Nnode,mesh.Nnode,10*mesh.Nel);
FEM.Mb              = sparse([],[],[],mesh.Nnode,mesh.Nnode, 4*mesh.Nbndf);

for i = 1:mesh.Nbndf
    ind             = mesh.node_map_bnd(:,i);
    dx              = mesh.node(:,ind(2)) - mesh.node(:,ind(1));
    FEM.Mb(ind,ind) = FEM.Mb(ind,ind) + norm(dx)*mesh.locmass_bnd;
end

for i = 1:mesh.Nel
    ind             = mesh.node_map(:,i);
    % assume const and identity Jacobian, note here the local stiffness
    % matrix does not contain the Jacobian since square grid is used
    dx              = mesh.node(:,ind(3)) - mesh.node(:,ind(1));
    detJ            = prod(abs(dx)); % iJ = diag(1./dx);
        
    FEM.W1(ind,i)   = mesh.w1;
    FEM.W2(ind,i)   = mesh.w2;
    FEM.W3(ind,i)   = mesh.w3;
   
    FEM.M(ind,ind)  = FEM.M(ind,ind) + mesh.locmass*detJ;    
end

FEM.MR              = chol(FEM.M);

% find bounding box
max_xy  = max(mesh.node, [], 2);
min_xy  = min(mesh.node, [], 2);

tol = 1E-10;
% find nodes on the left boundary 
FEM.bnd_ort     = {'l', 'b', 'r', 't'};
FEM.bnd_mat     = cell(4,1);
FEM.bnd_ind     = cell(4,1);

FEM.bnd_ind{1}  = find(abs(mesh.node(1,:) - min_xy(1)) < tol);
FEM.bnd_ind{2}  = find(abs(mesh.node(2,:) - min_xy(2)) < tol);
FEM.bnd_ind{3}  = find(abs(mesh.node(1,:) - max_xy(1)) < tol);
FEM.bnd_ind{4}  = find(abs(mesh.node(2,:) - max_xy(2)) < tol);

for i = 1:4
    FEM.bnd_mat{i}  = sparse(FEM.bnd_ind{i},FEM.bnd_ind{i},1,mesh.Nnode,mesh.Nnode);
end

% default penalty matrix for integrated value over boundary (the constraint vector)
c = sparse(mesh.node_map_bnd(1,:),1,1/mesh.Nbndf,mesh.Nnode,1) + ...
    sparse(mesh.node_map_bnd(2,:),1,1/mesh.Nbndf,mesh.Nnode,1); % projection is now c*c'
FEM.c_mean      = c*c'*mesh.Nel; 

% Precompute permutation vector for modified stiffness matrix to give 
% minimum bandwidth. 
Apen            = FEM.W1*FEM.W1' + FEM.W2*FEM.W2' + FEM.W3*FEM.W3' + FEM.c_mean;

FEM.p           = symamd(Apen);   % this is a MatLab function
FEM.r(FEM.p)    = 1:size(FEM.p,2); % inverse permutation

FEM.DoF         = mesh.Nnode;

for i = 1:mesh.Nbndf
    ind = mesh.node_map_bnd(:,i);
    dx  = mesh.node(:,ind(2)) - mesh.node(:,ind(1));
    FEM.Mb(ind,ind) = FEM.Mb(ind,ind) + norm(dx)*mesh.locmass_bnd;
end

FEM.penalty     = 1E8; % ( mesh.xyratio/max(mesh.Nside,mesh.Mside) )^(-2);
FEM.mesh        = mesh;

end