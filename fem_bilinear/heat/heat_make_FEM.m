function FEM = heat_make_FEM(mesh)
%HEAT_MAKE_FEM  
%
% makes the basic FEM structure 
%
% Tiangang Cui, 31/Oct/2012


FEM.W1  = sparse([],[],[],mesh.Nnode,mesh.Nel,4*mesh.Nel);
FEM.W2  = sparse([],[],[],mesh.Nnode,mesh.Nel,4*mesh.Nel);
FEM.W3  = sparse([],[],[],mesh.Nnode,mesh.Nel,4*mesh.Nel);

FEM.V1  = sparse([],[],[],mesh.Nnode,mesh.Nel,4*mesh.Nel);
FEM.V2  = sparse([],[],[],mesh.Nnode,mesh.Nel,4*mesh.Nel);
FEM.V3  = sparse([],[],[],mesh.Nnode,mesh.Nel,4*mesh.Nel);
FEM.V4  = sparse([],[],[],mesh.Nnode,mesh.Nel,4*mesh.Nel);

FEM.M   = sparse([],[],[],mesh.Nnode,mesh.Nnode,10*mesh.Nel);
FEM.Mb  = sparse([],[],[], mesh.Nnode, mesh.Nnode, 4*mesh.Nbndf);

for i = 1:mesh.Nbndf
    ind = mesh.node_map_bnd(:,i);
    dx  = mesh.node(:,ind(2)) - mesh.node(:,ind(1));
    FEM.Mb(ind,ind) = FEM.Mb(ind,ind) + norm(dx)*mesh.locmass_bnd;
end

for i = 1:mesh.Nel
    ind = mesh.node_map(:,i);
    dx  = mesh.node(:,ind(3)) - mesh.node(:,ind(1));
    detJ  = prod(abs(dx));    % iJ = diag(1./dx);
        
    FEM.W1(ind,i)   = mesh.w1;
    FEM.W2(ind,i)   = mesh.w2;
    FEM.W3(ind,i)   = mesh.w3;
    
    FEM.M(ind,ind)  = FEM.M(ind,ind) + mesh.locmass*detJ;  
    
    FEM.V1(ind,i)   = mesh.v1*sqrt(detJ);
    FEM.V2(ind,i)   = mesh.v2*sqrt(detJ);
    FEM.V3(ind,i)   = mesh.v3*sqrt(detJ);
    FEM.V4(ind,i)   = mesh.v4*sqrt(detJ);
end

FEM.const_detJ  = detJ;


FEM.DoF     = mesh.Nnode;

FEM.K       = FEM.W1*FEM.W1' + FEM.W2*FEM.W2' + FEM.W3*FEM.W3';
FEM.I       = speye(FEM.DoF);
FEM.ML      = spdiags(sum(FEM.M, 2), 0, FEM.DoF, FEM.DoF);
FEM.iML     = spdiags(1./sum(FEM.M, 2), 0, FEM.DoF, FEM.DoF);
FEM.iMLK    = FEM.iML*FEM.K;


FEM.Apen        = FEM.W1*FEM.W1' + FEM.W2*FEM.W2' + FEM.W3*FEM.W3' + FEM.M;

FEM.p           = symamd(FEM.Apen); % this is a MatLab function
FEM.r(FEM.p)    = 1:size(FEM.p,2); % inverse permutation

tol = 1E-10;
% find nodes on the left boundary 
FEM.bnd_ort     = {'l', 'b', 'r', 't'};
FEM.bnd_mat     = cell(4,1);
FEM.bnd_ind     = cell(4,1);

max_xy  = max(mesh.node, [], 2);
min_xy  = min(mesh.node, [], 2);

FEM.bnd_ind{1}  = find(abs(mesh.node(1,:) - min_xy(1)) < tol);
FEM.bnd_ind{2}  = find(abs(mesh.node(2,:) - min_xy(2)) < tol);
FEM.bnd_ind{3}  = find(abs(mesh.node(1,:) - max_xy(1)) < tol);
FEM.bnd_ind{4}  = find(abs(mesh.node(2,:) - max_xy(2)) < tol);

for i = 1:4
    FEM.bnd_mat{i}  = sparse(FEM.bnd_ind{i},FEM.bnd_ind{i},1,mesh.Nnode,mesh.Nnode);
end

FEM.penalty     = ( mesh.xyratio/max(mesh.Nside,mesh.Mside) )^(-2);
FEM.mesh        = mesh;

end

