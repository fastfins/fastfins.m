function prior = make_prior_mrf(mesh, k, cond, sigma)
%MAKE_PRIOR_MRF    
%
% makes the GMRF prior distribution for a rectangular grid 
%
% Q u = w, Q = inv(sigma) * (K(cond) + k I)
%
%%%%%%%%%%%%%%%%%%%% input: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mesh:     
% k:            variable controls the eigen spectrum of the prior
% sigma:        the standard deviation
% cond:         
%
% Tiangang Cui, 09/May/2014

mesh_matrix  = make_matrix(mesh);
prior = update_mrf(mesh_matrix, k, cond, 1/sigma);
prior.mesh_matrix = mesh_matrix;
prior.mesh = mesh;

tmp = prior.Misq*prior.K*prior.Misq;
tmp = (tmp + tmp')*0.5;
prior.chiAs = eig(tmp,'vector');

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = make_matrix(mesh)

out.ind_i   = zeros(4^2, mesh.Nel);
out.ind_j   = zeros(4^2, mesh.Nel);
out.Nnode  = mesh.Nnode;
out.Nel    = mesh.Nel;

detJs = zeros(1, mesh.Nel);

for i = 1:mesh.Nel
    ind     = mesh.node_map(:,i);
    
    % assume const and identity Jacobian
    dx      = mesh.node(:,ind(3)) - mesh.node(:,ind(1));
    detJ    = prod(abs(dx));
    
    ii      = repmat(ind', 4, 1);
    jj      = repmat(ind , 1, 4);
    
    out.ind_i(:,i)  = ii(:);
    out.ind_j(:,i)  = jj(:);
    detJs(i)        = detJ;
end

out.loc_xx  = mesh.loc.xx(:);
out.loc_yy  = mesh.loc.yy(:);
out.loc_xy  = mesh.loc.xy(:);

%temp_K = mesh.loc.xx(:)*cond(:,1)' + mesh.loc.yy(:)*cond(:,2)' + mesh.loc.xy(:)*cond(:,3)';
%K = sparse(ind_i(:), ind_j(:), temp_K(:), mesh.N_node, mesh.N_node);

temp_M      = mesh.locmass(:)*detJs;
M           = sparse(out.ind_i(:), out.ind_j(:), temp_M(:), mesh.Nnode, mesh.Nnode);
out.ML      = full(sum(M,2));

end