function [Ms, Ks] = build_rom_matrices(model, U, Vmu, Vk)
%set up the reduced order model for the 2nd order PDE
%
%inputs:
%  U:   nxr0, the projection basis for the states
%  Vk:  nxr1, the projection basis for the kappa field
%  Vmu: nxr1, the projection basis for the mu_a field

ne = model.mesh.num_elem;
nd = model.local_elem.num_dim;
nb = model.local_elem.num_node;
nq = model.local_elem.num_quad;
nU = size(U,2);

% evaluate grad of u for each i
gradu = cell(nd,1);
for di = 1:nd
    gradu{di} = zeros(nq*ne, nU);
end
u = zeros(nq*ne, nU);
%
for i = 1:nU
    tmp = reshape(U(model.mesh.node_map',i), nb, ne);
    for di = 1:nd
        for dj = 1:nd
            invJtij_duj = model.local_elem.grad_f_quad_pts{dj}*(tmp.*reshape(model.mesh.inv_Jt{di,dj},1,[]));
            gradu{di}(:,i) = gradu{di}(:,i) + reshape(invJtij_duj,[],1);
        end
    end
    u(:,i) = reshape(model.local_elem.f_quad_pts*tmp,[],1);
end
%
WdetJ = model.local_elem.quad_weights(:)*model.mesh.detJ(:)';
%
nmu = size(Vmu,2);
Ms = zeros(nU^2, nmu);
for i = 1:nmu
    mu_atq = model.local_elem.f_quad_pts*reshape(Vmu(model.mesh.node_map',i), nb, ne);
    tmp = (u.*mu_atq(:).*WdetJ(:))'*u;
    Ms(:,i) = tmp(:);
end
%
nK = size(Vk,2);
Ks = zeros(nU^2, nK);
for i = 1:nK
    k_atq = model.local_elem.f_quad_pts*reshape(Vk(model.mesh.node_map',i), nb, ne);
    tmp = zeros(nU);
    for di = 1:nd
        tmp = tmp + (gradu{di}.*k_atq(:).*WdetJ(:))'*(gradu{di});
    end
    Ks(:,i) = tmp(:);
end
%

end