function [As, inds] = build_poisson_rom_matrices(model, grad_s, Veta, Vb, U, Vx)
%set up the reduced order model for the 2nd order PDE
%
%inputs:
%  U:  nxr, the projection basis for the states
%  Vs: mx1, the projection basis for the tensor field
%      scalar: m = 1
%      vector: m = d
%      tensor: m = d+1

ne = model.mesh.num_elem;
nd = model.local_elem.num_dim;
nb = model.local_elem.num_node;
nq = model.local_elem.num_quad;
nU = size(U,2);

As = cell(3,1);
inds = cell(3,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for the diagonal terms
nK = size(Veta,2);
% assemble squared K basis
nKs = (nK+1)*nK/2;
As{1} = zeros(nU^2, nKs);
inds{1} = zeros(nKs, 1);
k = 0;
for i = 1:nK
    for j = i:nK
        tmp = zeros(nU);
        for di = 1:nd
            if i == j
                tmp = tmp + (grad_s{di}.*Veta(:,i).*model.WdetJ(:))'*(grad_s{di}.*Veta(:,j));
            else
                tmp = tmp + (grad_s{di}.*Veta(:,i).*model.WdetJ(:))'*(grad_s{di}.*Veta(:,j)) + ...
                    + (grad_s{di}.*Veta(:,j).*model.WdetJ(:))'*(grad_s{di}.*Veta(:,i));
            end
        end
        k = k + 1;
        As{1}(:,k) = tmp(:);
        inds{1}(k) = (j-1)*nK + i;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nB = size(Vb,2);
%
nBs = (nB+1)*nB/2;
As{2} = zeros(nU^2, nBs);
inds{2} = zeros(nBs, 1);
k = 0;
for i = 1:nB
    for j = i:nB
        left  = zeros(size(grad_s{di}));
        right = zeros(size(grad_s{di}));
        for di = 1:nd
            ind   = (1:nq*ne) + (di-1)*nq*ne;
            left  = left  + grad_s{di}.*Vb(ind,i);
            right = right + grad_s{di}.*Vb(ind,j);
        end
        if i == j
            tmp = (left.*model.WdetJ(:))'*right;
        else
            tmp = (left.*model.WdetJ(:))'*right + (right.*model.WdetJ(:))'*left;
        end
        k = k + 1;
        As{2}(:,k) = tmp(:);
        inds{2}(k) = (j-1)*nB + i;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% map states to bnd
ubnd = zeros(model.local_elem_bnd.num_quad*model.pa_bnd_mesh.num_elem, nU);
for i = 1:nU
    ui = U(:,i);
    ubnd(:,i) = reshape(model.local_elem_bnd.f_quad_pts*ui(model.pa_bnd_mesh.node_map'), [], 1);
end
%
nx = size(Vx,2);
beta = zeros(model.local_elem_bnd.num_quad*model.pa_bnd_mesh.num_elem, nx);
for i = i:nx
    beta_full = model.pa_bnd2domain*Vx(:,i);
    beta(:,i) = reshape(model.local_elem_bnd.f_quad_pts*beta_full(model.pa_bnd_mesh.node_map'), [], 1);
end
%
nxs = (nx+1)*nx/2;
As{3} = zeros(nU^2, nxs);
inds{3} = zeros(nxs, 1);
k = 0;
for i = 1:nx
    for j = i:nx
        if i == j
            tmp = (ubnd.*beta(:,i).*model.pa_bnd_WdetJ(:))'*(ubnd.*beta(:,j));
        else
            tmp = (ubnd.*beta(:,i).*model.pa_bnd_WdetJ(:))'*(ubnd.*beta(:,j)) + ...
                + (ubnd.*beta(:,j).*model.pa_bnd_WdetJ(:))'*(ubnd.*beta(:,i));
        end
        k = k + 1;
        As{3}(:,k) = tmp(:);
        inds{3}(k) = (j-1)*nx + i;
    end
end

end