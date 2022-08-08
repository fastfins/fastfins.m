function [He, Hb] = assemble_matrices(model, fstate, grad_s)
%assemble_stiff_scalar
%
% assemble the stiffness matrix for a scalar field viscosity
%
% Tiangang Cui, 10/Sep/2020
%

ne = model.mesh.num_elem;
nd = model.local_elem.num_dim;
nb = model.local_elem.num_node;
nq = model.local_elem.num_quad;

%
sqrt_eta = fstate.^(model.p_rate/4-1/2);
grad = cell(nd, 1);
for di = 1:nd
    grad{di} = model.grad{di}.*repmat(sqrt_eta, nb, 1);
end

cross = zeros(nq*nb, ne);
tmpb = fstate.^(model.p_rate/4-1);
for di = 1:nd
    bi = tmpb.*grad_s{di};
    cross = cross + model.grad{di}.*repmat(bi, nb, 1);
end

fill_e = zeros(nb^2, ne);
fill_b = zeros(nb^2, ne);
for ei = 1:ne
    loc = zeros(nb, nb);
    for di = 1:nd
        tmp = reshape(grad{di}(:,ei), nq, nb);
        loc = loc + tmp'*(model.local_elem.quad_weights(:).*tmp);
    end
    fill_e(:,ei) = loc(:)*model.mesh.detJ(ei);
    %
    tmp = reshape(cross(:,ei), nq, nb);
    loc = tmp'*(model.local_elem.quad_weights(:).*tmp);
    fill_b(:,ei) = loc(:)*model.mesh.detJ(ei)*(model.p_rate-2);
end

He = accumarray(model.fill_i, fill_e(:), model.fill_d, [], [], true);
Hb = accumarray(model.fill_i, fill_b(:), model.fill_d, [], [], true);

end
