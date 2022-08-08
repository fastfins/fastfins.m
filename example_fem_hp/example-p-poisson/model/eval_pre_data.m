function [fstate, grad_s] = eval_pre_data(model, state)
%eval_pre_data
%
% pre evaluate functions at quadrature points.
%
% Tiangang Cui, 10/Sep/2020

ne = model.mesh.num_elem;
nd = model.local_elem.num_dim;
nq = model.local_elem.num_quad;

% map grad to quad
grad_s = calc_grad(model.mesh, model.local_elem, state);

fstate = zeros(nq, ne);
for di = 1:nd
    fstate = fstate + grad_s{di}.^2;
end
fstate = fstate + model.epsilon;

end