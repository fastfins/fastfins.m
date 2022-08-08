function Jty = matvec_Jty(model, sol, dy)
%matvec_Jty
%
% Adjoint solve
%
% Tiangang Cui, 10/Sep/2020

rhs = -model.obs_operator'*dy;

% adjoint state
p = zeros(model.mesh.dof, size(dy,2));

switch model.solver
    case {'trust_region'}
    case {'line_search'}
        p(sol.p,:) = sol.L'\(sol.L\(rhs(sol.p,:)));
end

s_bnd_at_q = model.local_elem_bnd.f_quad_pts*sol.state(model.pa_bnd_mesh.node_map');

% adjoint state at boundary
Jty = zeros(size(model.pa_bnd2domain,2), size(dy,2));
for i = 1:size(dy,2)
    loc_p = reshape(p(model.pa_bnd_mesh.node_map', i), model.local_elem_bnd.num_node, []);
    p_bnd_at_q = model.local_elem_bnd.f_quad_pts*loc_p;
    %
    % calc matvec of adjoint model with dy
    fil = model.local_elem_bnd.f_quad_pts' * (p_bnd_at_q.*s_bnd_at_q.*model.pa_bnd_WdetJ);
    tmp = accumarray(model.pa_bnd_fill_vec(:), fil(:), [model.mesh.dof,1], [], [], false);
    Jty(:,i) = model.pa_bnd2domain' * tmp;
end

Jty = Jty.*sol.dxdu;

end
