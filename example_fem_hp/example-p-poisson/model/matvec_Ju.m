function Ju = matvec_Ju(model, sol, du)
%matvec_Ju
%
% Linearised forward solve
%
% Tiangang Cui, 10/Sep/2020

rhs     = zeros(model.mesh.dof, size(du,2));
dstate  = zeros(model.mesh.dof, size(du,2));
dbeta   = sol.dxdu.*du;
% map boundary condition to boundary quad
for i = 1:size(dbeta,2)
    dbeta_full  = model.pa_bnd2domain*dbeta(:,i);
    dbeta_at_q  = model.local_elem_bnd.f_quad_pts*dbeta_full(model.pa_bnd_mesh.node_map');
    %
    fill = model.local_elem_bnd.mass_at_q*(dbeta_at_q.*model.pa_bnd_WdetJ);
    Mbnd = accumarray(model.pa_bnd_fill_i, fill(:), model.fill_d, [], [], true);
    %
    %mass_botbnd = assemble_mass_botbnd(model, dbeta_at_q);
    rhs(:,i) = - Mbnd*sol.state;
end

switch model.solver
    case {'trust_region'}
    case {'line_search'}
        dstate(sol.p,:) = sol.L'\(sol.L\(rhs(sol.p,:)));
end
Ju = model.obs_operator*dstate;

end