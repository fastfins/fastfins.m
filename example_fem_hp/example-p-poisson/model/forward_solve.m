function sol = forward_solve(model, u)
%forward_solve
%
% Mock forward solve
%
% Tiangang Cui, 10/Sep/2020

if model.exp_param
    x = model.beta_weight.*exp(u);
    dxdu = x;
else
    x = u;
    dxdu = ones(size(u));
end

% mass marix for the unknown boundary condition
beta_full = model.pa_bnd2domain*x;
beta = model.local_elem_bnd.f_quad_pts*beta_full(model.pa_bnd_mesh.node_map');
fill = model.local_elem_bnd.mass_at_q*(beta.*model.pa_bnd_WdetJ);
Mbnd = accumarray(model.pa_bnd_fill_i, fill(:), model.fill_d, [], [], true);

switch model.solver
    case {'trust_region'}
        sol = solve_p_poisson_tr(model, Mbnd);
    case {'line_search'}
        sol = solve_p_poisson_ls(model, Mbnd);
end

sol.d = model.obs_operator*sol.state;
sol.dxdu = dxdu;
sol.qoi  = model.qoi_vec*sol.state;

end