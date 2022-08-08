function f = eval_energy(model, Mbnd, fstate, state)

% energy norm
% f = 0;
f1 = (1/model.p_rate)*sum( ( fstate(:).^(model.p_rate/2) ).*model.WdetJ(:) );

% forcing
f2 = - model.b'*state;

f3 = 0.5 * state'*(Mbnd*state);

f = f1 + f2 + f3;

%%%%
%{
beta_full = model.pa_bnd2domain*beta;
beta_at_q = model.local_elem_bnd.f_quad_pts*beta_full(model.pa_bnd_mesh.node_map');

% map state to boundary quad, without Jacobian and quadrature scaling
s_bnd_at_q = model.local_elem_bnd.f_quad_pts*state(model.pa_bnd_mesh.node_map');

f11 = (1/model.p_rate)*(model.local_elem.quad_weights(:)' * ( fstate.^(model.p_rate/2) ) * model.mesh.detJ(:) );
disp(f1-f11)

% boundary debug
f31 = sum(beta_at_q(:).*(s_bnd_at_q(:).^2).*model.pa_bnd_WdetJ(:));
disp(f3-0.5*f31);
%}
end