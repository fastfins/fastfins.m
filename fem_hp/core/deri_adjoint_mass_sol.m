function g = deri_adjoint_mass_sol(mesh, local_elem, lambda, state)

fill_i  = mesh.node_map';

loc_lambda = local_elem.f_quad_pts*lambda(mesh.node_map');
loc_state = local_elem.f_quad_pts*state(mesh.node_map');
%
fill = ( local_elem.f_quad_pts.*local_elem.quad_weights(:) )'...
    *( (loc_lambda.*loc_state).*mesh.detJ(:)' );

g = accumarray(fill_i(:), fill(:), [mesh.dof,1], [], [], false);

g = g(:);

end
