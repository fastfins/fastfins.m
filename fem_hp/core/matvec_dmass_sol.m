function g = matvec_dmass_sol(mesh, local_elem, df, state)

fill_i  = mesh.node_map';
%
loc_df = local_elem.f_quad_pts*df(mesh.node_map');
loc_state = local_elem.f_quad_pts*state(mesh.node_map');
fill  = (local_elem.f_quad_pts.*local_elem.quad_weights(:))'*(loc_df.*loc_state.*mesh.detJ(:)');

g = accumarray(fill_i(:), fill(:), [mesh.dof,1], [], [], false);

end
