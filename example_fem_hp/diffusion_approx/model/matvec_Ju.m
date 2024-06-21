function Ju = matvec_Ju(model, sol, du)

% transform to physical space x
dx = sol.dxdu.*du;
dk = sol.dkdx.*dx;

Ju = zeros(model.n_sensors, model.n_datasets*size(dx,2));
dstate = zeros(size(sol.state));

for i = 1:size(dx,2)
    t1 = matvec_dmass_sol(model.mesh, model.local_elem, dx, sol.state);
    t2 = matvec_dstiff_sol(model.mesh, model.local_elem, 'scalar', [], dk, sol.state);
    
    dstate(sol.p,:) = -sol.L'\(sol.L\(t1(sol.p,:)+t2(sol.p,:)));
    ind = (1:model.n_datasets) + (i-1)*model.n_datasets;
    Ju(:,ind) = model.obs_operator*dstate;
end

Ju = reshape(Ju, [], size(dx,2));

end