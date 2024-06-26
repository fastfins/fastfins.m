function Jty = matvec_Jty(model, sol, dy)

dy = reshape(dy, model.n_sensors, []);
t1 = model.obs_operator'*dy;  % for multiple RHS, size(dy, 2) = N * N_RHS
t2 = sol.mu_a.*t1;
lambda = zeros(size(sol.state,1), size(dy,2));
lambda(sol.p,:) = -sol.L'\(sol.L\t2(sol.p,:));

N   = size(dy,2)/model.n_datasets;
Jty = zeros(length(sol.dxdu(:)), N);

for k = 1:N
    bi = (k-1)*model.n_datasets;
    for j = 1:model.n_datasets
        g1 = deri_adjoint_mass_sol(model.mesh, model.local_elem, lambda(:,bi+j), sol.state(:,j));
        g2 = deri_adjoint_stiff_sol(model.mesh, model.local_elem, 'scalar', [], lambda(:,bi+j), sol.state(:,j));
        Jty(:,k) = Jty(:,k) + g1 + g2.*sol.dkdx + t1.*sol.state;
    end
end

% transform to physical space x
Jty = sol.dxdu.*Jty;

end

