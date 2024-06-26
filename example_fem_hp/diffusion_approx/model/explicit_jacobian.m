function J = explicit_jacobian(model, sol)
%
% compute the Jacobian
% Tiangang Cui, 10/Mar/2013

t = model.obs_operator'.*sol.mu_a;
lambda = zeros(size(sol.state,1), model.n_sensors);
lambda(sol.p,:) = -sol.L'\(sol.L\t(sol.p,:));

J = zeros(model.n_sensors*model.n_datasets, length(sol.dxdu(:)));

for k = 1:model.n_datasets
    bi = (k-1)*model.n_sensors;
    for j = 1:model.n_sensors
        J(bi+j,:) = deri_adjoint_stiff_sol(model.mesh, model.local_elem, 'scalar', [], lambda(:,j), sol.state(:,k));
        J(bi+j,:) = J(bi+j,:).*sol.dkdx(:)';
        J(bi+j,:) = J(bi+j,:) + deri_adjoint_mass_sol(model.mesh, model.local_elem, lambda(:,j), sol.state(:,k))';
        J(bi+j,:) = J(bi+j,:) + model.obs_operator(j,:).*sol.state(:,k)';
    end
end

% transform to physical space x
J = J.*sol.dxdu(:)';

end
