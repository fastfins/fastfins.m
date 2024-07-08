function sol = rom_solve_debug(model, prior_redu, rom, ur)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sol.mu_a  = exp(prior_redu.basis*ur+prior_redu.mean_u);
sol.kappa = 1./(2*(sol.mu_a + 0.2*model.mu_s)); 

% assemble mass matrix, interior
M = assemble_mass(model.mesh, model.local_elem,  model.fill_i, sol.mu_a, model.fill_d);
% assemble stiffness matrix, interior
K  = assemble_stiff(model.mesh, model.local_elem, model.grad, model.fill_i, sol.kappa, model.fill_d);
% assemble boundary condition
% A  = M + K + model.M_bnd*(2/pi);
% assemble forcing term, do nothing

Mr = rom.states'*(M*rom.states);
Kr = rom.states'*(K*rom.states);
Mbndr = rom.states'*(model.M_bnd*rom.states)*(2/pi);

Ar = Mr + Kr + Mbndr;

%
sol.state = Ar\(rom.states'*model.b);
% apply observation operator, which contains the mass matrix if needed
sol.d  = model.obs_operator*(sol.mu_a.*(rom.states*sol.state));

%
r_mu_a = rom.K_mu*exp(rom.redu_u_mu*ur+rom.redu_mean_mu);
r_k = rom.K_k*(1./(2*(exp(rom.redu_u_k*ur+rom.redu_mean_k)+0.2*rom.mu_s)));

Mr1 = reshape(rom.Ms*r_mu_a, rom.dof, rom.dof);
Kr1 = reshape(rom.Ks*r_k, rom.dof, rom.dof);

Ak = reshape(rom.Ms*r_mu_a + rom.Ks*r_k, rom.dof, rom.dof);
%norm(rom.Abnd + Ak - Ar)
s = (rom.Abnd+Kr1+Mr)\rom.b;

norm(Mr1 - Mr, 'fro')/norm(Ar)
%norm(Kr1 - Kr, 'fro')/norm(Ar)

%norm(s - sol.state)/norm(sol.state)
sol.d  = model.obs_operator*(sol.mu_a.*(rom.states*s));

end


