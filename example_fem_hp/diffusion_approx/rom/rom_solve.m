function sol = rom_solve(rom, vr)

r_mu_a = rom.K_mu*exp(rom.redu_u_mu*vr + rom.redu_mean_mu);
r_k = rom.K_k*(1./(2*(exp(rom.redu_u_k*vr+rom.redu_mean_k)+0.2*rom.mu_s)));

Ak = reshape(rom.Ms*r_mu_a + rom.Ks*r_k, rom.dof, rom.dof);
sol.state = (rom.Abnd + Ak)\rom.b;

% apply observation operator, which contains the mass matrix if needed
sol.d  = rom.obs_operator*(exp(rom.param_redu*vr+rom.mean_u).*(rom.states*sol.state));

end

