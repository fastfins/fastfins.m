
%%
i = 2;
model = models{i};
obs =  obss{i};
prior = priors{i};

%%
u = obs.true_u + randn(prior.dof,1)*1E-2;
tol = 1E-5;

Fp = zeros(obs.n_data, prior.dof);
Fm = zeros(obs.n_data, prior.dof);
for i = 1:prior.dof
    up = u;
    up(i) = up(i) + tol;
    um = u;
    um(i) = um(i) - tol;
    solp = forward_solve(model, up);
    Fp(:,i) = solp.d;
    solm = forward_solve(model, um);
    Fm(:,i) = solm.d;
end
Jd = (Fp-Fm)/(2*tol);

%%
sol = forward_solve(model, u);
J = explicit_jacobian(model, sol);
norm(J - Jd, 'fro')

Jt = zeros(size(J'));
for i = 1:obs.n_data
    e = zeros(obs.n_data,1);
    e(i) = 1;
    Jt(:,i) = matvec_Jty(model, sol, e);
end
norm(J - Jt', 'fro')

Je = zeros(size(J));
for i = 1:prior.dof
    e = zeros(prior.dof,1);
    e(i) = 1;
    Je(:,i) = matvec_Ju(model, sol, e);
end
norm(J - Je, 'fro')
