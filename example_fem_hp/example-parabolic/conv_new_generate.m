% conv_new_generate samples

rng(1024);

n_samples = 5E2;

qois = [];
obss = [];
us = [];

tic;
for i = 1:n_samples
    u = prior_random(prior, 1);
    sol = forward_solve(model, obs.true_u);

    us = [us, u(:)];
    qois = [qois, sol.q];
    obss = [obss, sol.d];
end
toc