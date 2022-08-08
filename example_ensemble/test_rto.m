
beta = 1; % temperature
ind = 1:obs.n_data; % generate all data
y = data_random(model, obs, prior, v, ind, beta);