function [model, obs, prior] = set_nonlinear_poisson(m,n)

%m = 20;
%n = 200;
%model.exp_flag is redundant, keep it for consistency with other codes

sd = ((1:m)'.^(-2)).*rand(m,1);

model.F = (rand(m)*(sd.*rand(m,n)));
model.exp_flag = true;
model.type = 'poisson';

D = (1:n).^(-1);
prior.D = 4*D(:);
prior.mean_u = zeros(n,1);
prior.dof = n;
prior.type = 'Field';

obs.true_u = prior.D.*randn(n,1);

%model.Is = 20;
%tmp = model.Is*exp(-model.F*exp(obs.true_u));
tmp = model.F*exp(obs.true_u);

obs.data  = poissrnd(tmp);
obs.n_data = m;
obs.like = 'poisson';
obs.jacobian = false;

end