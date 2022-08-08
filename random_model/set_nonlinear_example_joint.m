function [model, obs, prior] = set_nonlinear_example_joint(m,n)

%m = 20;
%n = 200;
%model.exp_flag is redundant, keep it for consistency with other codes

[Ql,~] = qr(randn(m),0);
[Qr,~] = qr(randn(n,m),0);
sd = 10.^linspace(1, -5, m);

model.F = Ql*(sd(:).*Qr');
model.exp_flag = true;

D = (1:n).^(-4);
prior.D = D(:);
prior.mean_u = zeros(n,1);
prior.dof = n;
prior.type = 'Field';

obs.true_u = prior.D.*randn(n,1);

obs.std = 1e-1;
if model.exp_flag
    obs.data = model.F*exp(obs.true_u) + obs.std*randn(m,1);
    model.type = 'exp';
else
    obs.data = model.F*obs.true_u + obs.std*randn(m,1);
    model.type = 'linear';
end
obs.n_data = m;
obs.like = 'normal';
obs.jacobian = false;

end