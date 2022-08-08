function [model, obs, prior] = set_nonlinear_example(m,n)

%m = 20;
%n = 200;
%model.exp_flag is redundant, keep it for consistency with other codes

[Ql,~] = qr(randn(m),0);
[Qr,~] = qr(randn(n,m),0);
%sd = 10.^linspace(2, 0, m);
sd = 1E2*(1:m).^(-1);
%sd = 10.^linspace(2, 0, m)/100;

model.F = Ql*(sd(:).*Qr');
model.exp_flag = true;
model.qoi = 1E-4;

D = (1:n).^(-1);
prior.D = 2*D(:);
prior.mean_u = zeros(n,1);
prior.dof = n;
prior.type = 'Field';

obs.true_u = prior.D.*randn(n,1);

if model.exp_flag
    obs.data = model.F*exp(obs.true_u) + randn(m,1);
    model.type = 'exp';
else
    obs.data = model.F*obs.true_u + randn(m,1);
    model.type = 'linear';
end
obs.n_data = m;
obs.std = 1;
obs.like = 'normal';
obs.jacobian = false;

end