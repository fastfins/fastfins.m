function [model, obs, prior] = set_nonlinear_2exp(m,n)

%m = 20;
%n = 200;
%model.exp_flag is redundant, keep it for consistency with other codes

[Ql,~] = qr(randn(m),0);
[Qr,~] = qr(randn(n,m),0);
sd = 10.^linspace(2, 0, m)/500;

model.F = Ql*(sd(:).*Qr');
model.exp_flag = true;

D = (1:n).^(-1);
prior.D = 2*D(:);
prior.mean_u = zeros(n,1);
prior.dof = n;
prior.type = 'Field';

obs.true_u = prior.D.*randn(n,1);

tmp = exp(model.F*exp(obs.true_u));

obs.std = mean(tmp)/10;
obs.data = tmp + randn(m,1)*obs.std;
model.type = '2exp';
obs.n_data = m;
obs.like = 'normal';
obs.jacobian = false;

end