function [U, s, V, sol] = svd_rand_WJ(model, obs, prior, sol, tol, nmax)
%SVD_WsolTENING_F_RAND
%
% randomized SVD for factorizing the forward model
% requires the sol by runing 
% [~,~,~, sol] = minus_log_post(model, prior, v);
%
% Tiangang Cui, 17/Jan/2014

if ~isstruct(sol)
    [~,~,~,~,sol] = minus_log_post(model, obs, prior, sol);
end

if prior.dof >= obs.n_data
    [U, s, V] = more_param(model, obs, prior, sol, tol, nmax);
else
    [U, s, V] = more_obs  (model, obs, prior, sol, tol, nmax);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U, s, V] = more_param(model, obs, prior, sol, tol, nmax)

% Let O*O = Gamma_{obs}^{-0.5}, L*L' = \Gamma_{pr}
% J  = O *F*diag(dxdu) *L
% Jt = L'*diag(dxdu)*F'*O'

 Nrand      = min(nmax, obs.n_data+5);
 
% apply J' to data space perturbation
 Juty       = zeros(prior.mesh_dof, Nrand);
 O          = randn(obs.n_data, Nrand);
 for i = 1:Nrand
     Juty(:,i)  = matvec_Jty(model, sol, O(:,i)./obs.std);
 end
 Jvty       = matvec_prior_Lt(prior, Juty);

% apply J to param space perturbation
[Qv, ~]     = qr(Jvty,0);
 Qu         = matvec_prior_L (prior, Qv);
 Jvx        = zeros(obs.n_data, Nrand);
for i = 1:Nrand
    Jvx(:,i)    = matvec_Ju(model, sol, Qu(:,i))./obs.std;
end

% Jvx = J*Qv = A * T * B'
% J   = J*Qv*Qv'
% J   = A * T * B'*Qv'
 
[A, T, B]  = svd(Jvx); % final svd
 
 dT         = diag(T);
 ind        = find(dT>=tol); % truncate
 U          = A(:,ind);
 V          = Qv*B(:,ind);
 s          = dT(ind);
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U, s, V] = more_obs(model, obs, prior, sol, tol, nmax)

% Let O*O = Gamma_{obs}^{-0.5}, L*L' = \Gamma_{pr}
% J  = O *F*diag(dxdu) *L
% Jt = L'*diag(dxdu)*F'*O'

 Nrand      = min(nmax, prior.dof+5);
 
% apply J to param space perturbation
 O          = randn(prior.dof, Nrand);
 Qu         = matvec_prior_L (prior, O);
 Jvx        = zeros(obs.n_data, Nrand);
for i = 1:Nrand
    Jvx(:,i)    = matvec_Ju(model, sol, Qu(:,i))./obs.std;
end
 
% apply J' to data space perturbation
[Qv, ~]     = qr(Jvx,0);
 Juty       = zeros(prior.mesh_dof, Nrand);
 for i = 1:Nrand
     Juty(:,i)  = matvec_Jty(model, sol, Qv(:,i)./obs.std);
 end
 Jvty       = matvec_prior_Lt(prior, Juty);

% Jvty  = J'*Qv = A * T * B'
% J'    = J'*Qv*Qv'
% J'    = A * T * B'*Qv' 
%       = V * T * U'

[A, T, B]  = svd(Jvty); % final svd
 
 dT         = diag(T);
 ind        = find(dT>=tol); % truncate
 V          = A(:,ind);
 U          = Qv*B(:,ind);
 s          = dT(ind);
 
end