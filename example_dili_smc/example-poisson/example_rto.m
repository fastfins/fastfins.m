% This script demonstrates running the RTO algorithm


% build the linearisation
model.explicit_ja = false;
rto_lin = rto_local_lin(model, obs, prior, vmap, 1E-2, obs.n_data-1);


R = randn(prior.dof, 100);
tic;
rto_out0 = rto_samples(model, obs, prior, rto_lin, R);
toc

% markov chain on the indices
[ind_chain0, acc0] = rto_metropolize(rto_out0.logrhos);
samples0 = rto_out0.us(:,ind_chain0);

% Notes:
%
% following the derivation of the subspace RTO, the RTO density is given by 
%   exp( - 0.5 | Q^T H(x) |^2) det( Q^T H(x) )
% In the code, det( Q^T H(x) ) is given by out.logds + out.logc, there is a
% bit complications in computing the term ( - 0.5 | Q^T H(x) |^2 ), it is 
% split into ( -HI.mlrto - 0.5*|v_perp|^2 ), here we only compute HI.mlrto
% in the code for the nonlinear transformation in the split RTO formula
%
% for the posterior, we have exp( -0.5 |F(x) - y|^2 - 0.5 |v|^2 ), where
% the term ( - 0.5 |v|^2 ) split into ( 0.5*|v_r|^2 - 0.5*|v_perp|^2 ), so
% in the ratio posterior/rto, we only need to compute
%   exp( -0.5 |F(x) - y|^2 - 0.5 |v_r|^2 + HI.mlrto - out.logds - out.logc)
% or 
%   exp( -0.5 |F(x) - y|^2 - 0.5 |v_r|^2 + out.logrtos ). 
%
% In cases the RTO density is defined using a model other thann the forward
% model used in the posterior, we do need to compute the terms
%   exp( -0.5 |F(x) - y|^2 - 0.5 |v_r|^2 )
% this is given by 
%   out.vrs  --- v_r parameter in the subspace
%   out.vs   --- full space parameter with a i.i.d. Gaussian prior
%   out.us   --- full space parameter with a correlated Gaussian prior
%
% us = sqrt(C)*vs 
%
%

model.explicit_ja = true;
tic;
rto_out1 = rto_samples(model, obs, prior, rto_lin, R);
toc
% markov chain on the indices
[ind_chain1, acc1] = rto_metropolize(rto_out1.logrhos);
samples1 = rto_out1.us(:,ind_chain1);

% trust region version, 
% modified forward model is 
%   t( ( F(x) - J(x-x_ref) - F(x_ref) )./d ; tau, e ) + J(x-x_ref) + F(x_ref)
% where d is a scaling reference d = |F(x_ref)| + 1
% the trust region function is defined as
%   t(x; tau, e) = x       for abs(x) < tau(1-e)
%   t(x; tau, e) = sign(x) for abs(x) > tau(1+e)
% 
% in the following setup, tau = 1 and e = 0.2
%
% this basically turns off the nonlinear forward model after its
% sufficiently far away from the linearisation relative to the reference 
% model output. This effectively make the modified forward model linear at 
% infinity. 
%
tic;
rto_out2 = rto_samples_trust(model, obs, prior, rto_lin, R, 1, 0.2);
toc

% markov chain on the indices
[ind_chain2, acc2] = rto_metropolize(rto_out2.logrhos);
samples2 = rto_out2.us(:,ind_chain2);
