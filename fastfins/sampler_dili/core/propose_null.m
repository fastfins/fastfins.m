function [alpha, next] = propose_null(proposal, ml_target, P, dt, curr, r_null)
% null space proposal for the Hessian separated MCMC
% Tiangang Cui, 25/Mar/2013
%
beta   = 2*sqrt(2*dt)/(2 + dt);
beta_p = (2 - dt)/(2 + dt);
% copy the y_sub
next.v_sub  = curr.v_sub;
% next y in null space, without ref point
next.v_null = beta_p*curr.v_null + beta*r_null;
%next.v_null = beta_p*(curr.v_null-kernel.ref_null) + beta*r_null + kernel.ref_null;
% next y
next.v      = P*next.v_sub + next.v_null;
% eval density
switch proposal
    case {'MALA'}
        [next.mlpt, next.mllkd, next.grad] = ml_target(next.v);
        % project the gradient
        next.grad_v_sub = P'*next.grad;
    otherwise
        [next.mlpt, next.mllkd] = ml_target(next.v);
end
% acceptance
alpha = curr.mllkd - next.mllkd;

end
