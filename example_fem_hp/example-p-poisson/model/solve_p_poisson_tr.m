function sol = solve_p_poisson_tr(model, Mbnd)
%solve_p_poisson_tr
%
% Solve the forward model using trust region and Newton iteration
%
% Tiangang Cui, 10/Sep/2020

if model.iter_disp == 0
    dd = 'off';
else
    dd = 'iter';
end
opt  = optimoptions('fminunc', 'Algorithm', 'trust-region', 'SpecifyObjectiveGradient',true,... 
    'HessianMultiplyFcn',@(HI,v) matvec_hessian(HI, v), ...
    'Display',dd, 'MaxIterations', model.max_newton_iter, 'StepTolerance', 1e-20, 'OptimalityTolerance', 1e-20, 'FunctionTolerance', 1e-20);

sol.state = fminunc_2020a(@(s) obj(model, Mbnd, s), model.b, opt);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, g, H] = obj(model, Mbnd, state)

%tic;
[fs, grad_s] = eval_pre_data(model, state);

% energy norm
f = eval_energy(model, Mbnd, fs, state);
%toc

[He, Hb] = assemble_matrices(model, fs, grad_s);

g = He*state + Mbnd*state - model.b;
    
% hessian
H = He + Hb + Mbnd;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
function Hv = matvec_hessian(hess, v)

Hv  = hess*v;

end