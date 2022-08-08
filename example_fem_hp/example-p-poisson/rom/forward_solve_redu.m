function sol = forward_solve_redu(model, states, u)
%forward_solve
%
% Mock forward solve
%
% Tiangang Cui, 10/Sep/2020

if model.exp_param
    x = model.beta_weight.*exp(u);
    dxdu = x;
else
    x = u;
    dxdu = ones(size(u));
end

% mass marix for the unknown boundary condition
beta_full = model.pa_bnd2domain*x;
beta = model.local_elem_bnd.f_quad_pts*beta_full(model.pa_bnd_mesh.node_map');
fill = model.local_elem_bnd.mass_at_q*(beta.*model.pa_bnd_WdetJ);
Mbnd = accumarray(model.pa_bnd_fill_i, fill(:), model.fill_d, [], [], true);


sol.state = solve_p_poisson_redu(model, Mbnd, states);
sol.d = model.obs_operator*sol.state;
sol.dxdu = dxdu;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function state = solve_p_poisson_redu(model, Mbnd, U)
%solve_p_poisson_ls
%
% Solve the forward model using line search and Newton iteration
%
% Tiangang Cui, 29/Oct/2016

state = model.b;
[fs, grad_s] = eval_pre_data(model, state);
cost = eval_energy(model, Mbnd, fs, state);

num_iter = 0;

if model.iter_disp >= 1
    fprintf('Iter    Cost----------------    L2(res)--------------    LS Iter\n');
end

while true  
    num_iter = num_iter + 1;
    % stiffness matrix with eta
    [He, Hb] = assemble_matrices(model, fs, grad_s);
    
    % residual
    res = He*state + Mbnd*state - model.b;
    res_l2 = sqrt(res'*model.mass*res);
    
    % residual condition check
    if res_l2 < model.res_tol 
        % print line search info, residual condition met
        % ls_num_iter = 0
        if model.iter_disp >= 1
            fprintf('%4d    %020.13E    %020.13E    %4d\n', num_iter, cost, res_l2, 0);
        end
        break;
    end
    
    % hessian: H = He + Hb + Mbnd;
    
    % line search direction
    %{
    [L,~,p] = chol(He + Hb + Mbnd, 'lower', 'vector');
    % solve
    update  = zeros(size(res));
    update(p,:) = - L'\(L\res(p,:));
    %}  
    
    A = U'*(He + Hb + Mbnd)*U;
    update = -U*(A\(U'*res));
    
    
    % line search
    [state, cost, ls_num_iter, ds] = armijo_search(model, Mbnd, state, cost, res, update);
        
    % print line search info
    if ds < model.res_tol
        if model.iter_disp >= 1
            fprintf('%4d    %020.13E    %020.13E    %4d\n', num_iter, cost, res_l2, ls_num_iter);
        end
        break;
    end
    
    % iteration check
    if num_iter >= model.max_newton_iter
        % print line search info, max iter reached
        if model.iter_disp >= 1
            fprintf('%4d    %020.13E    %020.13E    %4d\n', num_iter, cost, res_l2, ls_num_iter);
        end
        disp('>>>> Warning: Max Newton iterations reached');
        break;
    end
    
    if model.iter_disp >= 2
        fprintf('%4d    %020.13E    %020.13E    %4d\n', num_iter, cost, res_l2, ls_num_iter);
    end
    
    [fs, grad_s] = eval_pre_data(model, state);
end

state;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [state_out, cost_out, num_iter, ds] = armijo_search(model, Mbnd, state, cost, grad, update)
% ARMIJO_SEARCH    performs a line search satisfies the armijo condition
% [next, a] = armijo_search(opt_def, curr)
%
% given a direction of the line search, compute the next point that satify
% the Armijo condition
%
%   --  state:      current point
%   --  cost:       function value
%   --  grad:       gradient/residual
%   --  update:     search direction
%
% find cost(state + alpha*update) <= cost(state) + ftol*alpha*grad'*update

% directional derivative
dd = grad'*update;
% iteration
num_iter = 0;


if dd > 0
    disp('>>>> Error: Bad search direction');
    state_out = state;
    cost_out = cost;
    return;
end


alpha = 1;
c = 1E-4;
golden_ratio = 0.5*(sqrt(5) + 1);

while 1
    num_iter = num_iter + 1;
    state_new = state + alpha*update;
    [fs,~] = eval_pre_data(model, state_new);
    cost_new = eval_energy(model, Mbnd, fs, state_new);
    % check for suffcient decrease
    if cost_new - cost < c*alpha*dd
        state_out = state_new;
        cost_out = cost_new;
        ds = sqrt(update'*model.mass*update)*alpha;
        return;
    end
    
    if num_iter < model.max_line_search_iter
        alpha = alpha*(1-1/golden_ratio);
    else
        disp('>>>> Warning: Max line search iterations reached');
        if cost_new < cost 
            state_out = state_new;
            cost_out = cost_new;
        else
            state_out = state;
            cost_out = cost;
        end
        ds = sqrt(update'*model.mass*update)*alpha;
        return;
    end
    
end

end