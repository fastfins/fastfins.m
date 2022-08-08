function sol = rom_solve(rom, vr)

% assemble the boundary matrix
c = rom.x_K*(rom.beta_weight.*exp(rom.redu_u*vr + rom.redu_mean));
tmp = c*c';
x = tmp(rom.inds{3});
Mbnd = reshape(rom.As{3}*x, rom.dof, rom.dof);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

state = rom.b;
[sqrt_eta, b] = eval_pre_data_rom(rom, state);

num_iter = 0;

if rom.iter_disp >= 1
    fprintf('Iter    L2(res)--------------\n');
end

while true
    num_iter = num_iter + 1;
    % stiffness matrix with eta
    tmp = sqrt_eta*sqrt_eta';
    x = tmp(rom.inds{1});
    He = reshape(rom.As{1}*x, rom.dof, rom.dof);
    %
    tmp = b*b';
    x = tmp(rom.inds{2});
    Hb = reshape(rom.As{2}*x, rom.dof, rom.dof);
    
    % residual
    res = He*state + Mbnd*state - rom.b;
    res_l2 = norm(res);
    
    % residual condition check
    if res_l2 < rom.res_tol
        % print line search info, residual condition met
        % ls_num_iter = 0
        if rom.iter_disp >= 1
            fprintf('%4d    %020.13E\n', num_iter, res_l2);
        end
        break;
    end
    
    state = state - (He + Hb + Mbnd)\res;
    
    % iteration check
    if num_iter >= rom.max_newton_iter
        % print line search info, max iter reached
        if rom.iter_disp >= 1
            fprintf('%4d    %020.13E\n', num_iter, res_l2);
        end
        disp('>>>> Warning: Max Newton iterations reached');
        break;
    end
    
    if rom.iter_disp >= 2
        fprintf('%4d    %020.13E\n', num_iter, res_l2);
    end
    
    [sqrt_eta, b] = eval_pre_data_rom(rom, state);
end

sol.state = state;
sol.d = rom.obs_operator*sol.state;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
function [state_out, cost_out, num_iter, ds] = armijo_search(rom, Mbnd, state, cost, grad, update)
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
%}