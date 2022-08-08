function [xmin,solmin] = line_search(obj, hessian, opts, init)
% LINE_SEARCH    perform line seach, using either Wolfe condition or Armijo
%                backtracking, search direction is given by one of inexact
%                Newton, nonlinear CG or BFGS
%
% Tiangang Cui, 16/August/2020

% opt_def

% initialize
xmin = init;
[solmin,gmin] = obj(xmin);

f0 = solmin.f;

step = opts.line_step;

if strcmp(opts.solver, 'BFGS')
    Y = zeros(length(gmin), opts.bfgs_max_num);
    S = zeros(length(gmin), opts.bfgs_max_num);
    j = 0;
else
    Y = [];
    S = [];
end

fprintf('\n\nIteration \t Obj Value \t\t 1st KKT \t\t CG \t\n');
i = 0;
while 1
    % get the next search direction
    switch opts.solver
        case {'Newton'}
            H = hessian(solmin);
            p = -H\gmin;
            iter = 0;
        case {'Newton_CG'}
            [p,iter] = newton_cg(opts, hessian, solmin, gmin);
        case {'NCG'}
            % preconditioned conjugate gradient
        otherwise
            % preconditioend BFGS
            p = -two_loops(Y, S, j, gmin, 'iH');
            iter = 0;
    end
    
    % Armijo
    [xnew, step] = armijo_search(obj, opts, step, xmin, solmin.f, gmin, p);
    step = min(step, 1);
    % adjoint gradient
    [solnew,gnew] = obj(xnew);
    n_gnew = norm(gnew);
    
    % iteration
    i = i+1;
    
    info = iter_info(opts, f0, norm(xnew-xmin), abs(solnew.f-solmin.f), solnew.f, n_gnew, i, iter);
    if info > 0
        xmin = xnew;
        gmin = gnew;
        solmin = solnew;
        return;
    end
    
    if strcmp(opts.solver, 'BFGS')
        % update Y and S
        y = gnew - gmin;
        s = xnew - xmin;
        
        if y'*s > 1e-5
            j = j+1;
            n = min(j, opts.bfgs_max_num);
            if n > 1
                Y(:,2:n) = Y(:,1:(n-1));
                S(:,2:n) = S(:,1:(n-1));
            end
            Y(:,1) = y;
            S(:,1) = s;
        end
    end
    
    % swap state
    xmin = xnew;
    gmin = gnew;
    solmin = solnew;
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nx, nstep] = armijo_search(obj, opts, step, x, f, g, p)
% ARMIJO_SEARCH    performs a line search satisfies the armijo condition
% [next, a] = armijo_search(opt_def, curr)
%
% given a direction of the line search, compute the next point that satify
% the Armijo condition
%
% curr: information at current point, has
%   --    x: current point
%   --    f: function value
%   -- grad: gradient
%   --    p: search direction
%
% find f(x + a*p) <= f(x) + ftol*a*grad_f(x)'*p
%      abs( grad_f(x+a*p)' * p ) <= gtol grad_f(x)'*p
%
% the initial choice of the jump size as 1 for newton and quasi newton,
% there is no maximum jump size, since our problem is propoper.
%
% Tiangang Cui, 15/Mar/2013


f0 = f;    % initial function value
d0 = g'*p; % initial directional derivative
nstep = step;

if d0 > 0
    disp('Bad search direction');
    nx  = x;
    return;
end

golden_ratio = 0.5*(sqrt(5) + 1);

% start up, assign a initial step size and set the function evalaution
% counter to 1
na = step;
i = 1;

%figure
%plot(curr.p)

while 1
    % minus log post value, gradient, and gradient of the likelihood
    % na 
    nsol = obj(x + na*p);
    
    % check for suffcient decrease, if not met
    if nsol.f >= (f0 + opts.line_ftol*na*d0)
        if i <= opts.line_max_feval
            i = i+1;
            na = na*(1-1/golden_ratio);
            continue;
        else
            disp('Max line search iterations reached');
            nx  = x + na*p;
            break;
        end
    end
    % if met
    nx = x + na*p;
    break;
end

if i == 1
    nstep = step*2;
end

end
