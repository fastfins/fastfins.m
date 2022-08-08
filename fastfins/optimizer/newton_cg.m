function [x,i,pdef] = newton_cg(opts, hessmult, sol, g)
% NEWTON_CG_LINE   Conjugate gradient inexact Newton solver, without BFGS
%                  preconditioner
%
% delta > 0 gives the Steihaug CG
%
% Tiangang Cui, 16/August/2020

% trust region info

pdef = true;

x  = zeros(size(g));
r  = g;
ng = norm(g);

% congugate direction
cd = -r;
% r'Pr
top = r'*r;

forcing = min(opts.CG_forcing_tol, sqrt(ng));
forcing = max(0.1, forcing)*ng;

% start
i = 0;
while 1
    % Hessmult, apply to the congugate direction
    Hcd = hessmult(sol, cd);
    % curvature
    bottom = cd'*Hcd;
    
    % negative curvature, use the inexact Newton trick
    if bottom/(cd'*cd) < opts.CG_zero_tol
    % if bottom < 0
        pdef = false;
        if i == 0
            x = cd;
        else
            alpha = top/bottom;
            x = x + alpha*cd;
        end
        break;
    end
    
    % update x and residual
    alpha = top/bottom;
    x = x + alpha*cd;
    r = r + alpha*Hcd;
    
    i = i + 1;
    
    % recompute residual after certain number of iterations
    if mod(i, opts.CG_restart) == 0
        r = g + hessmult(sol, x);
    end
    
    % check residual
    if norm(r) < forcing
        break;
    end
    
    % check the iteration count
    if i >= opts.CG_max_iter
        disp('Maximum iteration reached, terminating CG')
        break;
    end
    
    ntop = r'*r;
    % update conjugate direction
    cd = - r + (ntop/top)*cd;
    
    % next iteration
    top = ntop;
end

end
