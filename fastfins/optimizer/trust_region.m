function [xmin,solmin] = trust_region(obj, hessian, opts, init)
% TRUST_REGION   Subspace inexact Newton trust region solver
%
% Tiangang Cui, 16/August/2020
% opt_def
% initialize
xmin = init;
[solmin,gmin] = obj(xmin);
n_gmin = norm(gmin);

f0 = solmin.f;

if strcmp(opts.solver, 'BFGS')
    Y = zeros(length(gmin), opts.bfgs_max_num);
    S = zeros(length(gmin), opts.bfgs_max_num);
    j = 0;
else
    Y = [];
    S = [];
end

% initial trust region
delta = 0.05*opts.TR_radius;

fprintf('\n\nIteration \t Obj Value \t\t 1st KKT \t\t CG \t\n');
i = 0;
while 1
    % subspace search
    switch opts.solver
        case {'Newton'}            
            H = hessian(solminmin);
            p = -H\gmin;
            % subspace is g and p
            u = -gmin/n_gmin;
            v = p - u*(u'*p);
            nv = norm(v);
            if nv > 1e-10
                v = v/nv;
                U = [u, v];
                B = U'*(H*U);
                a = U'*gmin;
                [s, quad_decrease] = twod_subspace(a, B, delta);
                jump = U*s;
            else
                d = -(delta/n_gmin)*gmin;
                a = 0.5*d'*(H*d);
                b = gmin'*d;
                [jump, quad_decrease] = cauchy(a, b, d);
            end
        case {'Newton_CG'}
            [p,iter,pdef] = newton_cg(opts, hessian, solmin, gmin);
            %if p'*curr.grad>0
            %    disp('bad');
            %    return;
            %end
            % subspace is g and p
            if pdef || iter > 1
                % solve the 2d problem
                u = -gmin/n_gmin;
                v = p - u*(u'*p);
                nv = norm(v);
            else
                nv = -1;
            end
            if nv > 1e-10
                v = v/nv;
                U = [u, v];
                B = U'*hessian(solmin, U);
                a = U'*gmin;
                [s, quad_decrease] = twod_subspace(a, B, delta);
                jump = U*s;
            else
                d = -(delta/n_gmin)*gmin;
                a = 0.5*d'*hessian(solmin, d);
                b = gmin'*d;
                [jump, quad_decrease] = cauchy(a, b, d);
            end
        otherwise
            p = -two_loops(Y, S, j, gmin, 'iH');
            iter = 0;
            
            u = -gmin/n_gmin;
            v = p - u*(u'*p);
            nv = norm(v);
            if nv > 1e-5
                v = v/nv;
                U = [u, v];
                tmp1 = two_loops(Y, S, j, U(:,1), 'H');
                tmp2 = two_loops(Y, S, j, U(:,2), 'H');
                B = U'*[tmp1, tmp2];
                a = U'*gmin;
                [s, quad_decrease] = twod_subspace(a, B, delta);
                jump = U*s;
            else
                d = -(delta/n_gmin)*gmin;
                a = 0.5*d'*two_loops(Y, S, j, d, 'H');
                b = gmin'*d;
                [jump, quad_decrease] = cauchy(a, b, d);
            end
    end
    
    
    % next point
    xnew = xmin + jump;
    [solnew,gnew] = obj(xnew);
    n_gnew = norm(gnew);
    
    i = i + 1;
    info = iter_info(opts, f0, norm(jump), abs(solnew.f-solmin.f), solnew.f, n_gnew, i, iter);
    if info > 0
        xmin = xnew;
        solmin = solnew;
        return;
    end
    
    % evaluate the gain
    rho = (solnew.f-solmin.f)/quad_decrease;
    
    if rho < 1E-10
        delta = delta*0.1;
    elseif rho < 0.25
        delta = 0.25*delta;
    else
        if norm(jump) > (delta - 1e-2)
            delta = min(opts.TR_radius, 2*delta);
        end
    end
    
    % update
    if rho > 0.1
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
        n_gmin = n_gnew;
        solmin = solnew;
    end
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [s, quad_decrease] = twod_subspace(a, B, delta)

% min  g'Vs + 0.5*s'*V'*H*V*s
% s.t. s'V'Vs <= delta^2
%
% min a'*s + 0.5 s'*B*s, s.t. s'*C*s <= delta^2

% C = U'*U; % C is identity

% solve the unconstrained problem
s = -B\a;

% check if this violate the constrain
if s'*s <= delta^2
    quad_decrease = a'*s + 0.5*s'*B*s;
    return;
end

% then reparametrize
% s1 = delta * cos t;
% s2 = delta * sin t;
%
% the model we want to minimize is
% min a1*s1 + a2*s2 + 0.5*(B11*s1^2 + B22*s2^2 + 2*B12s1*s2)
% min delta*a1*cost + delta*a2*sint + (0.5*delta^2) * (B11*cost^2 +
% B22^sint^2 + 2*B12*sintcost)
% t in [0, 2*pi)

% coarse scale searching
ts = linspace(0,2*pi,100);
fs = sub_con(a, B, delta, ts);
% find a starting point
[~, ind] = min(fs);
ct = ts(ind);
% Newton solve

for i = 1:100
    [df1, df2] = sub_con_df(a, B, delta, ct);
    nt = ct - df1/df2;
    if abs(nt-ct) < 1e-15
        break;
    else
        ct = nt;
    end
end
t = nt;
        
s = delta*[cos(t); sin(t)];
quad_decrease = a'*s + 0.5*s'*B*s;

end


function f = sub_con(a, B, delta, t)

f = delta*a(1)*cos(t) + delta*a(2)*sin(t) + ...
    (0.5*delta^2) * ( B(1,1)*cos(t).^2 + B(2,2)*sin(t).^2 + 2*B(1,2)*sin(t).*cos(t) );

end


function [df1, df2] = sub_con_df(a, B, delta, t)

df1 = -delta*a(1)*sin(t) + delta*a(2)*cos(t) + ...
    (0.5*delta^2) * ( - 2*B(1,1)*cos(t).*sin(t) + 2*B(2,2)*sin(t).*cos(t) - ...
    2*B(1,2)*sin(t).^2 + 2*B(1,2)*cos(t).^2 );

df2 = -delta*a(1)*cos(t) - delta*a(2)*sin(t) + ...
    (0.5*delta^2) * ( - 2*B(1,1)*(cos(t)^2 - sin(t)^2) + 2*B(2,2)*(cos(t)^2 - sin(t)^2) - ...
    4*B(1,2)*cos(t)*sin(t) - 4*B(1,2)*cos(t)*sin(t) );

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [jump, quad_decrease] = cauchy(a, b, d)

%d = -(delta/curr.n_grad)*curr.grad;
% p = alpha*d, alpha in [0, 1]
% the local model is
% g'*d*alpha = 0.5*alpha^2*d'*B*d
%
% min a*alpha^2 + b*alpha, alpha in [0, 1]
%
%b = curr.grad'*d;
%a = 0.5*d'*opts.hessian(curr.hessinfo, d);

if a > 0
    alpha = -b/(2*a);
    if alpha < 1 && alpha > 0
        jump = d*alpha;
        quad_decrease = a*alpha^2 + b*alpha;
    elseif (a+b) < 0
        jump = d;
        quad_decrease = a+b;
    else
        jump = 0;
        quad_decrease = 0;
        disp('Cauchy solve failed: not move')
    end
else
    if (a+b) < 0
        jump = d;
        quad_decrease = a+b;
    else
        jump = 0;
        quad_decrease = 0;
        disp('Cauchy solve failed: not move')
    end
end

end
