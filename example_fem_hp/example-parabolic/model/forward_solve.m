function sol = forward_solve(model, u)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(u) < model.mesh.dof
    u = reshape(u, 1, []); 
else
    u = reshape(u, model.mesh.dof, []); 
end
    
% determine kappa type
n  = size(u,2);
nd = model.local_elem.num_dim;
if n == 1
    sol.kappa_type = 'scalar';
    if model.exp_param
        x = exp(u);
        sol.dxdu = x;
    else
        x = u;
        sol.dxdu = ones(length(u(:)),1);
    end
    sol.kappa   = x;
elseif n == nd && nd > 1
    sol.kappa_type = 'vector';
    if model.exp_param
        x = exp(u);
        sol.dxdu = x(:);
    else
        x = u;
        sol.dxdu = ones(length(u(:)),1);
    end
    sol.kappa   = x;
elseif n == (nd+1)
    sol.kappa_type = 'tensor';
    x = u;
    dxdu = ones(size(u));
    if model.exp_param
        x(:,1) = exp(u(:,1));
        dxdu(:,1) = x(:,1);
    end
    sol.kappa   = x;
    sol.dxdu    = dxdu(:);
else
    error('kappa not implemented')
end

% assemble stiffness matrix
if model.sq_param
    Ak  = assemble_stiff_sq(model.mesh, model.local_elem, model.grad, model.fill_i, x, model.fill_d);
else
    Ak  = assemble_stiff(model.mesh, model.local_elem, model.grad, model.fill_i, x, model.fill_d);
end
% assemble boundary condition
A   = Ak + model.bmass + model.exp_thres*model.stiff;
% assemble forcing term, do nothing

% Crank-Nicolson
% P*u(t+1) = Q*u(t) + f
sol.Pc  = model.mass + A*(model.dt*0.5);
sol.Qc  = model.mass - A*(model.dt*0.5);
fc  = model.b*model.dt;
%
% initial time steps using implicit Euler with smaller time steps
sol.Pi  = model.mass + A*(model.dt/model.T_lead);
fi  = model.b*(model.dt/model.T_lead);

%
% cholesky 
%   A(p,p) = LL' <=> P'*A*P = LL' <=> A = PLL'P' as PP' = I
%   Au = f <=> (PL)(L'P')u = f <=> Ly = P'f, L'(P'u) = y
%   P'f = f(p)

[sol.Lc,~,sol.pc] = chol(sol.Pc, 'lower', 'vector');
[sol.Li,~,sol.pi] = chol(sol.Pi, 'lower', 'vector');
% solve

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the initial condition is given, the state matrix does not contain T=0
sol.init  = zeros(model.dof, model.T_lead);
sol.state = zeros(model.dof, model.T_nstep);
for t = 1:model.T_lead
    if t == 1
        g = model.mass*model.init + fi; 
    else
        g = model.mass*sol.init(:,t-1) + fi; 
    end
    sol.init(sol.pi,t) = sol.Li'\(sol.Li\g(sol.pi));
end

sol.state(:,1) = sol.init(:,model.T_lead);
for t = 2:model.T_nstep
    g = sol.Qc*sol.state(:,t-1) + fc; 
    sol.state(sol.pc,t) = sol.Lc'\(sol.Lc\g(sol.pc));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% apply observation operator, which contains the mass matrix if needed
sol.d  = model.obs_operator*sol.state(:, model.obs_ind);
sol.d  = sol.d(:);

sol.q = model.pred_operator*sol.state(:, model.pred_ind);
sol.q  = sol.q(:);

% apply qoi function
%{
if model.qoi_flag
    sol.qoi = model.phi*(Ak*sol.state(:, model.pred_ind)); % the const forcing term is not added, does not affect the MC
else
    sol.qoi = [];
end
%}

end


