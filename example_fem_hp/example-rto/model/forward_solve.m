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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cholesky 
%   A(p,p) = LL' <=> P'*A*P = LL' <=> A = PLL'P' as PP' = I
%   Au = f <=> (PL)(L'P')u = f <=> Ly = P'f, L'(P'u) = y
%   P'f = f(p)

[sol.L,~,sol.p] = chol(A,'lower', 'vector');
% solve
sol.state = zeros(size(model.b));
sol.state(sol.p,:)  = sol.L'\(sol.L\model.b(sol.p,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% apply observation operator, which contains the mass matrix if needed
sol.d  = model.obs_operator*sol.state;

% apply qoi function
if model.qoi_flag
    sol.Q   = -model.phi*Ak;
    sol.qoi = sol.Q*sol.state; % the const forcing term is not added, does not affect the MC
else
    sol.qoi = [];
end

end


