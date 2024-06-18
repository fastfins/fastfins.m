function sol = forward_solve(model, u)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(u) < model.mesh.dof
    u = reshape(u, 1, []); 
else
    u = reshape(u, model.mesh.dof, []); 
end
    
if model.exp_param
    x = exp(u);
    sol.dxdu = x;
    sol.dx2du2 = x;
else
    x = u;
    sol.dxdu = ones(length(u(:)),1);
    sol.dx2du2 = zeros(length(u(:)),1);
end
sol.mu_a  = x;
sol.kappa = 1./(2*(sol.mu_a + 0.2*(model.mu_s))); 

% assemble mass matrix, interior
M = assemble_mass(model.mesh, model.local_elem,  model.fill_i, sol.mu_a,  model.fill_d);
% assemble stiffness matrix, interior
K  = assemble_stiff(model.mesh, model.local_elem, model.grad, model.fill_i, sol.kappa, model.fill_d);
% assemble boundary condition
A  = M + K + model.M_bnd*(2/pi);
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


