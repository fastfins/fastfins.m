function A = build_stiff_matrix(model, u)

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
    if model.exp_kappa
        x = exp(u);
        sol.dxdu = x;
    else
        x = u;
        sol.dxdu = ones(size(u));
    end
    sol.kappa_b = [];
elseif n == nd && nd > 1
    sol.kappa_type = 'vector';
    if model.exp_kappa
        x = exp(u);
        sol.dxdu = x(:);
    else
        x = u;
        sol.dxdu = ones(size(u));
        sol.dxdu = sol.dxdu(:);
    end
    sol.kappa_b = [];
elseif n == (nd+1)
    sol.kappa_type = 'tensor';
    x = u;
    dxdu = ones(size(u));
    if model.exp_kappa
        x(:,1) = exp(u(:,1));
        dxdu(:,1) = x(:,1);
    end
    sol.dxdu    = dxdu(:);
    sol.kappa_b = u(:,2:end);
else
    error('kappa not implemented')
end

% assemble stiffness matrix
Ak  = assemble_stiff(model.mesh, model.local_elem, model.grad, model.fill_i, x, model.fill_d);
% assemble boundary condition
A   = Ak + model.bmass + model.exp_thres*model.stiff;
end