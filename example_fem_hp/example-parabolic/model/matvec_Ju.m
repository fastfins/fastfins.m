function Ju = matvec_Ju(model, sol, du)

% transform to physical space x
dx = sol.dxdu.*du;
Ju = zeros(model.n_sensors, model.n_datasets*size(dx,2));
dstate = zeros(model.dof, model.T_obs_nstep);

% Crank-Nicolson:
% P(k) du(t) = Q(k) du(t-1) + Q(dk) u(t-1) - P(dk) u(t)
%            = Q(k) du(t-1) - dt/2 A(dk) u(t-1) - dt/2 A(dk) u(t)
%            = Q(k) du(t-1) - dt/2 A(dk) (u(t-1) + u(t))
%
% Implicit Euler:
% P(k) du(t) = M du(t-1) - P(dk) u(t)
%            = M du(t-1) - dt A(dk) u(t)

for i = 1:size(dx,2)
    dk = reshape(dx(:,i), model.mesh.dof, []);
    for t = 1:model.T_lead
        if t == 1
            g = zeros(model.dof,1);
        else
            g = model.mass*dstate(:,1);
        end
        tmp = sol.init(:,t)*(model.dt/model.T_lead);
        if model.sq_param
            g = g - matvec_dstiff_sol_sq(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, dk, tmp);
        else
            g = g - matvec_dstiff_sol(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, dk, tmp);
        end
        dstate(sol.pi,1) = sol.Li'\(sol.Li\g(sol.pi));
    end
    %
    for t = 2:model.T_obs_nstep
        g = sol.Qc*dstate(:,t-1);
        %
        tmp = (sol.state(:,t-1) +  sol.state(:,t))*(model.dt*0.5);
        if model.sq_param
            g = g - matvec_dstiff_sol_sq(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, dk, tmp);
        else
            g = g - matvec_dstiff_sol(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, dk, tmp);
        end
        dstate(sol.pc,t) = sol.Lc'\(sol.Lc\g(sol.pc,:));
    end
    ind = (1:model.n_datasets) + (i-1)*model.n_datasets;
    Ju(:,ind) = model.obs_operator*dstate(:,model.obs_ind);
end

Ju = reshape(Ju, [], size(dx,2));

end
