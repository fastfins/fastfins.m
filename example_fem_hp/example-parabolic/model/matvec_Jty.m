function Jty = matvec_Jty(model, sol, dy)

dy  = reshape(dy, model.n_sensors, []);
N   = size(dy,2)/model.n_datasets;
Jty = zeros(length(sol.dxdu(:)), N);

lambda = zeros(model.dof, 1);
for k = 1:N
    ind = (k-1)*model.n_datasets + (1:model.n_datasets);
    dU  = zeros(size(model.obs_operator,2), model.T_obs_nstep);
    dU(:,model.obs_ind) = model.obs_operator'*dy(:,ind);
    lambda(:) = 0;
    for t = model.T_obs_nstep:-1:2
        if t == model.T_obs_nstep
            g = - dU(:,t);
        else
            g = sol.Qc*lambda - dU(:,t);
        end
        lambda(sol.pc) = sol.Lc'\(sol.Lc\g(sol.pc));
        %
        tmp = (sol.state(:,t-1) +  sol.state(:,t))*(model.dt*0.5);
        if model.sq_param
            Jty(:,k) = Jty(:,k) + deri_adjoint_stiff_sol_sq(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, lambda, tmp);
        else
            Jty(:,k) = Jty(:,k) + deri_adjoint_stiff_sol(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, lambda, tmp);
        end
    end
    %
    for t = model.T_lead:-1:1
        if t == model.T_lead
            g = sol.Qc*lambda - dU(:,1);
        else
            g = model.mass*lambda;
        end
        lambda(sol.pi) = sol.Li'\(sol.Li\g(sol.pi));
        %
        tmp = sol.init(:,t)*(model.dt/model.T_lead);
        if model.sq_param
            Jty(:,k) = Jty(:,k) + deri_adjoint_stiff_sol_sq(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, lambda, tmp);
        else
            Jty(:,k) = Jty(:,k) + deri_adjoint_stiff_sol(model.mesh, model.local_elem, sol.kappa_type, sol.kappa, lambda, tmp);
        end
    end
end

% transform to physical space x
Jty = sol.dxdu.*Jty;

end

