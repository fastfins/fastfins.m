function logd =  rto_func_d(model, obs, curr_lin, HI, trust_flag)
%
% Tiangang Cui, August, 1, 2017

if model.explicit_ja
    A   = HI.A + eye(curr_lin.nrank);
else
    % compute determinant
    if  curr_lin.dof >= obs.n_data
        m   = size(curr_lin.LPhi, 2);
        %du      = curr_lin.LPhi;
        JPhi    = zeros(obs.n_data, m);
        for i = 1:m
            JPhi(:,i)   = (matvec_Ju(model, HI, curr_lin.LPhi(:,i))./obs.std);
        end
        
        if trust_flag
            JPhi    = HI.diag(:).*(JPhi - curr_lin.linref) + curr_lin.linref;
        end
        
        A   = curr_lin.linref'*JPhi;
    else
        m   = size(curr_lin.Psi, 2);
        JtPsi   = zeros(curr_lin.nrank, m);
        for i = 1:m
            if trust_flag
                d   = HI.diag(:).*curr_lin.Psi(:,i);
                JtPsi(:,i)  = curr_lin.LPhi'*(matvec_Jty(model, HI, d./obs.std)) ...
                    + curr_lin.linref'*(curr_lin.Psi(:,i) - d);
            else
                gx          = matvec_Jty(model, HI, curr_lin.Psi(:,i)./obs.std);
                JtPsi(:,i)  = curr_lin.LPhi'*gx;
            end
        end        
        
        A   = curr_lin.s(:).*JtPsi';
    end
    A   = A + eye(curr_lin.nrank);
end

if cond(A) > 1E20
    disp('Bad condition number, fail to evaluate the determinant of Jacobian')
end

s       = svd(A);
logd    = sum(log(s));

end