function Ju = matvec_Ju(model, sol, du)
%

Ju = - sol.d(:).*(model.F*(sol.dxdu(:).*du));
    
end

% I*exp(-F*exp(u))
% J = I*exp(-F*x)*(-F)*exp(u)