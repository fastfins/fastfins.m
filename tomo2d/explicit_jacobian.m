function Ju = explicit_jacobian(model, sol)
%

Ju = -sol.d(:).*(model.F.*sol.dxdu(:)');
    
end

% I*exp(-F*exp(u))
% J = I*exp(-F*x)*(-F)*exp(u)