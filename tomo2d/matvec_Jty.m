function Jty = matvec_Jty(model, sol, dy)
%

dy = reshape(dy, length(sol.d(:)), []);
Jty = -sol.dxdu(:).*model.F'*( dy.*sol.d(:) );
    
end

% I*exp(-F*exp(u))
% J = I*exp(-F*x)*(-F)*exp(u)