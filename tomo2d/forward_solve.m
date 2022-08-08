function sol = forward_solve(model, u)

[x, sol.dxdu] = u2x( model.func, u );
sol.d = model.Is * exp( -model.F*x );

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x, dxdu] = u2x( func, u )
%U2X
% Transform the computation parameter u ( with prior N(m, C) ) to the
% physical parameter x used by the PDE
%
% Tiangang Cui, 17/Jan/2014
%
% function [x, dxdu, dx2du2] = prior_u2x( func, u )
% For the full Hessian

switch func.type
    case{'exp'}
        temp    = exp(u);
        x       = temp  + func.log_thres;
        dxdu    = temp;
        %dx2du2  = temp;
    case{'erf'}
        x       = erf(u).*func.erf_scale + func.erf_shift;
        dxdu    = func.erf_scale.*2./sqrt(pi).*exp(-u(:).^2);
        %dx2du2  = -2*dxdu.*u(:);
    otherwise
        x       = u;% + prior.shift;
        dxdu    = ones(size(x));
        %dx2du2  = zeros(size(x));
end

end
