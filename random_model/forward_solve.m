function sol = forward_solve(model, u)

switch model.type
    case {'linear'}
        sol.dxdu = ones(size(u));
        sol.d = model.F*u;
    case {'exp'}
        x = exp(u);
        sol.dxdu = x;
        sol.d = model.F*x;
        % add a qoi
        sol.qoi = 0.5*model.qoi*norm(sol.d)^2;
    case {'2exp'}
        x = exp(u);
        sol.dxdu = x;
        sol.d = exp(model.F*x);
    case {'poisson'}
        x = exp(u);
        sol.dxdu = x;
        %sol.d = model.Is*exp(-model.F*x);
        sol.d = model.F*x;
end


end