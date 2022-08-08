function Jty = matvec_Jty(model, sol, dy)

switch model.type
    case {'linear'}
        Jty = model.F'*dy;
    case {'exp'}
        Jty = model.F'*dy;
        Jty = sol.dxdu.*Jty;
    case {'2exp'}
        Jty = model.F'*(sol.d.*dy);
        Jty = sol.dxdu.*Jty;
    case {'poisson'}
        %Jty = - model.F'*(sol.d.*dy);
        %Jty = sol.dxdu.*Jty;
        Jty = model.F'*dy;
        Jty = sol.dxdu.*Jty;
end

end