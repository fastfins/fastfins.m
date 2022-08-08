function Qty = matvec_Qty(model, sol)

switch model.type
    case {'exp'}
        tmp = sol.dxdu.*model.F';
        Qty = model.qoi*tmp*sol.d(:);
    case {'linear', '2exp', 'poisson'}
end

end