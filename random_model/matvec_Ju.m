function Ju = matvec_Ju(model, sol, du)

switch model.type
    case {'linear'}
        Ju = model.F*du;
    case {'exp'}
        Ju = model.F*(sol.dxdu.*du);
    case {'2exp'}
        Ju = sol.d(:).*(model.F*(sol.dxdu.*du));
    case {'poisson'}
        %Ju = - sol.d.*(model.F*du);
        Ju = model.F*du;
end

end