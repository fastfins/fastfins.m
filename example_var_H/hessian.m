function H = hessian(model, prior, sol, d)

switch model.type
    case {'linear'}
        H = model.F'*model.F;
    case {'exp'}
        J = model.F.*sol.dxdu(:)';
        H = J'*J;
        for i = 1:numel(d)
            H = H + diag(J(i,:)) * d(i);
        end
    case {'2exp'}
        tmp = model.F.*sol.dxdu(:)';
        J = sol.d(:).*tmp';
        H = J'*J;
        for i = 1:numel(d)
            H = H + sol.d(i) * ( tmp(i,:)'*tmp(i,:) + diag(tmp(i,:)) ) * d(i);
        end
    case {'poisson'}
        disp('not implemented')
end

H = prior.D(:)'.*H.*prior.D(:);

end


