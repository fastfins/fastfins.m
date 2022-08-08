function lap = get_laplace(hessian, HI)

[lap.V, d]    = hessian(HI);

% w.r.t. llkd
ipseudo       = (d(:)+1).^(-1)   - 1;
ipseudo_half  = (d(:)+1).^(-0.5) - 1;

lap.V_ihalf   = lap.V.*ipseudo_half(:)';  % V*( (Lambda+I)^(-0.5) - I ), for random number
lap.V_half    = lap.V.*sqrt(d(:))';       % V*Lambda^(0.5), for density
lap.V_iH      = lap.V.*ipseudo(:)';       % V*( (Lambda+I)^(-1) - I ), for scale the gradient

lap.ldet_half = 0.5*sum(log(d+1));  % The determinant

end
