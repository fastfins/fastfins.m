function delta = delta_sample(n, c0, alpha, beta, left, right)
% GAMSAMP_FUN. This function computes a sample from the unnormalized
% probablity density function g(gamma)

g = @(rho) ( rho + 0.5*n*rho - c0*exp(rho) ) ...
    + alpha*log(exp(rho) - left) + beta*log(right - exp(rho));

rhos    = reshape(linspace(log(left+eps),log(right-eps),5000), 1, []);
%
hs      = [rhos(2)-rhos(1), rhos(2:end)-rhos(1:end-1), rhos(end)-rhos(end-1)];
weight  = 0.5*(hs(2:end)+hs(1:end-1));
%
logg    = arrayfun(g,rhos);
gpmf    = exp(logg - max(logg)).*weight;
gCDF    = cumsum(gpmf)/sum(gpmf);
delta   = exp(min(rhos(gCDF>rand)));

end

%{

g = @(rho) ( rho + 0.5*n*rho - c0*exp(rho) ) ...
    + alpha*log(exp(rho) - left) + beta*log(right - exp(rho));

rhos    = reshape(linspace(log(left+eps),log(right-eps),5000), 1, []);
%
hs      = [rhos(2)-rhos(1), rhos(2:end)-rhos(1:end-1), rhos(end)-rhos(end-1)];
weight  = 0.5*(hs(2:end)+hs(1:end-1));
%
logg    = arrayfun(g,rhos);
gpmf    = exp(logg - max(logg));
gCDF    = cumsum(gpmf)/sum(gpmf).*weight;
delta   = exp(min(rhos(gCDF>rand)));

%}

%{

g = @(rho) 0.5*n*log(rho) - c0*rho ...
    + alpha*log(rho - left) + beta*log(right - rho);

rhos    = exp(reshape(linspace(log(left+eps),log(right-eps),5000), 1, []));
%
hs      = [rhos(2)-rhos(1), rhos(2:end)-rhos(1:end-1), rhos(end)-rhos(end-1)];
weight  = 0.5*(hs(2:end)+hs(1:end-1));
%
logg    = arrayfun(g,rhos);
gpmf    = exp(logg - max(logg)).*weight;
gCDF    = cumsum(gpmf)/sum(gpmf);
delta   = min(rhos(gCDF>rand));

%}