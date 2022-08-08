function gamma = gamma_sample_2d(chiAs, c1, c2, alpha, beta, left, right)
% GAMSAMP_FUN. This function computes a sample from the unnormalized
% probablity density function g(gamma)

% c1 = 0.5*delta*u'*M*u;
% c2 = delta*u'*K*u;

g = @(rho) ( rho + sum( log(chiAs(:) + exp(rho)) ) - c1*exp(2*rho) - c2*exp(rho) ) ...
    + alpha*log(exp(rho) - left) + beta*log(right - exp(rho));

rhos    = reshape(linspace(log(left+eps),log(right-eps),5000), 1, []);
%
hs      = [rhos(2)-rhos(1), rhos(2:end)-rhos(1:end-1), rhos(end)-rhos(end-1)];
weight  = 0.5*(hs(2:end)+hs(1:end-1));
%
logg    = arrayfun(g,rhos);
gpmf    = exp(logg - max(logg)).*weight;
gCDF    = cumsum(gpmf)/sum(gpmf);
gamma   = exp(min(rhos(gCDF>rand)));

end

%{

g = @(rho) ( rho + sum( log(chiAs(:) + exp(rho)) ) - c1*exp(2*rho) - c2*exp(rho) ) ...
    + alpha*log(exp(rho) - left) + beta*log(right - exp(rho));

rhos    = reshape(linspace(log(left+eps),log(right-eps),5000), 1, []);
%
hs      = [rhos(2)-rhos(1), rhos(2:end)-rhos(1:end-1), rhos(end)-rhos(end-1)];
weight  = 0.5*(hs(2:end)+hs(1:end-1));
%
logg    = arrayfun(g,rhos);
gpmf    = exp(logg - max(logg)).*weight;
gCDF    = cumsum(gpmf)/sum(gpmf);
gamma   = exp(min(rhos(gCDF>rand)));

%}

%{

g = @(rho) sum( log(chiAs(:) + rho), 1 ) - c1*rho.^2 - c2*rho ...
    + alpha*log(rho - left) + beta*log(right - rho);

rhos    = exp( reshape(linspace(log(left+eps),log(right-eps),5000), 1, []) );
%
hs      = [rhos(2)-rhos(1), rhos(2:end)-rhos(1:end-1), rhos(end)-rhos(end-1)];
weight  = 0.5*(hs(2:end)+hs(1:end-1));
%
logg    = arrayfun(g,rhos);
gpmf    = exp(logg - max(logg)).*weight;
gCDF    = cumsum(gpmf)/sum(gpmf);
gamma   = min(rhos(gCDF>rand));

%}