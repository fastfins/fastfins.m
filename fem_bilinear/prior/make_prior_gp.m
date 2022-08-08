function prior = make_prior_gp(mesh, scale, power, sigma)
%MAKE_PRIOR_GP
% makes the Gaussian prior distribution for a given mesh and a correlation length
%
%%%%%%%%%%%%%%%%%%%% input: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mesh:   
% scale:        the tensor for scaling correlation length
% power:        the power term in the kernel
% sigma:        the standard deviation
%
% Tiangang Cui, 12/May/2012

n = size(mesh.node,2);
C = zeros(n);

for i   = 1:n
    x = mesh.node(:,i);
    d = (mesh.node - repmat(x,1,n));
    C(:,i)  = sum((scale*d).*d, 1);
end

if power == 2
    C = (exp(-0.5*C) + 1e-10*eye(n));
else
    C = exp(-0.5*C.^(0.5*power));
end

prior.C  = C*sigma^2;
prior.RC = chol(prior.C);
%prior.Q = inv(prior.C);
prior.cov_type  = 'GP';
prior.type = 'Field';
prior.mesh = mesh;
prior.DoF  = n;

%{
[V D] = eigs(C,300);

d = diag(D);

ind = d>1e-12;

Cv = V(:,ind);
Ce = d(ind);
%}

end

% -0.5*(|x-y|/s)^p*scale

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%