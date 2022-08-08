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

n = size(mesh.nodes,1);
C = zeros(n);

for i   = 1:n
    x = mesh.nodes(i,:);
    d = (mesh.nodes - repmat(x,n,1));
    C(:,i)  = sum((d*scale).*d, 2);
end

if power == 2
    C = (exp(-0.5*C) + 1e-10*eye(n));
else
    C = exp(-0.5*C.^(0.5*power));
end

prior.C = C*sigma^2;
prior.L = chol(prior.C, 'lower');
prior.cov_type = 'GP';
prior.type = 'Field';
prior.dof  = n;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%