function prior = make_prior_conv(mesh, centers, radii, weights, kernel_func)

n  = size(centers, 2);
np = size(mesh.node, 2);
prior.basis = zeros(np, n);

for i = 1:n
    x = centers(:,i);
    d = mesh.node - repmat(x,1,np);
    s = sqrt(sum(d.^2, 1))/radii(i);
    prior.basis(:,i) = kernel_func(s(:))*weights(i);
end

prior.type = 'Basis';
prior.mesh = mesh;
prior.DoF  = n;

end

