
categ = [-1, 0, 1];

m = length(categ);

r = randn(prior.dof, m);

phi = prior_cov_l(prior, r);

g = exp(phi)./sum(exp(phi),2);

R = mnrnd(1,g);

u = R*categ';

trisurf(model.mesh.node_tri, model.mesh.nodes(:,1), model.mesh.nodes(:,2), u, 'edgecolor', 'none')