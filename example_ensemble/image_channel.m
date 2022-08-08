function u = image_channel(mesh, low, high, width)

channel = @(x) sin(10*x-1)./(x-0.1)/50 + 0.4;

y = channel(mesh.nodes(:,1));
mask = abs(mesh.nodes(:,2) - y) < width;

u = ones(mesh.dof,1)*low;
u(mask) = high;

end