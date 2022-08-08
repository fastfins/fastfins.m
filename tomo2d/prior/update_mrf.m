function prec = update_mrf(mesh, k, cond, delta)
%UPDATE_MRF
%
% Q u = w, Q = inv(sigma) * (K(cond) + k I)
%
% Tiangang Cui, 09/May/2014

if length(cond) == 3 % for homogeneous case
    temp    = repmat(cond(:)', mesh.Nel, 1);
else
    temp    = cond;
end
temp_K  = mesh.loc_xx(:)*temp(:,1)' + mesh.loc_yy(:)*temp(:, 2)' + mesh.loc_xy(:)*temp(:,3)';
prec.K  = sparse(mesh.ind_i(:), mesh.ind_j(:), temp_K(:), mesh.Nnode, mesh.Nnode);

prec.M    = spdiags(mesh.ML,        0,mesh.Nnode,mesh.Nnode);
prec.Minv = spdiags(mesh.ML.^(-1),  0,mesh.Nnode,mesh.Nnode);
prec.Misq = spdiags(mesh.ML.^(-0.5),0,mesh.Nnode,mesh.Nnode);

%P  = (K + c) + spML*(k^2);
prec.P      = prec.K + prec.M*k;
prec.per    = symamd(prec.P); % this is a MatLab function

% prec.R and prec.sca are used for evaluating the 
prec.R      = chol(prec.P(prec.per,prec.per)); % upper triangular cholesky of the P matrix
prec.sca    = sqrt(mesh.ML(prec.per))/sqrt(delta); % sqrt of the scale matrix

prec.RQ     = sqrt(delta)*prec.Misq*prec.P; % prec.RQ' * prec.RQ = Q, sqrt of the precision
prec.Q      = prec.P*(prec.Minv*delta)*prec.P;  % precision matrix
prec.dof    = mesh.Nnode;
prec.cov_type = 'MRF';
prec.type     = 'Field';


end