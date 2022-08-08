function prec = update_mrf(mesh, k, cond, sigma)
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
temp_K      = mesh.loc_xx(:)*temp(:,1)' + mesh.loc_yy(:)*temp(:, 2)' + mesh.loc_xy(:)*temp(:,3)';
K           = sparse(mesh.ind_i(:), mesh.ind_j(:), temp_K(:), mesh.Nnode, mesh.Nnode);

nor         = gamma(1)/(gamma(2)*(k^2)*(4*pi)); % normalizing const of kernel
tmp         = sqrt(nor)/sigma;                  % std of the covariance

spML        = spdiags(mesh.ML,        0,mesh.Nnode,mesh.Nnode);
%sqspML     = spdiags(mesh.ML.^0.5,   0,mesh.Nnode,mesh.Nnode);
ispML       = spdiags(mesh.ML.^(-1),  0,mesh.Nnode,mesh.Nnode);
isqspML     = spdiags(mesh.ML.^(-0.5),0,mesh.Nnode,mesh.Nnode);

%P  = (K + c) + spML*(k^2);
P           = K + spML*(k^2);
prec.per    = symamd(P); % this is a MatLab function

% prec.R and prec.sca are used for evaluating the 
prec.R      = chol(P(prec.per,prec.per)); % upper triangular cholesky of the P matrix
prec.sca    = sqrt(mesh.ML(prec.per))/tmp; % sqrt of the scale matrix

prec.RQ     = tmp*isqspML*P; % prec.RQ' * prec.RQ = Q, sqrt of the precision
prec.Q      = P*(ispML*(tmp^2))*P;  % precision matrix
prec.DoF    = mesh.Nnode;
prec.cov_type = 'MRF';
prec.type     = 'Field';

end