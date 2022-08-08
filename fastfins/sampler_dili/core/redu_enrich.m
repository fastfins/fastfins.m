function [gsvd, P, S] = redu_enrich(hessian, options, gsvd, v)
% Reinitialize the reduced subspace
% Tiangang Cui, 01/Oct/2013

% first compute the local subspace, then determine if use the iterative
% update or the full update

switch options.hess_type
    case {'Eig'}
        [V_loc,d_loc] = hessian(v, options.local_trunc_tol, options.local_max_dim);
        s_loc = sqrt(d_loc);
    case {'SVD'}
        [~,s_loc,V_loc] = hessian(v, options.local_trunc_tol, options.local_max_dim);
end

% set redu_def.iterative_update_number to redu_def.gsvd_Nmax if full update
% is used all the time

iter_SVD_flag = true;
%iter_SVD_flag = false;

if iter_SVD_flag
    gD      = gsvd.S(:).^2;
    T       = gsvd.V'*V_loc;
    [Q,R]   = qr(V_loc - gsvd.V*T, 0);
    tmp1    = [T; R].*s_loc(:)';
    tmp2    = tmp1*tmp1' + diag([gD*gsvd.n; zeros(length(s_loc), 1)]);
    [Phi,D,~] = svd(tmp2, 0); % this line may cause numerical instability
    [S_glo,ind] = sort( sqrt(diag(D)) / sqrt(gsvd.n + 1), 'descend' );
    V_glo   = [gsvd.V, Q]*Phi(:, ind);
else % full update
    Vmat          = [(gsvd.V.*gsvd.S(:)')*sqrt(gsvd.n), V_loc.*s_loc(:)'];
    [V_glo,S_glo] = svd(Vmat/sqrt(gsvd.n+1),'econ'); % global SVD
    S_glo         = diag(S_glo);
    
end

r1 = min(options.global_max_dim, sum(S_glo>=options.global_trunc_tol));

% this is the svd
gsvd.V  = V_glo(:,1:r1);
gsvd.S  = S_glo(1:r1);
gsvd.n  = gsvd.n + 1;
gsvd.dof = r1;

r2 = sum(gsvd.S.^2>=options.lis_trunc_tol);
P  = gsvd.V(:,1:r2);
S  = gsvd.S(1:r2);

end

