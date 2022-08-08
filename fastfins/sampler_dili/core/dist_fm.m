function d = dist_fm(P1, S1, P2, S2)
% compute the distance of the Gaussians
% Tiangang Cui, 08/August/2013

% construct the common basis
[V,~] = qr([P1, P2], 0);

% project the basis onto the common basis
Q1  = (V'*P1).*S1(:)';
Q2  = (V'*P2).*S2(:)';

% construct the kernel on the projected space
A1  = Q1*Q1' + eye(size(Q1, 1));
A2  = Q2*Q2' + eye(size(Q2, 1));

% generalized eigenvalues
lambda  = eig(A1, A2);
d   = sum(log(lambda).^2);

end

%{
D1 = ( gsvd1.S.^2 + 1 ).^(-1)   - 1;
D2 = ( gsvd2.S.^2 + 1 ).^(0.5) - 1;

A = gsvd2.P'*gsvd1.P;
% B = A*diag(D1)*A';
B = scale_cols(A, D1)*A';

% E = D2*B*D2 + B*D2 + D2*B + D2*D2 + 2*D2 + B;
E = scale_rows(scale_cols(B,D2), D2) + scale_cols(B,D2) + scale_rows(B, D2) + diag(D2.^2 + 2*D2) + B;

lambda = eig(E) + 1;

d = sum( log(lambda).^2 );
%}

%A1 = (P2'*P1);
%A2 = (P1'*P2);
%r1 = (1 - diag(A1'*A1)).*param1.S(ind1)/max(param1.S(ind1));
%r2 = (1 - diag(A2'*A2)).*param2.S(ind2)/max(param2.S(ind2));

%{

A = P2'*P1;
r1 = (1 - diag(A'*A)).*param1.S/max(param1.S);
r2 = (1 - diag(A*A')).*param2.S/max(param2.S);

d = sqrt(sum(r1.^2)) + sqrt(sum(r2.^2));

%}

%{
% trim the dimension
N  = min(param1.dof, param2.dof);
A  = P2(:,1:N)'*P1(:,1:N);
r1 = (1 - diag(A'*A)).*param1.S(1:N)/param1.S(1);
r2 = (1 - diag(A*A')).*param2.S(1:N)/param2.S(1);

d  = sqrt(sum(r1.^2)) + sqrt(sum(r2.^2));
%}