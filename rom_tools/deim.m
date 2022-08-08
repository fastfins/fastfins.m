function [e,B,K] = deim(U, nU, reg_factor)
% DEIM intepolation / regression 
%
% U is the orthognal basis, number of columns equals the target DEIM rank
% multiplies with the rate
%
% rate the over sampling rate (= 1 means intepolation)
%
% deim_type determines which method to use
% 
% Tiangang Cui, 17/August/2015


if nU > size(U, 2)
    sprintf('WARNING: number of DEIM sampling points should not be more than the matrix dimension');
    nU = size(U, 2);
    nP = nU;
elseif reg_factor < 1
    sprintf('WARNING: DEIM sampling factor should be >= 1\n');
    nP = nU;
else
    nP = min(ceil(nU*reg_factor), size(U, 2));
end


[e, B, K] = deim_qr(U(:,1:nP), nU, nP);
%[e, B, K] = deim_classic(U, nU, nP);


end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [e, B, K] = deim_qr(U, nU, nP)
%
% nU is the number of POD baiss
% nP is the number of intepolation points

[~,~,p] = qr(U', 'vector');

e = p(1:nP);              % DEIM points
K = pinv(U(e,1:nU));
B = U(:,1:nU)*K;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [e, B, K] = deim_classic(U, nU, nP)
%
% nU is the number of POD baiss
% nP is the number of intepolation points

e       = zeros(nP,1);          % selection index
[~,ri]  = max(abs(U(:,1)));     % extract max elements from each columns
e(1)    = ri;


out_nP          = nP;
for i = 2:nP
    ind         = 1:(i-1);
    A           = U(e(ind), ind);
    if cond(A)  > 1E20
        out_nP  = i-1;
        break;
    end
    tmp         = A\U(e(ind),i);
    res         = U(:,i) - U(:,ind)*tmp;
    
    [~,ri]      = max(abs(res));
    e(i)        = ri;
end


if out_nP < nP
    sprintf('WARNING: large cond, iteration stopped at %i, with target rank %i and target number of points %i\n', out_nP, nU, nP);
end

ind     = 1:out_nP;             % indices of DEIM points
out_nU  = min([out_nP, nU]);    % rank of DEIM basis

e = e(ind);               % DEIM points
K = pinv(U(e,1:out_nU));
B = U(:,1:out_nU)*K;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

