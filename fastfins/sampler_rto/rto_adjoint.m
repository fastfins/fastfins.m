function Ttr = rto_adjoint(model, obs, rto_lin, HI, res, trust_flag)
%Compute the RTO adjoint operation. (dTdx)' * y
%
% Tiangang Cui, 20, Nov, 2019

tmp     = rto_lin.d(:).*res;
a       = rto_lin.linref*tmp;

m       = size(res, 2);
gvr     = zeros(rto_lin.nrank, m);
for i = 1:m
    if trust_flag
        d   = HI.diag(:).*a(:,i);
        gvr(:,i)    = rto_lin.LPhi'*(matvec_Jty(model, HI, d./obs.std)) ...
            + rto_lin.linref'*(a(:,i) - d);
    else
        gvr(:,i)    = rto_lin.LPhi'*(matvec_Jty(model, HI, (a(:,i)./obs.std)));
    end
end

Ttr     = gvr + tmp;

%{
tmp     = rto_lin.d(:).*res;
a       = rto_lin.Psi*(rto_lin.s(:).*tmp);

m       = size(res, 2);
gvr     = zeros(rto_lin.nrank, m);
for i = 1:m
    gvr(:,i)    = rto_lin.LPhi'*(matvec_Jty(model, HI, (a(:,i)./obs.std)));
end

Ttr     = gvr + tmp;
%}

end