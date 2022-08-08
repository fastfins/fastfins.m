function [f, g, HI] = rto_obj(model, obs, rto_lin, uperp, epsilon, vr)
% compute 0.5*||T(v) - epsilon||^2, its gradient and data for computing
% matrix vector product of the hessian
%
% Tiangang Cui, August, 1, 2017

u           = rto_lin.LPhi*vr + uperp;
HI          = forward_solve(model, u);
misfit      = (HI.d - obs.data)./obs.std;   % G(v) - d
HI.mllkd    = 0.5*sum(misfit(:).^2);

tmp         = rto_lin.linref'*misfit(:);
Tvr         = rto_lin.d(:).*(tmp + vr);
res         = Tvr - epsilon;

f           = 0.5*sum(res(:).^2);

HI.mlrto    = 0.5*sum(Tvr(:).^2);
HI.res      = HI.d - obs.data;

if model.explicitJ
    Ju      = explicit_J(model, HI);
    Jvr     = Ju*rto_lin.LPhi./obs.std;
    HI.A    = rto_lin.linref'*Jvr;
    HI.J    = rto_lin.d(:).*( HI.A + eye(rto_lin.nrank) );
    g       = HI.J'*res;
else
    g       = rto_adjoint(model, obs, rto_lin, HI, res, false);
end

end