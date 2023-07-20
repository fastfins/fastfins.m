function [f, g, HI] = rto_obj_trust(model, obs, rto_lin, uperp, epsilon, vr, trust_size, smooth_ratio)
% compute 0.5*||T(v) - epsilon||^2, its gradient and data for computing
% matrix vector product of the hessian
%
% Tiangang Cui, 20 Nov 2019

u           = rto_lin.LPhi*vr + uperp;
HI          = forward_solve(model, u);
Fnorm       = HI.d./obs.std;


misfit      = Fnorm - rto_lin.data;   % G(v) - d
HI.mllkd    = 0.5*sum(misfit(:).^2);

% Lxr = Fref + J(vr - xr)
linfor      = rto_lin.Fref + rto_lin.linref*(vr - rto_lin.xref);
% F(xr) - Lxr, normalised
fdiff       = Fnorm - linfor;

% trust region function and gradient
d           = abs(rto_lin.Fref) + 1;
[tr_f,tr_g] = trust_region_rto(fdiff./d, trust_size, smooth_ratio);
tr_f        = tr_f.*d;
HI.diag     = tr_g;

% modified forward model
ftilde      = tr_f + linfor;
mis         = ftilde - rto_lin.data;

% apply the Q^t matrix
tmp         = rto_lin.linref'*mis(:);
Tvr         = rto_lin.d(:).*(tmp + vr);
res         = Tvr - epsilon;

f           = 0.5*sum(res(:).^2);

HI.mlrto    = 0.5*sum(Tvr(:).^2);

if model.explicit_ja
    Ju      = explicit_jacobian(model, HI);
    Jvr     = Ju*rto_lin.LPhi./obs.std;
    Jnew    = HI.diag(:).*(Jvr - rto_lin.linref) + rto_lin.linref;
    
    HI.A    = rto_lin.linref'*Jnew;
    HI.J    = rto_lin.d(:).*( HI.A + eye(rto_lin.nrank) );
    g       = HI.J'*res;
else
    g       = rto_adjoint(model, obs, rto_lin, HI, res, true);
end

end
