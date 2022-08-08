function out = rto_one_sample(model, obs, prior, rto_lin)
%GET_MAP_MATLAB
%
% Runs optimization algorithms to get the MAP estitimate
%
% Tiangang Cui, August, 1, 2017

r       = randn(prior.dof, 1);
epsilon = rto_lin.Phi'*r;
vperp   = r - rto_lin.Phi*epsilon;
uperp   = matvec_prior_L (prior, vperp) + prior.mean_u;

vr_init = rto_lin.Phi'*rto_lin.v;
% notify-detailed

opt_HM = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,... 
    'HessianMultiplyFcn',@(HI,dv) rto_hessmult(model, obs, rto_lin, HI, dv), 'Display','off', ...
    'MaxIterations', 500, 'OptimalityTolerance', 1E-8, 'StepTolerance', 1E-8);

[vrnew,fnew,~,~,~,out] = fminunc_2020a(@(vr) rto_obj(model, obs, rto_lin, uperp, epsilon, vr), vr_init, opt_HM);

%out.output  = HI.d;
out.logd = rto_func_d(model, obs, rto_lin, out);
out.vr   = vrnew;

out.v    = rto_lin.Phi*vrnew + vperp;
out.u    = rto_lin.LPhi*vrnew + uperp;

%l2res       = 0.5*sum(r.^2, 1);
out.lp   = -0.5*sum(out.v.^2, 1);

% new.l2res - new.mlp has cancellations, as vperp is the same
% new.l2res - new.mlp = 0.5*( sum(epsilon.^2) - sum(vr.^2));
% abs(l2res + out.lp - 0.5*( sum(epsilon.^2, 1) - sum(vrnew.^2, 1) ))

out.logrho  = 0.5*( sum(epsilon.^2, 1) - sum(vrnew.^2, 1) ) - out.mllkd - out.logd + rto_lin.logc;

%0.5*sum(epsilon.^2, 1) 
%0.5*sum(vrnew.^2, 1) 

if abs(fnew) > 1E-6
    disp(fnew)
    %out.logrho  = - rto_lin.logd + rto_lin.logc;
end

end

