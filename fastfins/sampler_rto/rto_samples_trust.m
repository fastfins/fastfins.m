function out = rto_samples_trust(model, obs, prior, rto_lin, rs, trust_size, smooth_ratio)
%GET_MAP_MATLAB
%
% Runs optimization algorithms to get the MAP estitimate
%
% Tiangang Cui, 22 Nov 2019

n = size(rs, 2);

epsilons    = rto_lin.Phi'*rs;
vperps      = rs - rto_lin.Phi*epsilons;
uperps      = matvec_prior_L (prior, vperps) + repmat(prior.mean_u, 1, n);

vr_init     = rto_lin.Phi'*rto_lin.v;
%vr_init     = zeros(size(epsilon,1),1);
% notify-detailed


opt_HM = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,...
    'HessianMultiplyFcn',@(HI,dv) rto_hessmult(model, obs, rto_lin, HI, dv, true), 'Display','off', ...
    'MaxIterations', 500, 'OptimalityTolerance', 1E-5, 'StepTolerance', 1E-5);

out.logrtos = zeros(1, n);
out.logds = zeros(1, n);
out.llkds = zeros(1, n);
out.qois  = zeros(1, n);
out.ni  = zeros(1, n);
out.nf  = zeros(1, n);
vrs     = zeros(rto_lin.nrank, n);

for i = 1:n
    [vrnew,fnew,~,output,~,HI] = fminunc_2023a(@(vr) rto_obj_trust(model, obs, rto_lin, uperps(:,i), epsilons(:,i), vr, trust_size, smooth_ratio), vr_init, opt_HM);
    
    %disp(i)
    if abs(fnew) > 1E-4
        disp(i)
        disp(fnew)
    end
    
    out.logds(i) = rto_func_d(model, obs, rto_lin, HI, true);
    out.llkds(i) = -HI.mllkd;
    out.logrtos(i) = -HI.mlrto;
    vrs(:,i)     = vrnew;
    
    out.qois(i)  = HI.qoi;
    out.ni(i) = output.iterations;
    out.nf(i) = output.funcCount;
end

out.vrs     = vrs;
out.vs      = rto_lin.Phi*vrs + vperps;
out.us      = rto_lin.LPhi*vrs + uperps;
out.lps     = -0.5*sum(out.vs.^2, 1);

out.logrtos = out.logrtos + out.logds + rto_lin.logc;
out.logrhos = -0.5*sum(vrs.^2, 1) + out.llkds - out.logrtos;

end

