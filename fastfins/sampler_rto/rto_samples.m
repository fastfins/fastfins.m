function out = rto_samples(model, obs, prior, rto_lin, rs)
%GET_MAP_MATLAB
%
% Runs optimization algorithms to get the MAP estitimate
%
% Tiangang Cui, August, 1, 2017

n = size(rs, 2);

epsilons    = rto_lin.Phi'*rs;
vperps      = rs - rto_lin.Phi*epsilons;
uperps      = matvec_prior_L (prior, vperps) + repmat(prior.mean_u, 1, n);

vr_init     = rto_lin.Phi'*rto_lin.v;
%vr_init     = zeros(size(epsilon,1),1);
% notify-detailed


opt_HM = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,...
    'HessianMultiplyFcn',@(HI,dv) rto_hessmult(model, obs, rto_lin, HI, dv, false), 'Display','off', ...
    'MaxIterations', 500, 'OptimalityTolerance', 1E-8, 'FunctionTolerance', 1E-10, 'StepTolerance', 1E-8);

out.logrtos = zeros(1, n);
out.logds = zeros(1, n);
out.llkds = zeros(1, n);
out.qois  = zeros(1, n);
out.res   = zeros(obs.n_data, n);
out.ni  = zeros(1, n);
out.nf  = zeros(1, n);
vrs     = zeros(rto_lin.nrank, n);

for i = 1:n
    [vrnew,fnew,~,output,~,HI] = fminunc_2023a(@(vr) rto_obj(model, obs, rto_lin, uperps(:,i), epsilons(:,i), vr), vr_init, opt_HM);
    
    %disp(i)
    if abs(fnew) > 1E-4
        disp(i)
        disp(fnew)
    end
    
    out.logds(i) = rto_func_d(model, obs, rto_lin, HI, false);
    out.llkds(i) = -HI.mllkd;
    out.logrtos(i) = -HI.mlrto;
    vrs(:,i)     = vrnew;
    out.res(:,i) = HI.res(:);
    out.qois(i)  = HI.qoi;
    out.ni(i) = output.iterations;
    out.nf(i) = output.funcCount;
end

out.vrs     = vrs;
out.vs      = rto_lin.Phi*vrs + vperps;
out.us      = rto_lin.LPhi*vrs + uperps;

%l2res       = 0.5*sum(rs.^2, 1);
out.lps     = -0.5*sum(out.vs.^2, 1);

% new.l2res - new.mlp has cancellations, as vperp is the same
% new.l2res - new.mlp = 0.5*( sum(epsilon.^2) - sum(vr.^2));
% abs(l2res + out.lps - 0.5*( sum(epsilons.^2, 1) - sum(vrs.^2, 1) ))

out.logrtos = out.logrtos + out.logds + rto_lin.logc;
out.logrhos = -0.5*sum(vrs.^2, 1) + out.llkds - out.logrtos;

% for exact forward model, out.logrhos gives the log of omega, where 
% omega = post/rto
% for an approximate model, one can use out.logrhos - out.llkds + exact_llkd(out.vs)
% to get the omega

end

