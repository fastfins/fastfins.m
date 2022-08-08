function [vmap, HI] = rto_cond_map(model, obs, prior, v_init)
%GET_MAP_MATLAB   
%
% Runs optimization algorithms to get the MAP estitimate
%
% Tiangang Cui, August, 1, 2017

% full hessian, with log transformation

opt_HM = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,... 
    'HessianMultiplyFcn',@(HI,v) matvec_hessian(model, obs, prior, HI, v), 'Display','off', ...
    'MaxIterations', 500);

[vmap,~,~,~,~,HI] = fminunc_2018a(@(v) obj(model, obs, prior, v), v_init, opt_HM);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, g, HI] = obj(model, obs, prior, v)
%MINUS_LOG_POST_MAP
% 
% Compute ths log posterior for optimization
%
% Tiangang Cui, 23/Mar/2013

u           = matvec_prior_L (prior, v) + prior.mean_u;
[x, dxdu]   = u2x( prior, u );

HI      = forward_solve(model, x);
misfit  = (HI.d - obs.data)./obs.std;
mllkd   = 0.5*sum(misfit(:).^2); % minus log-likelihood


f       = mllkd + 0.5*sum(v(:).^2);       % minus log posterior

% gx    = adjoint(model, HI, misfit);
gx      = matvec_Jty(model, HI, misfit./obs.std);
gmllkd  = matvec_prior_Lt(prior, dxdu.*gx);
g       = gmllkd + v;

HI.dxdu  = dxdu;
HI.mllkd = mllkd;
HI.mlpt  = f;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w = matvec_hessian(model, obs, prior, HI, dv)

du      = matvec_prior_L(prior, dv); % transform to u
temp    = repmat(HI.dxdu,1,size(dv,2));
dx      = temp.*du;

wx      = zeros(size(dv));
for i = 1:size(dv,2)
     % transform to physical space x
    tmp     = matvec_Jx (model, HI, dx(:,i))./(obs.std.^2);
    wx(:,i) = matvec_Jty(model, HI, tmp);
end

wu  = wx.*temp; % transform the matvec to u
w   = matvec_prior_Lt(prior, wu) + dv; % transform back to v

end

