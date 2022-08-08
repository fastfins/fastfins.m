function vmap = get_map_2018a(ml_target, matvec_hessian, v_init)
%GET_MAP_MATLAB   
%
% Runs optimization algorithms to get the MAP estitimate
%
% Tiangang Cui, 04/May/2012

opt  = optimoptions('fminunc', 'Algorithm', 'trust-region', 'SpecifyObjectiveGradient',true,... 
    'HessianMultiplyFcn',@(HI,v) hessian(matvec_hessian, HI, v), ...
    'Display','iter', 'MaxIterations', 100);

vmap = fminunc_2018a(@(v) obj(ml_target, v), v_init, opt);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f, g, hessinfo] = obj(ml_target, v)
%MINUS_LOG_POST_MAP
% 
% Compute ths log posterior for optimization
%
% Tiangang Cui, 23/Mar/2013

[f, dummy, g, dummy, hessinfo] = ml_target(v);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w = hessian(matvec_hessian, HI, dv)

w = matvec_hessian(HI, dv) + dv;

end
