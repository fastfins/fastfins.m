function Hdv = rto_hessmult(model, obs, rto_lin, HI, dv, trust_flag)
%Hessmuly with the RTO objective function
%
%Tiangang Cui, Nov 11, 2015

if model.explicitJ
    Hdv = HI.J'*(HI.J*dv);
else
    Tv  = rto_linfor (model, obs, rto_lin, HI, dv, trust_flag);
    Hdv = rto_adjoint(model, obs, rto_lin, HI, Tv, trust_flag);
end

end