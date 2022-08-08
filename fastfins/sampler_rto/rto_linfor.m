function Tv = rto_linfor(model, obs, rto_lin, HI, dvr, trust_flag)
%Compute the RTO linearized operation. dTdx * x
%
% Tiangang Cui, August, 1, 2017

m       = size(dvr, 2);
% du      = matvec_prior_L (prior, dv);
du      = rto_lin.LPhi*dvr;

Jv      = zeros(obs.n_data, m);
for i = 1:m
    Jv(:,i) = matvec_Ju(model, HI, du(:,i))./obs.std;
end

if trust_flag
    tmp = rto_lin.linref*dvr;
    Jv  = HI.diag(:).*(Jv - tmp) + tmp;
end

a   = rto_lin.linref'*Jv;
Tv  = rto_lin.d(:).*(a + dvr);

end

