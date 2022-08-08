function [f, g] = trust_region_rto(x, tau, e)


tr = @(r, tau, e) r.*(abs(r)<=tau*(1-e)) ...
    + tau*sign(r).*(abs(r)>=tau*(1+e)) + ...
    (-(r-tau).^2/(4*e*tau) + (r-tau)/2 + tau - tau*e/4).*(abs(r)>tau*(1-e) & abs(r)<tau*(1+e));


%{
r = x/tau;
f = ones(size(x));
g = ones(size(x));
     
ind1 = abs(r)<=1;
ind2 = abs(r)> 1;

if sum(ind1) > 0
    f(ind1) = r(ind1);
end

if sum(ind2) > 0
    f(ind2) = sign( r(ind2) );
    g(ind2) = 0;
end

f = f*tau;
%}


r = x/tau;
f = ones(size(x));
g = ones(size(x));

ar = -1/(4*e);
br =  1/2;
cr =  1-e/4;

al =  1/(4*e);
bl =  1/2;
cl =  e/4-1;
     
     
ind1 = abs(r)<=(1-e);
ind2 = r>(1-e) & r<(1+e);
ind3 = r<(e-1) & r>(-1-e);
ind4 = abs(r)>= (1+e);

if sum(ind1) > 0
    f(ind1) = r(ind1);
end

if sum(ind2) > 0
    t = r(ind2) - 1;
    f(ind2) = ar*t.^2 + br*t + cr;
    g(ind2) = 2*ar*t + br;
end

if sum(ind3) > 0
    t = r(ind3) + 1;
    f(ind3) = al*t.^2 + bl*t + cl;
    g(ind3) = 2*al*t + bl;
end


if sum(ind4) > 0
    f(ind4) = sign( r(ind4) );
    g(ind4) = 0;
end

f = f*tau;


end