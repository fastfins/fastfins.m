function w = bilinear_weights_rec(bl, tr, x)

hx  = tr(1) - bl(1);
hy  = tr(2) - bl(2);

s1  = ( x(1) - bl(1) ) / hx;
s2  = ( x(2) - bl(2) ) / hy;

w = [(1-s1).*(1-s2), s1.*(1-s2), s1.*s2, (1-s1).*s2];

end