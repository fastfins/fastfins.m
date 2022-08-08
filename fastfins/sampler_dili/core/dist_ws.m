function d = dist_ws(P1, S1, P2, S2)
% compute the distance of the Gaussians
% Tiangang Cui, 05/Mar/2014

s1 = (S1.^2/sum(S1.^2)).^(0.25);
s2 = (S2.^2/sum(S2.^2)).^(0.25);

tmp = (P1.*s1(:)')'*(P2.*s2(:)');

d = sqrt(1 - sum(tmp(:).^2));

end

%disp(d);