function [sqrt_eta, b] = eval_pre_data_rom(rom, ur)
%eval_pre_data
%
% pre evaluate functions at quadrature points.
%
% Tiangang Cui, 10/Sep/2020

%{
grad_se = zeros(rom.eta_ne, rom.nd);
for di = 1:rom.nd
    grad_se(:,di) = rom.eta_U{di}*ur;
end
%}
grad_se = reshape(rom.eta_U*ur, [], rom.nd);
sqrt_eta = rom.eta_K*( (sum(grad_se.^2,2)+rom.epsilon).^(rom.p_rate/4-1/2) );

pb = zeros(sum(rom.b_ne), 1);
for di = 1:rom.nd
    %{
    grad_sbi = zeros(rom.b_ne{di}, rom.nd);
    for dj = 1:rom.nd
        grad_sbi(:,dj) = rom.b_U{di,dj}*ur;
    end
    %}
    ind = (1:rom.b_ne(di)) + sum(rom.b_ne(1:(di-1)));
    grad_sbi = reshape(rom.b_U{di}*ur, [], rom.nd);
    pb(ind) = ( (sum(grad_sbi.^2,2)+rom.epsilon).^(rom.p_rate/4-1) ).*grad_sbi(:,di);
end
b = rom.b_K*pb;

end