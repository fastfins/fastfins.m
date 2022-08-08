function [kernel, stat] = build_kernel(proposal, stat, sigma)
% Compute the operators for proposals duirng the LIS update
% Tiangang Cui, 25/Mar/2013
%

stat.M = stat.sum/stat.n;
stat.C = stat.cross/(stat.n-1) - stat.M(:)*stat.M(:)' + eye(length(stat.M))*1E-5;
%

[V,D] = eig(stat.C);
[d,i] = sort(abs(diag(D)), 'descend');
V     = V(:,i);

kernel.dt = exp(sigma)/sqrt(sum(d));
switch proposal
    case {'MALA', 'RW'}
        kernel.L  = V*((d(:).^(0.5)*sqrt(2*kernel.dt)).*V');
        kernel.C  = V*((d(:)*kernel.dt).*V');
    case {'OW_Prior'}
        tmp       = ( (0.5*kernel.dt)*d + 1 ).^(-1); % eigen values of the B operator
        DB        = sqrt(2*kernel.dt)*sqrt(d).*tmp;
        DA        = (1 - 0.5*kernel.dt*d).*tmp;
        kernel.B  = V*(DB(:).*V');
        kernel.A  = V*(DA(:).*V');
    case {'OW_Post'}
        tmp       = 1/(1 + 0.5*kernel.dt);
        %DB        = sqrt(2*kernel.dt)*tmp*sqrt(d);
        %kernel.B  = scale_cols(V, DB)*V';
        kernel.a  = (1 - 0.5*kernel.dt)*tmp;
        kernel.b  = sqrt(2*kernel.dt)*tmp;
        kernel.D  = V*((d(:).^(-0.5)).*V'); % inv sqrt of the covariance
        kernel.B  = V*((sqrt(d(:))*kernel.b).*V');
end

kernel.ref = stat.ref;

end

