function options = lis_options(varargin)
%
%Set the options for dimension reduction. If no parameters passed in, it
%returns the default options. If pass in a given option, it will modify the
%option according to the modification parameters passed in.
%
%Parameters are:
%
%to be added
%
%Tiangang Cui, August, 2019

defaultOptions  = struct([]);
%
defaultMethod   = 'Laplace';
expectedMethod  = {'Prior', 'Laplace', 'DILI', 'MALA'};
%
defaultHessType = 'Eig';
expectedHessType  = {'Eig', 'SVD'};
%
defaultMaxHess  = 500;
defaultMinHess  = 100;
%
defaultLocalMaxDIM    = 100;
defaultLocalTruncTol  = 1E-2;
defaultLISConvTol     = 1E-2;
%
defaultGlobalTruncTol = 1E-3;
defaultGlobalConvTol  = 1E-2;
defaultGlobalMaxDim   = 500;
%
p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
%
addOptional (p,'options',  defaultOptions, @(x) isstruct(x));
addParameter(p,'method',   defaultMethod,  @(x) any(validatestring(x,expectedMethod)));
addParameter(p,'hess_type',defaultHessType,@(x) any(validatestring(x,expectedHessType)));
%
addParameter(p,'max_hess', defaultMaxHess, validScalarPosNum);
addParameter(p,'min_hess', defaultMinHess, validScalarPosNum);
addParameter(p,'local_max_dim',   defaultLocalMaxDIM,   validScalarPosNum);
addParameter(p,'local_trunc_tol', defaultLocalTruncTol, validScalarPosNum);

addParameter(p,'global_max_dim',  defaultGlobalMaxDim,  validScalarPosNum);
addParameter(p,'global_trunc_tol',defaultGlobalTruncTol,validScalarPosNum);
addParameter(p,'global_conv_tol', defaultGlobalConvTol, validScalarPosNum);
addParameter(p,'lis_trunc_tol',   defaultLISConvTol,    validScalarPosNum);
%
p.KeepUnmatched = false;
parse(p,varargin{:});
tmp     = p.Results;
options = tmp.options;
tmp     = rmfield(tmp, 'options');

if isempty(options)
    options = tmp;
else
    list = {'method', 'hess_type', 'max_hess', 'min_hess', 'local_max_dim', 'local_trunc_tol', ...
        'global_max_dim', 'global_trunc_tol', 'global_conv_tol', 'lis_trunc_tol'};
    default_list = cellstr(p.UsingDefaults);
    for i = 1:length(list)
        if ~ismember(list{i}, default_list)
            options.(list{i}) = tmp.(list{i});
        end
    end
end

end