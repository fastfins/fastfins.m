function options = rto_hyper_options_old(varargin)
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
defaultLkd  = 'Gauss';
expectedLkd = {'Gauss', 'Poisson'};
%
defaultLambdaA = 1;
defaultDeltaA  = 1;
defaultLambdaT = 1E-4;
defaultDeltaT  = 1E-4;
%
defaultGammaL  = 1E-5;
defaultGammaR  = 1E1;
defaultGammaP  = 4;
%
p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
%
addOptional (p,'options',  defaultOptions, @(x) isstruct(x));
addParameter(p,'lkd_type', defaultLkd,  @(x) any(validatestring(x,expectedLkd)));
%
addParameter(p,'lambda_a', defaultLambdaA, validScalarPosNum);
addParameter(p,'lambda_t', defaultLambdaT, validScalarPosNum);
addParameter(p,'delta_a',  defaultDeltaA,  validScalarPosNum);
addParameter(p,'delta_t',  defaultDeltaT,  validScalarPosNum);
addParameter(p,'gamma_l',  defaultGammaL,  validScalarPosNum);
addParameter(p,'gamma_r',  defaultGammaR,  validScalarPosNum);
addParameter(p,'gamma_p',  defaultGammaP,  validScalarPosNum);
%
p.KeepUnmatched = false;
parse(p,varargin{:});
tmp     = p.Results;
options = tmp.options;

tmp     = rmfield(tmp, 'options');

if isempty(options)
    options = tmp;
else
    list = {'lkd_type', 'lambda_a', 'lambda_t', 'delta_a', 'delta_t', ...
        'gamma_l', 'gamma_r', 'gamma_p'};
    default_list = cellstr(p.UsingDefaults);
    for i = 1:length(list)
        if ~ismember(list{i}, default_list)
            options.(list{i}) = tmp.(list{i});
        end
    end
end

end
