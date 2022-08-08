function options = rto_hyper_options(varargin)
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
defaultName  = 'name';
%
defaultType  = 'gamma';
expectedType = {'gamma', 'beta'};
%
defaultAlpha = 1;
defaultBeta  = 1E-4;
%
defaultLeft  = 0;
defaultRight = inf;
%
p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
%
addOptional (p,'options',  defaultOptions, @(x) isstruct(x));
addParameter(p,'name', defaultName);
addParameter(p,'type', defaultType,  @(x) any(validatestring(x,expectedType)));
%
addParameter(p,'alpha', defaultAlpha, validScalarPosNum);
addParameter(p,'beta',  defaultBeta,  validScalarPosNum);
addParameter(p,'left',  defaultLeft,  validScalarPosNum);
addParameter(p,'right', defaultRight, validScalarPosNum);
%
p.KeepUnmatched = false;
parse(p,varargin{:});
tmp     = p.Results;
options = tmp.options;

tmp     = rmfield(tmp, 'options');

if isempty(options)
    options = tmp;
else
    list = {'name', 'type', 'alpha', 'beta', 'left', 'right'};
    default_list = cellstr(p.UsingDefaults);
    for i = 1:length(list)
        if ~ismember(list{i}, default_list)
            options.(list{i}) = tmp.(list{i});
        end
    end
end

end
