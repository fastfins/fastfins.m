function options = tomo_options(varargin)

defaultOptions = struct([]);
%
defaultMeshSize = 40;
%
defaultAngle = 2*pi/3;
defaultNSourc = 10;
defaultNDetec = 40;
defaultDWidth = 2;
%
defaultIs = 10;

defaultFType  = 'exp';
expectedFType = {'none', 'exp', 'erf'};
defaultLogThres = eps;
defaultErfScale = 0.5;
defaultErfShift = 0.5;
%

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
%
addOptional (p,'options',  defaultOptions,  @(x) isstruct(x));
addParameter(p,'mesh_size',defaultMeshSize, validScalarPosNum);
addParameter(p,'angle',    defaultAngle,    validScalarPosNum);
addParameter(p,'n_sourc',  defaultNSourc,   validScalarPosNum);
addParameter(p,'n_detec',  defaultNDetec,   validScalarPosNum);
addParameter(p,'d_width',  defaultDWidth,   validScalarPosNum);
addParameter(p,'Is',       defaultIs,       validScalarPosNum);

addParameter(p,'f_type',    defaultFType,    @(x) any(validatestring(x,expectedFType)));
addParameter(p,'log_thres', defaultLogThres, validScalarPosNum);
addParameter(p,'erf_scale', defaultErfScale, validScalarPosNum);
addParameter(p,'erf_shift', defaultErfShift, @(x) isnumeric(x) && isscalar(x));
%
%
p.KeepUnmatched = false;
parse(p,varargin{:});
tmp     = p.Results;
options = tmp.options;
tmp     = rmfield(tmp, 'options');

if isempty(options)
    options = tmp;
else
    list = {'mesh_size','angle', 'n_sourc', 'n_detec','d_width' ...
            'Is', 'f_type', 'log_thres', 'erf_scale', 'erf_shift'};
    default_list = cellstr(p.UsingDefaults);
    for i = 1:length(list)
        if ~ismember(list{i}, default_list) || ~isfield(options, list{i})
            options.(list{i}) = tmp.(list{i});
        end
    end
end

end