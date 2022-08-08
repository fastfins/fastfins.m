function options = bilinear_options(varargin)

defaultOptions  = struct([]);
%
defaultMeshSize = 40;
defaultXYRatio  = 1;
%
defaultTestCase = 'EIT';
expectedTestCase = {'EIT', 'Laplace', 'Heat', 'RD'};
defaultGMRES    = false;
defaultBeta     = 1;
%
defaultObsLoc   = [];
%
defaultEssBCFlag = true(4,1);
defaultEssBCVals = zeros(4,1);
%
defaultQOI   = @NOP;
defaultForce = @(xs) zeros(size(xs,2),1);
%
defaultObsTstart = inf;
defaultObsTfinal = inf;
defaultObsNtime  = 1;
defaultTinit     = inf;
defaultTnormal   = inf;
defaultNtransit  = 1;
defaultNnormal   = 1;
%

p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
%
addOptional (p,'options',  defaultOptions, @(x) isstruct(x));
addParameter(p,'mesh_size',defaultMeshSize, validScalarPosNum);
addParameter(p,'mesh_size_x',defaultMeshSize, validScalarPosNum);
addParameter(p,'mesh_size_y',defaultMeshSize, validScalarPosNum);
addParameter(p,'xyratio',  defaultXYRatio,  validScalarPosNum);
%
addParameter(p,'test_case',defaultTestCase,@(x) any(validatestring(x,expectedTestCase)));
addParameter(p,'gmres',    defaultGMRES,   @(x) islogical(x) && isscalar(x));
addParameter(p,'beta',     defaultBeta,    validScalarPosNum);
%
addParameter(p,'obs_locs', defaultObsLoc,  @(x) isnumeric(x));
%
addParameter(p,'ess_bc_vals', defaultEssBCVals, @(x) isnumeric(x) && length(x) == 4);
addParameter(p,'ess_bc_flag', defaultEssBCFlag, @(x) islogical(x) && length(x) == 4);
addParameter(p,'qoi_flux',    defaultQOI,  @(x) isa(x, 'function_handle'));
addParameter(p,'force_func',  defaultForce,@(x) isa(x, 'function_handle'));
%
addParameter(p,'obs_tstart',  defaultObsTstart, validScalarPosNum);
addParameter(p,'obs_tfinal',  defaultObsTfinal, validScalarPosNum);
addParameter(p,'obs_ntime',   defaultObsNtime,  validScalarPosNum);
addParameter(p,'tinit',       defaultTinit,     validScalarPosNum);
addParameter(p,'tnormal',     defaultTnormal,   validScalarPosNum);
addParameter(p,'ntransit',    defaultNtransit,  validScalarPosNum);
addParameter(p,'nnormal',     defaultNnormal,   validScalarPosNum);

%
p.KeepUnmatched = false;
parse(p,varargin{:});
tmp     = p.Results;
options = tmp.options;
tmp     = rmfield(tmp, 'options');

if isempty(options)
    options = tmp;
else
    list = {'mesh_size','mesh_size_x','mesh_size_y','xyratio','test_case',...
        'gmres','beta','obs_locs','ess_bc_vals','ess_bc_flag','qoi_flux','force_func',...
        'obs_tstart','obs_tfinal','obs_ntime','tinit','tnormal','ntransit','nnormal'};
    default_list = cellstr(p.UsingDefaults);
    for i = 1:length(list)
        if ~ismember(list{i}, default_list) || ~isfield(options, list{i})
            options.(list{i}) = tmp.(list{i});
        end
    end
end

if isempty(options.obs_locs)
    error('observation operator is not provided')
end

if strcmp(options.test_case, 'Heat')
    if isinf(options.obs_tstart) || isinf(options.obs_tfinal) || isinf(options.tinit) || isinf(options.tnormal)
        error('need to specify time steps for transit problems')
    end
end

end