function options = map_options(varargin)

defaultOptions  = struct([]);
%
%
defaultMaxIter  = 100;
defaultKKTtol   = 1E-5;
defaultJumpSizeTol = 1E-10;
defaultFvalTol  = 1E-10;
%
defaultSolver   = 'Newton_CG';
expectedSolver  = {'Newton', 'Newton_CG', 'BFGS', 'NCG'};
%
defaultTRradius = 100;
defaultLineMaxFeval = 20;
defaultLineFtol = 1E-4;
defaultLineGtol = 0.99;
defaultLineStep = 0.0001;
%
defaultBFGSmaxNum   = 100;
%
defaultCGrestart    = 50;
defaultCGforcingTol = 0.5;
defaultCGmaxIter    = 100;
defaultCGzeroTol    = 1E-2;


p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
%
addOptional (p,'options',   defaultOptions, @(x) isstruct(x));
%
addOptional (p,'max_iter',      defaultMaxIter, validScalarPosNum);
addOptional (p,'first_KKT_tol', defaultKKTtol,  validScalarPosNum);
addOptional (p,'jump_size_tol', defaultJumpSizeTol, validScalarPosNum);
addOptional (p,'fval_tol',      defaultFvalTol, validScalarPosNum);
%
addParameter(p,'solver',    defaultSolver, @(x) any(validatestring(x,expectedSolver)));
%
addOptional (p,'TR_radius', defaultTRradius,validScalarPosNum);
%
addOptional (p,'line_max_feval', defaultLineMaxFeval, validScalarPosNum);
addOptional (p,'line_ftol', defaultLineFtol, validScalarPosNum);
addOptional (p,'line_gtol', defaultLineGtol, validScalarPosNum);
addOptional (p,'line_step', defaultLineStep, validScalarPosNum);
%
addOptional (p,'bfgs_max_num', defaultBFGSmaxNum, validScalarPosNum);
%
addOptional (p,'CG_restart', defaultCGrestart,  validScalarPosNum);
addOptional (p,'CG_forcing_tol', defaultCGforcingTol, validScalarPosNum);
addOptional (p,'CG_max_iter', defaultCGmaxIter, validScalarPosNum);
addOptional (p,'CG_zero_tol', defaultCGzeroTol, validScalarPosNum);

%
p.KeepUnmatched = false;
parse(p,varargin{:});
tmp     = p.Results;
options = tmp.options;
tmp     = rmfield(tmp, 'options');

if isempty(options)
    options = tmp;
else
    list = {'max_iter', 'first_KKT_tol', 'jump_size_tol', 'fval_tol', 'solver', ...
        'TR_radius', 'line_max_feval', 'line_ftol', 'line_gtol', 'bfgs_max_num', ...
        'CG_restart', 'CG_forcing_tol', 'CG_max_iter', 'CG_zero_tol'};
    default_list = cellstr(p.UsingDefaults);
    for i = 1:length(list)
        if ~ismember(list{i}, default_list)
            options.(list{i}) = tmp.(list{i});
        end
    end
end

end

