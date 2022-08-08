function options = mcmc_options(varargin)
%
%Set the options for mcmc. If no parameters passed in, it returns the 
%default options. If pass in a given option, it will modify the option 
%according to the modification parameters passed in.
%
%Parameters are:
%
%to be added
%
%Tiangang Cui, August, 2019

defaultOptions  = struct([]);
%
defaultNstep    = 1E5;
defaultSbatch   = 1;
%
defaultAdapt    = true;
defaultNbatch   = 50;
defaultRate     = inf;
defaultSigma    = -3;
defaultDt       = 1;
%
defaultProposal = 'OW_Prior';
expectedProposal  = {'OW_Prior', 'OW_Post', 'RW', 'MALA', 'Newton'};
defaultGibbs    = false;
defaultNumPM    = 3;
defaultLISflag  = false;
%


p = inputParser;
validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
validScalar       = @(x) isnumeric(x) && isscalar(x);
%
addOptional (p,'options',  defaultOptions, @(x) isstruct(x));
addParameter(p,'nstep',    defaultNstep,   validScalarPosNum);
addParameter(p,'sbatch',   defaultSbatch,  validScalarPosNum);
%
addParameter(p,'proposal', defaultProposal,@(x) any(validatestring(x,expectedProposal)));
addParameter(p,'using_gibbs', defaultGibbs,   @(x) islogical(x) && isscalar(x));
addParameter(p,'num_pm',   defaultNumPM,   validScalarPosNum);
addParameter(p,'lis_flag', defaultLISflag, @(x) islogical(x) && isscalar(x));
%
addParameter(p,'adapt',    defaultAdapt,   @(x) islogical(x) && isscalar(x));
addParameter(p,'nbatch',   defaultNbatch,  validScalarPosNum);
addParameter(p,'rate',     defaultRate,    validScalarPosNum);
addParameter(p,'sigma',    defaultSigma,   validScalar);
addParameter(p,'dt',       defaultDt,      validScalarPosNum);

%
p.KeepUnmatched = false;
parse(p,varargin{:});
tmp     = p.Results;
options = tmp.options;
tmp     = rmfield(tmp, 'options');

if isempty(options)
    options = tmp;
else
    list = {'nstep', 'sbatch', 'adapt', 'nbatch', 'rate', 'sigma', 'dt' ...
        'proposal', 'using_gibbs', 'num_pm', 'lis_flag'};
    default_list = cellstr(p.UsingDefaults);
    for i = 1:length(list)
        if ~ismember(list{i}, default_list) || ~isfield(options, list{i})
            options.(list{i}) = tmp.(list{i});
        end
    end
end

% post process
if isinf(options.rate)
    switch options.proposal
        case{'MALA', 'Newton'}
            options.rate = 0.58;
        case{'RW'}
            options.rate = 0.23;
        case{'OW_Prior', 'OW_Post'}
            options.rate = 0.3;
    end
end

end