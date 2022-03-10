function timeStampStripe7(da, genotype, label, Options)
% TIMESTAMPSTRIPE7(da, genotype, label, Options)
%   Find appropriate timestamp for data stained by 'label' for all fixed
%   pullbacks compared to master timeline using leading edge of stripe7.
%   This is basically a wrapper for queriedSample's method timeStamp.m
%   
% Parameters
% ----------
% da : dynamicAtlas class instance
%   the dynamicAtlas for which we find timestamps
% genotype : str
% stain : str
% Options : struct with optional fields
%   save_fancy : optional bool, default=true
%   overwrite : optional bool, default=false
%       overwrite previous results
%   hands_on : optional bool, default=true
%       use interactive domain selection for chisquare fitting
%   cdf_min : optional float, default=0.01
%       intensity cumulative distribution function minimum 
%       cutoff
%   cdf_max : optional float, default=0.999 
%       intensity cumulative distribution function maximum 
%       cutoff
%   sigma : optional float, default=20 
%       smoothing used in the stripe ID in iLastik
%   optimize_trans : optional bool, default=true
%       allow translation of the stripe7 curve in fitting to reference
%       curves
%   timelineLabel : optional str, default=<same as label>
%       allows the <label> fixed data to  be compared to 
%       For ex, label=Runt (fixed data is Runt), but
%       masterTimeLineLabel=Eve (live data is Eve).
%   timelineLeadingTrailing : optional str, default='leading'
%       whether to use leading edge of timeLine stripe7 data or trailing
%       edge
%
% Returns
% -------
%
% Outputs
% -------
% dynamicAtlas.path/genotype/label/embryoID/timematch_<label>_<timelineLabel>stripe7_chisq.mat
% dynamicAtlas.path/genotype/label/embryoID/timematch_<label>_<timelineLabel>stripe7_chisq.txt
% dynamicAtlas.path/genotype/label/embryoID/timematch_<label>_<timelineLabel>stripe7_ssr.mat
% dynamicAtlas.path/genotype/label/embryoID/timematch_<label>_<timelineLabel>stripe7_ssr.txt
%
% NPMitchell 2020

% Unpack Options for what to perform
save_fancy = true; 
overwrite = false ;
preview = false ;
overwriteStripe = false ;
hands_on = true ;
cdf_min = 0.01 ;
cdf_max = 0.999 ;
% sigma = smoothing used in the stripe ID in iLastik
sigma = 20 ;
optimize_trans = true ;
% default label to use against timeline 
timelineLabel = label ;
% default label to use against timeline 
timelineLeadingTrailing = 'leading' ;

if hands_on
    disp(['TimeStamping fixed ' genotype ' ' label ' data with inspection'])
else
    disp(['TimeStamping fixed ' genotype ' ' label ' data with automatic params'])
end

% Unpack Options
if nargin > 3
    if isfield(Options, 'save_fancy')
        save_fancy = Options.save_fancy ;
    end
    if isfield(Options, 'overwrite')
        overwrite = Options.overwrite ;
    end
    if isfield(Options, 'preview')
        preview = Options.preview ;
    end
    if isfield(Options, 'overwriteStripe')
        overwriteStripe = Options.overwriteStripe ;
    end
    if isfield(Options, 'hands_on')
        hands_on = Options.hands_on ;
    end
    if isfield(Options, 'cdf_min')
        cdf_min = Options.cdf_min ;
    end
    if isfield(Options, 'cdf_max')
        cdf_max = Options.cdf_max ;
    end
    if isfield(Options, 'sigma')
        sigma = Options.sigma ;
    end
    if isfield(Options, 'optimize_trans')
        optimize_trans = Options.optimize_trans ;
    end
    if isfield(Options, 'timeLineLabel')
        timelineLabel = Options.timelineLabel ;
    end
    if isfield(Options, 'optimize_trans')
        timelineLeadingTrailing = Options.timelineLeadingTrailing ;
    end

end

% Unpack more options
prepend = da.lookup(genotype).prepend ;
exten = da.lookup(genotype).exten ;
if strcmp(label, 'Runt')
    timerfn = 'timematch_RuntNanobody_stripe7.mat' ;
elseif strcmp(label, 'Eve')
    timerfn = 'timematch_EveYFP_stripe7.mat' ;
else
    error('Label is not Runt or Eve. What stripe7 is this? Code for it here.')
end
step = 1 ;
sigmastep = sprintf('_sigma%03d_step%03d', sigma, step) ;
dt = 1 ;

% Build reference Directory for where the timeline is stored (with the
% timeline's stripe 7 curves)
refDir = dir(fullfile(da.path, 'timing', genotype, label, [da.timeLineMethod 'corr_ss*'])) ;
refDir = fullfile(refDir(1).folder, refDir(1).name) ;


%% Build the lookuptable for just non-dynamic data (assumes dynamic data) 
% If we are comparing Runt to Runt, then compare only fixed label=Runt to 
% live timelineLabel=Runt, so here just grab the fixed data. 
% If, however, we are comparing label=Runt to timelineLabel=Eve (for ex), 
% then grab ALL data of label (Runt).
if strcmp(label, timelineLabel)
    qs = da.findStaticGenotypeLabel(genotype, label) ;
else
    qs = da.findGenotypeLabel(genotype, label) ;
end
% note: similar to lum = da.lookup(genotype).map(label) ;


opts = struct() ;
opts.method = 'stripe7' ;
opts.save_fancy = save_fancy ; 
opts.overwrite = overwrite ;
opts.preview = preview ;
opts.overwriteStripe = overwriteStripe ;
opts.hands_on = hands_on ;
opts.cdf_min = cdf_min ;
opts.cdf_max = cdf_max ;
opts.optimize_trans = optimize_trans ;
opts.timelineLabel = timelineLabel ;
opts.timelineLeadingTrailing = timelineLeadingTrailing ;
opts.prepend = prepend ;
opts.exten = exten ;
opts.timerfn = timerfn ;
opts.sigmastep = sigmastep ;
opts.dt = dt ;
opts.refDir = refDir ;
qs.timeStamp(da, genotype, label, opts)
