function timeStampStripe7(da, genotype, label, Options)
% TIMESTAMPSTRIPE7(da, genotype, label, Options)
%   Find appropriate timestamp for data stained by 'label' for all fixed
%   pullbacks compared to master timeline using leading edge of stripe7
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
    lum = da.findStaticGenotypeLabel(genotype, label).meta ;
else
    lum = da.findGenotypeLabel(genotype, label).meta ;
end
% note: similar to lum = da.lookup(genotype).map(label) ;

%% Plotting
[colors, names] = define_colors ;
colors = colors ./ vecnorm(colors, 2, 2) ;
yellow = colors(3, :) ;
green = colors(5, :) ;
sky = colors(6, :) ;
stripecolor = green ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assumes we have already run EnsembleStripe to make smoothed images 

%% Process probabilities in iLastik -- train on stripe7 versus not a stripe

%% Identify stripes if not done so already
disp('First, we identify stripe7 in all fixed samples to be timestamped')
stripefn = 'Runt_stripe7curve.mat' ;
probfn = [sigmastep(2:end), filesep, 'smooth', filesep, label, '_smooth', ...
    sigmastep, '_bin_Probabilities.h5'] ;

for kk = 1:length(lum.folders) 

    % Extract metadata from lookupMap lum
    filename = lum.names{kk} ;
    embryoDir = lum.folders{kk} ;
    embryoID = lum.embryoIDs{kk} ;

    % Compute stripe if not on disk
    stripefnkk = fullfile(embryoDir, stripefn) ;
    recomputed_stripe = false ;
    if ~exist(stripefnkk, 'file') || overwrite
        if ~exist(stripefnkk, 'file')
            disp(['No stripe on disk. Building ' stripefnkk])
        else
            disp(['Overwriting stripe on disk: ' stripefnkk])
        end

        pfn = fullfile(embryoDir, probfn) ;
        disp([' -> Opening ' pfn])

        try
            dat = h5read(pfn, '/exported_data') ;
        catch
            disp('No h5 Probabilities file found on disk! Creating subsampled tif for iLastik input.')
            pathtofile = fullfile(lum.folders{kk}, 'sigma020_step001', 'smooth');
            dummynamein = [label '_smooth_sigma020_step001.tif'];
            dummynameout = [label '_smooth_sigma020_step001_bin.tif'];
            im = imread(fullfile(pathtofile, dummynamein));
            im2 = imresize(im, 0.5);
            imwrite(im2,fullfile(pathtofile, dummynameout));
            disp([ 'Created ' dummynameout]);
            userin = input(' ---> Create Prob file in ilastik. Shall MATLAB wait? [yes/no]', 's') ;
            if contains(userin, 'y') || contains(userin, 'Y')
                disp('Waiting for user to create iLastik output for stripe probabilities')
                input('Press Enter when finished.', 's') ;
            else
                error('Exiting so user can train in iLastik')
            end
            % re-read the data
            dat = h5read(pfn, '/exported_data') ;
        end    
        
        % % Used to crop half of the image and detect stripe7 only in the
        % % posterior half
        % midx = round(0.5 * size(dat, 2)) ;
        % midy = round(0.5 * size(dat, 3)) ;
        % dcrop = squeeze(dat(2, midx:end, :)) ;
        % addx = midx ;
        % addy = 0 ;

        % Extract the leading curve
        maskfn = fullfile(embryoDir, 'Runt_mask.tif');  

        thres = 0.5 ;
        minmaxsz = [5e3, 1e8] ;
        curv = extractStripeEdges(dat, maskfn, thres, minmaxsz, embryoID) ;
        % depth = max(curv) - min(curv) ;
        stripe7curve = curv ; 
        stripe7curve_frac = curv ;
        stripe7curve_frac(:, 1) = double(stripe7curve_frac(:, 1)) / double(size(dat, 2)) ;
        stripe7curve_frac(:, 2) = double(stripe7curve_frac(:, 2)) / double(size(dat, 3)) ;
        %%% MAKE A NOTE OF THIS!!! 
        % added the following line so that the stripe7curv_frac data is
        % properly scaled.  This will NOT be helpful for all data but here,
        % because the binary was made on an image that cut off the anterior
        % 50% of the embryo, we have to add .5 back to the stripe 7 curve
        % fraction coordinates (stripe7curv_frac)
        % stripe7curve_frac = stripe7curve_frac + [0.5 , 0];
        save(stripefnkk, 'stripe7curve', 'stripe7curve_frac')

        % resize stripe7 curv
        curv = stripe7curve_frac ;
        recomputed_stripe = true ;
    end

    % Save fancy image
    fancyImFn = fullfile(embryoDir, [label 'stripe7_rgb_' embryoID '.png']) ;
    if save_fancy && (~exist(fancyImFn, 'file') || overwrite ) 

        % Load stripe if we didn't just compute it
        if ~recomputed_stripe
            disp(['Loading stripe7curve from ' stripefnkk])
            load(stripefnkk, 'stripe7curve_frac') ;
            curv = stripe7curve_frac ;
        end

        disp(['Building overlay image: ' fancyImFn])
        close all
        % Load image and all other images of this embryo
        estruct = da.findEmbryo(embryoID).meta ;

        nlabels = length(estruct.names) ;
        % Make RGB image
        alllabels = '' ;
        for qq = 1:nlabels
            % load this channel/label
            test = fullfile(estruct.folders{qq}, estruct.names{qq});
            badtest = '/Users/mattlefebvre/Desktop/WT_data_server/WT/Runt/201907031428/MAX_Cyl1_2_000000_c2_rot_scaled_view1.tif';
            if strcmp(test,badtest)
                hi = 1
            end
            imq = double(imread(fullfile(estruct.folders{qq}, estruct.names{qq}))) ;

            if qq == 1
                combined = zeros(size(imq, 1), size(imq, 2), 3) ;
            end

            % append to all labels for fancy image saving
            alllabels = [alllabels estruct.labels{qq} '_'] ;

            % Adjust intensity
            [f,x] = ecdf(imq(:));
            f1 = find(f>cdf_min, 1, 'first');
            f2 = find(f<cdf_max, 1, 'last');
            lim = [x(f1) x(f2)];
            imq = mat2gray(imq, double(lim));
            combined(:, :, 1) = combined(:, :, 1) + imq * colors(qq, 1) ;
            combined(:, :, 2) = combined(:, :, 2) + imq * colors(qq, 2) ;
            combined(:, :, 3) = combined(:, :, 3) + imq * colors(qq, 3) ;
        end
        % correct for oversaturation
        combined(combined > 1) = 1.0 ;
        combined = combined / max(combined(:)) ;
        wDV = size(combined, 1) ;
        wAP = size(combined, 2) ;

        % Now make the figure
        fig = figure('visible', 'off') ;
        imshow(combined)
        hold on;
        imwrite(combined, fullfile(embryoDir, [alllabels 'rgb_' embryoID '.png']))
        plot(curv(:, 1)*wAP, curv(:, 2)*wDV, 'o-', 'color', stripecolor)
        % plot(curv' + midx, (1:length(curv)) + midy - width, 'o-', 'color', yellow)
        title([embryoID ': Identification of stripe 7'])
        saveas(fig, fancyImFn)
        close all
    end   
end

%% Now match to time domain Runt Nanobody or EveYFP data
close all

% Load reference curves
tmp = load(fullfile(refDir, 'curve7stats_collapsed_filtered.mat')) ;
if strcmp(timelineLeadingTrailing, 'leading')
    LEX = tmp.LXs ;
    LES = tmp.LSs ;
    LEY = zeros(size(LEX)) ;
    for qq=1:size(LEX, 2)
        LEY(:, qq) = (1:size(LEX, 1)) / size(LEX, 1) ;
    end
elseif strcmp(timelineleadingTrailing, 'trailing')
    TEX = tmp.TXs ;
    TES = tmp.TSs ;
    TEY = zeros(size(TEX)) ;
    for qq=1:size(TEX, 2)
        TEY(:, qq) = (1:size(TEX, 1)) / size(TEX, 1) ;
    end
else
    error("Options.timelineLeadingTrailing must be 'leading' or 'trailing'")
end

% Go through each embryo and find time
for kk = 1:length(lum.folders)
    embryoDir = lum.folders{kk} ;
    embryoID = lum.embryoIDs{kk} ;
    eDirs4ID = da.findEmbryo(embryoID).meta ;

    % assert that the current expt is dynamic
    assert(lum.nTimePoints(kk) == 1)

    % Only compute if not done already or overwrite
    matfn = ['timematch_' label '_' timelineLabel 'stripe7_chisq.mat'] ;
    txtfn = ['timematch_' label '_' timelineLabel 'stripe7_chisq.txt'] ;
    matfn_ssr = ['timematch_' label '_' timelineLabel 'stripe7_ssr.mat'] ;
    txtfn_ssr = ['timematch_' label '_' timelineLabel 'stripe7_ssr.txt'] ;
    not_on_disk = ~exist(fullfile(embryoDir, matfn), 'file')  ;
    if overwrite || not_on_disk
        % Which adjusted curve
        load(fullfile(embryoDir, stripefn), 'stripe7curve_frac') ;
        curv = stripe7curve_frac ;

        % Match to dynamic data curv
        ssds = [] ; 
        tstamps = [] ;
        options = optimset('MaxIter', 25, 'TolX', 1e-2) ; 
        % Other options: ('FunctionTolerance', 1e-7, 'Display','iter','PlotFcns',@optimplotfval);

        % Consider each timepoint, compare to this curve to reference curves
        msg = [embryoID ': Considering ' label ' curve'] ;
        disp(msg)

        % Cut into segments before optimizing the placement of the curve7
        disp('Cutting curve into leading and trailing segments')
        edgeskk = cutCurveIntoLeadingTrailingSegments(curv); 
        % Each element of lekk is a cell. edgeskk{1} is the (AP, DV) coords of
        % the leading edge. edgeskk{2} is the (AP, DV) coords of the trailing
        % edge.
        lekk = edgeskk{1} ;  % leading edge (AP, DV)
        tekk = edgeskk{2} ;  % trailing edge (AP, DV)
        
        % Debug
        % figure; plot(lekk(:, 1), lekk(:, 2), '.')
        % hold on;
        % plot(LEX(:, 1))

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Chisq
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute chisquared(t) for this curve
        disp('Computing chisquared')
        if strcmp(timelineLeadingTrailing, 'leading')
            [chisq, chisqn, ~, ssr] = chisquareCurves(lekk, LEX', LEY', (LES.^2)', ...
                true, optimize_trans) ;
        else
            [chisq, chisqn, ~, ssr] = chisquareCurves(tekk, TEX', TEY', (TES.^2)', ...
                true, optimize_trans) ;
        end
        % if allow_rotation
        %     guess = [0, 0, 0]; 
        %     opt = fminsearch(@(x0)curveDifferenceRotTrans(x0, cvfit, xyd), guess, options) ;
        % else
        %     guess = [0, 0]; 
        %     opt = fminsearch(@(x0)curveDifferenceTrans(x0, cvfit, xyd), guess, options) ;
        % end

        % Fit chisq to parabola
        % Convert into match and uncertainty, using at least 4 and at most 10
        % timepoints to fit the parabolic minimum.
        chisqn = fillmissing(chisqn, 'movmedian', 5) ;
        if hands_on
            % Load the data to show the time matching
            load(fullfile(embryoDir, stripefn), 'stripe7curve_frac') ;
            if strcmp(timelineLeadingTrailing, 'leading')
                dataStruct.refcurves = cat(3, LEX, LEY) ;
            elseif strcmp(timelineLeadingTrailing, 'trailing')
                dataStruct.refcurves = cat(3, TEX, TEY) ;
            else
                error("Options.timelineLeadingTrailing must be 'leading' or 'trailing'")
            end
            dataStruct.dataframe = imread(fullfile(embryoDir, lum.names{kk})) ;
            header_preamble = [embryoID ': '] ;
            [matchtime, matchtime_unc, fit_coefs, nidx, pidx] = ...
                chisqMinUncInteractiveDomain(chisqn, 4, 10, ssr, ...
                            header_preamble, dataStruct) ;
            % Save fit figure
            figfn = fullfile(embryoDir, ...
                ['timematch_' label '_' timelineLabel 'stripe7_chisq_fit_interactivedomain.png']) ;
            saveas(gcf, figfn)
            clf
        else         
            [tmatch, unc, fit_coefs] = chisqMinUncertainty(chisqn, 4, 10) ;
        end
        matchtime_minutes = matchtime * dt ;
        matchtime_unc_minutes = matchtime_unc * dt ;

        % Plot the result from ChiSq unc
        timedense = 1:0.1:length(chisqn) ;
        % minimimum of the fit is cstar
        cstar = fit_coefs.ystar ;
        close all ;
        fig = figure('visible', 'off') ;
        plot(chisqn, '.-')
        hold on;
        yy = polyval(fit_coefs.p, timedense, fit_coefs.S, fit_coefs.mu) ;
        plot(timedense, yy, '--')
        errorbar(matchtime, cstar, matchtime_unc, 'horizontal')
        ylim([0, 30])
        xlims = xlim() ;
        xlim([0 xlims(2)])
        ylabel('\chi^2 / N')
        xlabel('time [min]')
        title([embryoID, ...
            ': a = ' num2str(matchtime_minutes), ...
            ' \pm ', num2str(matchtime_unc_minutes), '  min'])

        % Save into all embryoID matching dirs
        for qq=1:length(eDirs4ID.folders)
            edir = eDirs4ID.folders{qq} ;
            fitfn = [ 'timematch_' label '_' timelineLabel 'stripe7_chisq_fit.png'] ;
            disp(['Saving ' fitfn ' into ' edir])
            figfn = fullfile(edir, fitfn) ;
            saveas(fig, figfn)
        end
        clf

        % Save Chisq timing as mat and txt for all matching embryoID dirs
        for qq=1:length(eDirs4ID.folders)
            edir = eDirs4ID.folders{qq} ;
            disp(['Saving ' matfn ' into ' edir])

            fnmat = fullfile(edir, matfn) ;
            fntxt = fullfile(edir, txtfn) ;
            save(fnmat, 'matchtime', 'matchtime_unc', 'matchtime_minutes', 'matchtime_unc_minutes')    

            disp(['Saving matchtime to ', fntxt])
            dlmwrite(fntxt, [matchtime_minutes, matchtime_unc_minutes])
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SSR 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Filter to get minimum of SSR
        % windowWidth = 3;  % for quick smoothing
        % kernel = ones(windowWidth, 1) / windowWidth;
        % ssdsm = filter(kernel, 1, ssds);
        [~, minID] = min(ssr) ;
        optID = max(1, minID-6):min(minID+6, length(ssr)) ;
        % pentad = ssdsm(optID) ;    
        pentad = ssr(optID) ;    

        % Fit to parabola
        % plot(tstamps, ssds)
        % plot(ssds)
        % plot(tstamps(optID), ssds(optID))    
        % Fit to y = a x^2 + b x + c
        [p, S] = polyfit(optID(:), pentad(:), 2) ;
        a = p(1) ;
        b = p(2) ;
        % [y_fit,delta] = polyval(p, trange(pentx) * dt, S);
        ci = polyparci(p, S, 0.6827) ;
        a_unc = a - ci(1, 1) ;
        b_unc = b - ci(1, 2) ;
        % Find minimum
        if a > 0
            matchtime = -b / (2 * a) ;
            dmdb = -1 / (2 * a) ;
            dmda = 0.5 * b / a^2 ; 
            matchtime_unc = sqrt((dmda * a_unc)^2 + (dmdb * b_unc)^2) ;
        else
            % the quadratic fit is so bad that it is upside down
            matchtime = minID ;
            matchtime_unc = NaN ;
        end
        matchtime_minutes = matchtime * dt ;
        matchtime_unc_minutes = matchtime_unc * dt ;

        % save an image of the five-pt fit
        close all
        fig = figure('visible', 'off') ;
        % errorbar(trange(optID) * dt / 60, pentad, pentad_unc)
        hold on;
        plot(optID * dt, pentad, 's')
        ttt = min(optID * dt):0.1:max(optID * dt) ;
        tttm = ttt * dt ; % convert to minutes
        plot(tttm, polyval(p, ttt, S), '-')
        plot(matchtime_minutes, polyval(p, matchtime, S), 'o', 'color', green)
        ploterr(matchtime_minutes, polyval(p, matchtime, S), matchtime_unc_minutes, [])
        ylabel('Sum of squared distances')
        xlabel('time [min]')
        title([embryoID, ...
            ': a = ' num2str(matchtime_minutes), ...
            ' \pm ', num2str(matchtime_unc_minutes), '  min'])
        ylims = ylim;
        plot((1:length(ssr))*dt, ssr, 'k.--')
        ylim(ylims)
        fitfn_ssr = [ 'timematch_' label '_' timelineLabel 'stripe7_ssr_fit.png'] ;
        fn = fullfile(embryoDir, fitfn_ssr) ;
        disp(['Saving figure to ' fn])
        saveas(fig, fn)
        close all

        % Save timing as mat and txt
        for qq=1:length(eDirs4ID.folders)
            edir = eDirs4ID.folders{qq} ;
            disp(['Saving ' matfn_ssr ' into ' edir])
            fnmat = fullfile(edir, matfn_ssr) ;
            fntxt = fullfile(edir, txtfn_ssr) ;
            save(fnmat, 'matchtime', 'matchtime_unc', 'matchtime_minutes', 'matchtime_unc_minutes')

            disp(['Saving matchtime to ', fntxt])
            dlmwrite(fntxt, [matchtime_minutes, matchtime_unc_minutes])
        end
    else
        % Check that timestamp exists in all embryoID dirs
        disp(['Timestamp exists in ' label ' embryoID directory'])
        disp('Checking that timestamp exists in all embryoID directories')
        timestamp_loaded = false ;
        for qq=1:length(eDirs4ID.folders)
            edir = eDirs4ID.folders{qq} ;
            fnmat0 = fullfile(embryoDir, matfn) ;
            fnmat = fullfile(edir, matfn) ;
            fntxt = fullfile(edir, txtfn) ;
            if ~exist(fnmat, 'file')
                disp([' > Copying ' matfn ' into ' edir])
                % Load the timestamp if we haven't done so
                if ~timestamp_loaded
                    load(fnmat0, 'matchtime', 'matchtime_unc', ...
                        'matchtime_minutes', 'matchtime_unc_minutes')  
                end
                disp([' > > Saving ' matfn ' into ' edir])
                save(fnmat, 'matchtime', 'matchtime_unc', ...
                    'matchtime_minutes', 'matchtime_unc_minutes')
                
                % Now save it in the other embryoID directory (other label)
                save(fnmat, 'matchtime', 'matchtime_unc', 'matchtime_minutes', 'matchtime_unc_minutes')    
                disp([' >  > Copying matchtime into ', fntxt])
                dlmwrite(fntxt, [matchtime_minutes, matchtime_unc_minutes])
            end
        end
        
        disp('Checking that ssr timestamp exists in all embryoID directories')
        timestamp_ssr_loaded = false ;
        for qq=1:length(eDirs4ID.folders)
            edir = eDirs4ID.folders{qq} ;
            fnmat0_ssr = fullfile(embryoDir, matfn_ssr) ;
            fnmat_ssr = fullfile(edir, matfn_ssr) ;
            fntxt_ssr = fullfile(edir, txtfn_ssr) ;
            if ~exist(fnmat_ssr, 'file')                
                disp([' > Copying ' matfn_ssr ' into ' edir])
                % Load the timestamp if we haven't done so
                if ~timestamp_ssr_loaded
                    load(fnmat0_ssr, 'matchtime', 'matchtime_unc', ...
                        'matchtime_minutes', 'matchtime_unc_minutes')  
                end
                disp([' > > Saving timematch_stripe7_ssr.mat into ' edir])
                save(fnmat_ssr, 'matchtime', 'matchtime_unc', ...
                    'matchtime_minutes', 'matchtime_unc_minutes')
                disp([' > > Copying matchtime to ', fntxt_ssr])
                dlmwrite(fntxt_ssr, [matchtime_minutes, matchtime_unc_minutes])
            end
        end
    end    
end
