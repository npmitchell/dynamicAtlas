function [i_tau0j_tau0jrelaxed, itau0jInfo, NLKLBLInfo] = ...
    makeMasterTimeLinePIV(da, genotype, label, Options)
% MAKEMASTERTIMELINEPIV(genotype, label, Options)
%   Build a master timeline for this genotype based on the dynamic
%   pullbacks contained in genotype/label/ using PIV correlation
%   
% Parameters
% ----------
% da : dynamicAtlas class instance
% genotype : str
% stain : str
% Options : struct with optional fields
%   preview
%   overwrite
%
% Returns
% -------
%
% Outputs
% -------
% dynamicAtlas.path/timing/genotype/label/<method>corr_ss%02d/
% dynamicAtlas.path/timing/genotype/label/stripe7corr_ss%02d/
% dynamicAtlas.path/timing/genotype/label/timeline_ss%02d_<method>corr/
% dynamicAtlas.path/genotype/label/embryoID/timematch_curve7_dynamic.mat
% dynamicAtlas.path/genotype/label/embryoID/timematch_curve7_dynamic.txt
%
% NPMitchell 2022

%% 0. Preparations
% General options
preview = false ;           % display intermediate results
overwrite = false ;         % overwrite previous results

%default values for display
log_display = 0;
curves_only = false ;
piv_computation_method = 'pivlab' ;
corr_method = 'pullbackPathlineDisplacement' ;

%sets so that underscore method doesn't mess up titles
set(0,'DefaultTextInterpreter','none')

% save_snaps = false;
% Unpack Options
if nargin > 3
    if isfield(Options, 'preview')
        preview = Options.preview ;
    end
    if isfield(Options, 'overwrite')
        overwrite = Options.overwrite ;
    end
    
    if isfield(Options, 'save_images')
        save_images = Options.save_images ;
    end
    
    if isfield(Options, 'piv_computation_method')
        piv_computation_method = Options.piv_computation_method ;
    end
    
    if isfield(Options, 'piv_extension')
        piv_extension = Options.piv_extension ;
    end
    
    if isfield(Options, 'curves_only')
        curves_only = Options.curves_only ;
    end
    
    if isfield(Options, 'log_display')
        log_display = Options.log_display ;
    end
   
end

outdir = fullfile(da.path, 'timing', genotype, label) ;
corrOutDir = fullfile(outdir, sprintf([corr_method 'Corr'])) ;

% build expts, a list of dynamic experiments to use for building master
% timeline
% First get the index of the requested genotype in the dynamicAtlas
genoIndex = 1 ;
found = false ;
while ~found
    if strcmp(da.genotypes{genoIndex}, genotype)
        found = true ;
    else
        genoIndex = genoIndex + 1;
    end
end

% Now pull all dynamic embryos from this stain
thismap = da.lookup(genotype) ;
methods(thismap)
dynamic_embryos = thismap.findDynamicLabel(label).meta ;

% Collate dynamic experiments into expts varible. 
for ii = 1:length(dynamic_embryos.embryoIDs)
    expts{ii} = fullfile(dynamic_embryos.folders{ii}) ;
    exptIDs{ii} = fullfile(dynamic_embryos.embryoIDs{ii}) ;
end

% Determine the master timeline experiment and name its index 'hard'
% First search for the master timeline experiment via indicator txt file
hard = 0 ;
for ii = 1:length(expts)
    if ~isempty(dir(fullfile(expts{ii}, 'master_timeline_designee.txt'))) && hard == 0
        hard = ii ;
        disp(['Found master timeline designation in embryo: ' dynamic_embryos.embryoIDs{ii}])
    end
end
% Handle case where no experiment was named master timeline
if hard == 0
    disp('No dynamic experiment has yet been designated as master timeline')
    disp("This designation signifies that which experiment's dt between frames remains fixed in the master timeline")
    disp('Choose which experiment receives this designation:')
    for ii = 1:length(expts)
        disp([exptIDs{ii} ': ' dynamic_embryos.embryoIDs{ii}])
    end
    
    % Request user to fix one master timeline from list of available
    hard_is_set = false ;
    while ~hard_is_set
        hard = input('Enter choice index: ') ;
        if isnumeric(hard)
            hard_is_set = true ;
        else
            disp('Choice entered is non-numeric!')
        end
    end
    
    % Now save designation in txt file
    fntmp = fullfile(expts{hard}, 'master_timeline_designee.txt') ;
    fid = fopen(fntmp, 'wt');
    % make header
    fprintf(fid, 'this embryo is designated as the master timeline'); 
    fclose(fid);
end

% build some strings for I/O
mipfn = [thismap.prepend thismap.exten] ;
pivfn = ['*.mat'];
extn = ['_' corr_method] ;

% Aesthetics
colorset = define_colors() ;
blue = colorset(1, :) ;
red = colorset(2, :) ;
yellow = colorset(3, :) ;
purple = colorset(4, :) ;
green = colorset(5, :) ;
sky = colorset(6, :) ;
colors = [yellow; sky] ;
gray = [0.5, 0.5, 0.5] ;


%% now get embryoIDs to iterate over
qs = da.findGenotypeLabel(genotype, label ) ;
exptIDs = qs.meta.embryoIDs  ;

% Ok, let's clear the queried Sample from memory now that we know which
% embryos are in the desired set.
clearvars qs 

% For each embryo, let's load its PIV and deform along the PIVpathlines
XYs = cell(length(exptIDs), 1) ;
t0Vs = zeros(length(exptIDs), 1) ;
for ii = 1:length(exptIDs)
    disp(['qq = ' num2str(ii)])
    eID = exptIDs{ii} ;
    qs = da.findEmbryo(eID) ;
    
    %
    options = struct() ;
    options.method = 'pivlab' ;
    qs.ensurePIV(options)
    piv = qs.getPIV(options) ;    
    piv.x = piv.X0 ;
    piv.y = piv.Y0 ;
    piv.vx = piv.vx{1} ;
    piv.vy = piv.vy{1} ;
    t0V = dlmread(fullfile(qs.meta.folders{1}, 't0V.txt')) ;
    options = struct() ;
    options.preview = true ;
    
    if ii == 1
        X0 = piv.X0 / piv.rescaleFactor ;
        Y0 = piv.Y0 / piv.rescaleFactor ; 
    else
        try
            assert(all(all(piv.x / piv.rescaleFactor == X0)))
            assert(all(all(piv.y / piv.rescaleFactor == Y0)))
        catch
            error(['X0 and Y0 are not equal across datasets indicated '...
                'for building master timeline! Be sure PIV was '...
                'computed correctly for ' eID])
        end
    end
    [XX, YY] = qs.getPullbackPathlines(piv, piv.x, piv.y, t0V, options) ;
    pbPathlines = cat(4, XX, YY) ;
    XYs{ii} = pbPathlines ;
    t0Vs(ii) = t0V ;
    
    
    % preview results
    % if ~exist(fullfile(qs.meta.folders{1}, 'pullbackPathlines_images'), 'dir')
    %     mkdir(fullfile(qs.meta.folders{1}, 'pullbackPathlines_images'))
    % end
    % 
    % for tidx = 1:size(XX, 1)
    %     fnout = fullfile(qs.meta.folders{1}, 'pullbackPathlines_images', ...
    %         sprintf('pullbackPathlines_%06d.png', tidx)) ;
    % 
    %     if ~exist(fnout, 'file')
    %         xtmp = squeeze(XX(tidx, :, :)) ;
    %         ytmp = squeeze(YY(tidx, :, :)) ;
    %         scatter(xtmp(:), ytmp(:), 5)
    %         axis equal
    %         title(['t = ' num2str(tidx)])
    %         xlim([0, 1740])
    %         zlim([0, 2050])
    %         axis off
    %         FF = getframe(gca) ;
    %         FF = FF.cdata ;
    %         imwrite(FF, fnout)
    %     end
    % end
    disp('next...')
end

%% IV. Extract cross correlations between experiments
corrImOutDir = fullfile(corrOutDir, 'correspondence_images') ;
corrPathOutDir = fullfile(corrOutDir, 'correspondence_paths') ;
log_display = false ;
extn = ['_' corr_method] ;
corrDatOutDir = fullfile(corrOutDir, 'crosscorrelation') ;
if ~exist(corrDatOutDir, 'dir')
    mkdir(corrDatOutDir)
end
if ~exist(corrPathOutDir, 'dir')
    mkdir(corrPathOutDir)
end
if ~exist(corrImOutDir, 'dir')
    mkdir(corrImOutDir)
end

for ii = 1:length(exptIDs)
    disp(['qq = ' num2str(ii)])
    for jj = 1:length(exptIDs)
        disp(['pp = ' num2str(jj)])
        
        % Define the correlation matrix filename
        ijstr = [ '_' exptIDs{ii} '_' exptIDs{jj} extn ] ;
        
        
        cpathfn = fullfile(corrPathOutDir, ['cpath' ijstr '.mat']) ;
        disp(['Seeking cpathfn = ' cpathfn])        
        
        cfn = fullfile(corrDatOutDir, ['corr' ijstr '.mat']) ;
        disp(['and seeking cfn = ' cfn])
            
        % Decide to compute the correlations or not
        if ii == jj && (~exist(cpathfn, 'file') || overwrite)                
            
            % Self-correspondence
            tpath = [(1:size(XYs{ii}, 1))', (1:size(XYs{ii}, 1))'] ;
            time_correspondences = tpath;
            save(cpathfn, 'time_correspondences') ;
        elseif ~exist(cpathfn, 'file') || overwrite || ~exist(cfn, 'file')

            if exist(cfn, 'file') && ~overwrite
                load(cfn, 'corrM', 'slopeM')
            else
                disp('computing cross-correlation matrix')
                corrM = zeros(size(XYs{ii}, 1), size(XYs{jj}, 1)) ;
                slopeM = zeros(size(XYs{ii}, 1), size(XYs{jj}, 1)) ;
                for tidxQ = 1:size(XYs{ii}, 1)
                    for tidxP = 1:size(XYs{jj}, 1)

                        % If this is X0 and Y0, we will get NaNs unless we mask
                        if tidxQ - t0Vs(ii) == 0 
                            if tidxP - t0Vs(jj) == 0
                                corrM(tidxQ, tidxP) = 1 ;
                                slopeM(tidxQ, tidxP) = 1 ;
                            else
                                corrM(tidxQ, tidxP) = 0 ;
                                slopeM(tidxQ, tidxP) = 0 ;
                            end
                        elseif tidxP - t0Vs(jj) == 0
                            corrM(tidxQ, tidxP) = 0 ;
                            slopeM(tidxQ, tidxP) = 0 ;
                        else
                            % displacement for eIDq 
                            qx = squeeze(XYs{ii}(tidxQ, :, :, 1)) - X0 ;
                            qy = squeeze(XYs{ii}(tidxQ, :, :, 2)) - Y0 ;

                            % displacement for eIDp
                            px = squeeze(XYs{jj}(tidxP, :, :, 1)) - X0 ;
                            py = squeeze(XYs{jj}(tidxP, :, :, 2)) - Y0 ;
                            corrM(tidxQ, tidxP) = corr([qx(:); qy(:)], [px(:); py(:)]) ;
                            mm = polyfit([qx(:); qy(:)], [px(:); py(:)], 1) ;
                            slopeM(tidxQ, tidxP) = sin(atan2(mm(1),1) * 2) ;
                            if mod(tidxQ - t0Vs(ii), 30)== 0 && mod(tidxP - t0Vs(jj), 30) == 0
                                clf
                                plot(qx(:), px(:), '.'); hold on; plot(qy(:), py(:), '.')
                                title(['m = ', num2str(slopeM(tidxQ,tidxP))])
                                pause(0.01)
                            end
                        end
                    end
                end
                
                save(cfn, 'corrM', 'slopeM')
            end

            % Now fast march through the matrix
            clf
            imagesc(1:size(XYs{ii}, 1), 1:size(XYs{jj}, 1), corrM .* slopeM)
            title([exptIDs{ii} ' to ' exptIDs{jj}])
            xlabel([exptIDs{ii} ' timeline [min]'])
            ylabel([exptIDs{jj} ' timeline [min]'])
            colorbar
            caxis([-1,1])
            colormap(cbkry(256))
            pause(0.01)

            %% CORRESPONDENCE PAIR
            if ~exist(cpathfn, 'file')
                cij = corrM .* slopeM;

                % Output is one float for each ROW in the image
                % imagesc(cij), or each COLUMN of imagesc(cij')
                options = struct() ;
                options.clims = [-1,1] ;
                [tpath, tpath0, cij_exponent, Woffset] = shortestPathInImage(cij, options) ;
                
                msg = 'keep this path and use for timeline? [Y/n]' ;
                title(msg)
                reply = input(msg, 's') ;
                if contains(lower(reply), 'n')
                    abortPath = true ;
                else
                    abortPath = false ;
                end
                % If the resulting correspondence path is good, save path and
                % move on. 
                if ~abortPath
                    
                    [uncertainties, uncs, movmeanWindowForUnc] = ...
                        estimateCorrespondencePathUncertainties(tpath, ...
                        tpath0, exptIDs{ii}, exptIDs{jj}, movmeanWindowForUnc) ;
                    title(['Correspondences between ' exptIDs{ii} ' and ' exptIDs{jj} ])
                    figfn = fullfile(corrImOutDir, ['correspondence' ijstr '_' corr_method '.pdf']) ;
                    disp(['Saving ' figfn])
                    saveas(gcf, figfn) ;
                    
                    % Save Woffset and cij_exponent        
                    cpath_param_fn = fullfile(corrPathOutDir, ['cpath_params' ijstr '.mat']) ;
                    save(cpath_param_fn, 'Woffset', 'cij_exponent', 'movmeanWindowForUnc')

                    % Save the path
                    time_correspondences = tpath;
                    readme = ['time_correspondences(i, j) gives the correspondence between exptID=', ...
                        exptIDs{ii} ' and ' exptIDs{jj} ', with (i,j) being integer time index in timeline ',...
                        exptIDs{ii} ' and float timestamp in timeline ' exptIDs{jj} '. ', ...
                        'Uncertainties are the time-smoothed |differences| between the maxima in cij (guesses) ', ...
                        'and the correspondence path, with the temporal smoothing of the normed differences smoothed over ',...
                        'movmeanWindowForUnc timesteps.'] ;
                    save(cpathfn, 'time_correspondences', 'uncertainties', 'uncs', 'tpath0', 'movmeanWindowForUnc', 'readme') ;
                end
                
            end
        end
        
    end
end
disp('Done with correspondence curves')


%% only continue if not only drawing correspondence curves but doing timeline as well
if (curves_only == 0)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HERE expts is rebuilt from disk if network of timeline correspondences 
    % was previously saved
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% VII. Minimize difference between timing points
    % {tj} = {dtj} + frame
    % Fit line of best fit to each ridge extraction using time-time-corr cell
    
    timelineDir = fullfile(outdir, sprintf('timeline_pivDisplacementCorr')) ;
    if ~exist(timelineDir, 'dir')
        mkdir(timelineDir)
    end
    cpathFnBase = fullfile(corrPathOutDir, ['cpath_%s_%s' extn '.mat']) ;
    
    options = struct('overwrite', overwrite, ...
        'useOffDiagonalsOnly', false) ;
    [ttc, expts, exptIDs] = ...
        buildTimeTimeCorrespondences(expts, exptIDs, hard, ...
        timelineDir, cpathFnBase, dynamic_embryos, options) ;

    %% VIII. Relax timepoints to reference curve (time of dataset #hard)
    % This is assigning an intial guess to the the entire timeline for tau0.
    % This will be updated in the next sections so its OK that its just an
    % intial guess and only includes for example TTC{1}{hard} and not TTC
    % {hard}{1}
    % so-called 'hard reference' is the master timeline
    options = struct() ;
    options.overwrite = false ;
    options.maxDistReference = 20 ;
    options.timelineStiffness = 6 ;
    options.saveResults = false ;
    featureMatchedStr = ['_' corr_method] ;
    [i_tau0j_tau0jrelaxed, itau0jInfo, NLKLBLInfo] =...
        relaxPairwiseCorrespondenceNetwork(ttc, hard, ...
        expts, exptIDs, timelineDir, featureMatchedStr, options) ;
end


%% Check consistency: look at t0Vs for these timepoints
tauShort = i_tau0j_tau0jrelaxed ;
for eIDind = 1:length(exptIDs)
    eID = exptIDs{eIDind} ;
    qs = da.findEmbryo(eID );
    t0V = dlmread(fullfile(qs.meta.folders{1}, 't0V.txt')) ;
    tstamp = tauShort(tauShort(:, 1) == eIDind, 2:3) ;
    tstamp = tstamp(t0V, :) ;
    if eIDind == 1
        t0Vs = tstamp ;
    else
        t0Vs = cat(1, t0Vs, tstamp);
    end
end

t0Vs_preRelax = t0Vs(:, 1) ;
t0Vs_relaxed = t0Vs(:, 2) ;
t0Vmean = mean(t0Vs_relaxed) ;
t0Vstd = std(t0Vs_relaxed) ;
save(fullfile(timelineDir, 't0Vs.mat'), 't0Vs_preRelax', 't0Vs_relaxed', ...
    't0Vmean', 't0Vstd')

close all
figure('Units', 'centimeters', 'position', [0,0,6,6])
h1 = histogram(t0Vs_preRelax, length(exptIDs)) ;
hold on;
h2 = histogram(t0Vs_relaxed, length(exptIDs)) ;
legend({'all-to-one', 'all-to-all'})
xlabel('t0V after timeline construction')
ylabel('occurrence')
title(['$t_0^V =' num2str(t0Vmean) '\pm' num2str(t0Vstd) '$'],...
    'interpreter', 'latex')
saveas(gcf, fullfile(timelineDir, 't0Vs.pdf'))
pause(1)


%% XIIb. Save timestamps in the original embryo folders MINUS t0V (mean)
maxv = -Inf ;
clf
figure('Units', 'centimeters', 'position', [0,0,6,6])
for ii = 1:length(expts)
    tjr = i_tau0j_tau0jrelaxed(i_tau0j_tau0jrelaxed(:, 1) == ii, 3) - t0Vmean;
    tjr0 = i_tau0j_tau0jrelaxed(i_tau0j_tau0jrelaxed(:, 1) == ii, 2) - t0Vmean;

    % For the moment, set uncertainty = 1 min --> todo: update this
    matchtime = tjr * 60 ;
    matchtime_unc = ones(size(tjr)) * 60 ;
    matchtime_minutes = tjr ;
    matchtime_unc_minutes = ones(size(tjr)) ;
    % Save Chisq timing as mat and txt
    fnmat = fullfile(expts{ii}, ['timematch' featureMatchedStr '.mat']) ;
    fntxt = fullfile(expts{ii}, ['timematch' featureMatchedStr '.txt']) ;
    save(fnmat, 'matchtime', 'matchtime_unc', 'matchtime_minutes', 'matchtime_unc_minutes')    

    disp(['Saving matchtime to ', fntxt])
    size(matchtime_minutes)
    size(matchtime_unc_minutes)
    dlmwrite(fntxt, cat(2, matchtime_minutes, matchtime_unc_minutes))


    %% plot it
    plot((1:length(tjr))+tjr(1), tjr, '.', 'markersize', 0.2, 'color', colorset(ii, :))
    hold on;
    plot((1:length(tjr0))+tjr0(1), tjr0, 'o', 'markersize', 1.3, 'color', colorset(ii, :))
    maxv = max(maxv, max(tjr(:))) ;
end

plot(1:maxv, 1:maxv, 'k--')
axis equal
axis square
xlabel('timestamp in reference timeline [min]')
ylabel('consensus timestamp')
saveas(gcf, fullfile(timelineDir, sprintf('relaxation_results_linearplot_subtract_t0V.pdf')))


%% Additional characterization: look at peak velocity wrt t0V
qs = da.findGenotypeLabel(genotype, label) ;
piv = qs.getPIV(struct('method', 'pivlab')) ;

map2exptID = zeros(length(exptIDs), 1) ;
for ii = 1:length(exptIDs)
    map2exptID(ii) = find(strcmpi(qs.meta.embryoIDs, exptIDs{ii})); 
end
vrms = {} ;
t0Ns = [] ;
t0N_indices = [] ;
close all
figure('Units', 'centimeters', 'position', [0,0,6,6])
for ii = 1:length(piv.vx)
    timestamps = i_tau0j_tau0jrelaxed(i_tau0j_tau0jrelaxed(:, 1) == ii, 3) - t0Vmean;
    timestamps = timestamps(1:end-1) ;
    sqspeed = (piv.vx{map2exptID(ii)}).^2 + (piv.vy{map2exptID(ii)}).^2 ;
    vv = sqrt(mean(mean(sqspeed, 2), 3)) ;
    vv = vv(1:end-1) ;
    vrms{ii} = vv ;
    timestamps_all{ii} = timestamps ;
    
    % t0_nikolas is max rms speed
    [~, t0N] = max(gradient(vv)) ;
    t0N_indices = [t0N_indices, t0N] ;
    t0Ns = [t0Ns, timestamps(t0N)] ;
    
    plot(timestamps, vv)
    hold on;
end
xlim([-20,80])
ylabel('$v_{rms}$ [pix/min]', 'interpreter', 'latex')
xlabel('morphological time relative to $t_0^V$ [a.u. or min]', 'interpreter', 'latex')
saveas(gcf, fullfile(timelineDir, 'rms_velocity_curves_in_timeline_t0V.pdf')) 

% make mean and std of t0N
t0Nmean = mean(t0Ns) ;
t0Nstd = std(t0Ns) ;

% Save the t0N data
save(fullfile(timelineDir, 't0Ns.mat'), 't0Ns', 't0N_indices', ...
    't0Nmean', 't0Nstd', 'vrms', 'exptIDs', 'timestamps_all')


clf
histogram(t0Ns, length(exptIDs))
xlabel('$t_0^N$ [max rms speed]', 'interpreter', 'latex')
ylabel('count')
title(['$t_0^N =' num2str(t0Nmean) '\pm' num2str(t0Nstd) '$'],...
    'interpreter', 'latex')
saveas(gcf, fullfile(timelineDir, 'rms_velocity_curves_timestamps_in_timeline.pdf')) 

%% Plot T0N-based timelines
clf
for ii = 1:length(piv.vx)
    timestamps = i_tau0j_tau0jrelaxed(i_tau0j_tau0jrelaxed(:, 1) == ii, 3) - t0Vmean - t0Ns(ii);
    timestamps = timestamps(1:end-1) ;
    
    % t0_nikolas is max rms speed
    plot(timestamps, vrms{ii})
    hold on;
end
ylabel('$v_{rms}$  [pix/min]', 'interpreter', 'latex')
xlabel('morphological time relative to $t_0^N$ [a.u. or min]', 'interpreter', 'latex')
xlim([-20-t0Nmean,80-t0Nmean])
saveas(gcf, fullfile(timelineDir, 'rms_velocity_curves_in_timeline_t0N.pdf')) 

%% Plot raw with t0N-based rigid time. 
clf
for ii = 1:length(piv.vx)
    tt = i_tau0j_tau0jrelaxed(i_tau0j_tau0jrelaxed(:, 1) == ii, 3) ;
    timestamps = (1:length(tt)) - t0N_indices(ii) ;
    timestamps = timestamps(1:end-1) ;
    
    % t0_nikolas is max rms speed
    plot(timestamps, vrms{ii})
    hold on;
end
ylabel('$v_{rms}$  [pix/min]', 'interpreter', 'latex')
xlabel('original time relative to $t_0^N$ [a.u. or min]', 'interpreter', 'latex')
xlim([-20-t0Nmean,80-t0Nmean])
saveas(gcf, fullfile(timelineDir, 'rms_velocity_curves_in_timeline_t0N_rigidTime.pdf')) 



