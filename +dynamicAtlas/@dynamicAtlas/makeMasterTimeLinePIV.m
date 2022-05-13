function makeMasterTimeLinePIV(da, genotype, label, Options)
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
%   thres
%   ssfactor
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

%% -I. Declaration
corr_method = 'PIV' ; % realspace for correlation 

%% 0. Preparations
% General options
preview = false ;           % display intermediate results
overwrite = false ;         % overwrite previous results
% thres = 0.5 ;
ssfactor = 4 ;

%default values for display
log_display = 0;
curves_only = false ;
piv_computation_method = 'pivlab' ;
corr_method = 'pullbackPathlines' ;

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
    if isfield(Options, 'ssfactor')
        ssfactor = Options.ssfactor ;
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
corrOutDir = fullfile(outdir, sprintf([corr_method 'corr'])) ;

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
extn = [sprintf('_ss%02d', ssfactor) '_' corr_method] ;
substr = sprintf('_ss%02d', ssfactor) ;
mipsearchstr = da.lookup(genotype).prepend ;
ssmipfnBase = [ mipsearchstr substr ] ; 
ssmipfn = [mipsearchstr substr '.tif'];

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
                                disp('check here')
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
            colormap(bwr)
            pause(0.01)



            %% CORRESPONDENCE PAIR
            if ~exist(cpathfn, 'file')
                
                
                
                
                cij = corrM .* slopeM;
                cij_orig = cij ;

                % find the path (rough path via maxima) 
                tpath = zeros(size(cij, 1), 2) ;
                tpath(:, 1) = 1:size(cij, 1) ;
                for tq = 1:size(cij, 1)
                    [~, ind] = max(cij(tq, :)) ;
                    tpath(tq, 2) = ind ;
                end

                % Use correlation heatmap
                close all
                %displays log plot if indicated
                if (log_display == 1)
                    imagesc(log(cij))
                else
                    imagesc(cij)
                end
                hold on;
                plot(tpath(:, 2), tpath(:, 1), 'o') 
                ntpguess = min(size(cij)) ;
                xtmp = (1:ntpguess) + tpath(1, 2) ;
                ytmp = (1:ntpguess) + tpath(1, 1) ;
                plot(xtmp, ytmp, 'k--')
                xtmp = tpath(end, 2) - (0:ntpguess-1) ;
                ytmp = tpath(end, 1) - (0:ntpguess-1) ;
                plot(xtmp, ytmp, 'k--')
                axis equal
                axis tight
                msg = 'Does path look ok? Enter=yes, Backspace=no, n/Delete=No correspondence' ;
                xlabel(['time, dataset ', exptIDs{jj}])
                ylabel(['time, dataset ', exptIDs{ii}])
                title(msg)
                disp(msg)
                good_button = false ;
                abort = false ;
                while ~good_button
                    button = waitforbuttonpress() ;
                    if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                        move_on = true ;
                        good_button = true ;
                    elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                        move_on = false;
                        good_button = true ;
                        % Save this 
                        title(['Initial correspondence between ' exptIDs{ii} ' and ' exptIDs{jj} ])
                        set(gca,'fontsize', 12);
                        saveas(gcf, fullfile(corrImOutDir, ['correspondence' ijstr '_initial.png']))
                    elseif button && (strcmp(get(gcf, 'CurrentKey'), 'delete') || strcmp(get(gcf, 'CurrentKey'), 'n') )
                        abort = true ;
                        good_button = true ;
                        move_on = true ;
                    end
                end

                % Shortest path 
                Woffset = 0. ;
                cij_exponent = 1. ; 
                while ~move_on
                    close all
                    %displays log plot if indicated
                    if (log_display == 1)
                        imagesc(log(cij_orig))
                    else
                        imagesc(cij_orig)
                    end
                    axis equal
                    hold on;
                    plot(tpath(:, 2), tpath(:, 1), 'ro') 
                    plot(tpath(1, 2), tpath(1, 1), 'ks')
                    plot(tpath(end, 2), tpath(end, 1), 'k^')
                    xtmp = (1:ntpguess) + tpath(1, 2) ;
                    ytmp = (1:ntpguess) + tpath(1, 1) ;
                    plot(xtmp, ytmp, 'k--')
                    xtmp = tpath(end, 2) - (0:ntpguess-1) ;
                    ytmp = tpath(end, 1) - (0:ntpguess-1) ;
                    plot(xtmp, ytmp, 'k--')
                    msg = 'Do endpoints look ok? Enter=yes, Backspace=no/select' ;
                    title(msg)
                    disp(msg)
                    xlabel(['time, dataset ', exptIDs{jj}])
                    ylabel(['time, dataset ', exptIDs{ii}])
                    startpt = tpath(1, :) ;
                    % Guess start/endpt to be first/last tpath points
                    endpt = tpath(end, :) ;
                    button = waitforbuttonpress() ;
                    if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                        move_on = true ;
                        disp('set startpt/endpt:')
                        startpt
                        endpt
                    elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                        move_on = false;
                        msg = 'Start timept is ok? Enter=yes, Backspace=no' ;
                        scatter(startpt(2), startpt(1), 100, 'r', 'filled')
                        disp(msg)
                        title(msg)
                        axis equal
                        button = waitforbuttonpress() ;
                        if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                            disp('Great, how about end point?')
                        elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                            % New guess for startpoint is maximum along column
                            [~, startpt(1, 1)] = max(cij(:, 1)) ;
                            scatter(startpt(2), startpt(1), 100, 'r', 'filled')
                            msg = 'Ok, how about startpt now? Enter=yes, Backspace=no' ;
                            disp(msg)
                            title(msg)
                            button = waitforbuttonpress() ;
                            if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                                disp('Great, how about end point?')
                            elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                                msg = 'Find the startpoint manually by clicking' ;
                                disp(msg)
                                title(msg)
                                startpt = ginput(1) ;
                                startpt = [startpt(2) startpt(1)] ;
                            end
                        end
                        msg = 'End timept is ok? Enter=yes, Backspace=no' ;
                        scatter(startpt(2), startpt(1), 100, 'g', 'filled')
                        scatter(endpt(2), endpt(1), 100, 'r', 'filled')
                        disp(msg)
                        title(msg)
                        button = waitforbuttonpress() ;
                        if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                            disp('Great, all done.')
                        elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                            % New guess for startpoint is maximum along column
                            [~, endpt(1, 1)] = max(cij(:, end)) ;
                            endpt(1, 2) = size(cij, 2) ;
                            scatter(endpt(2), endpt(1), 100, 'r', 'filled')
                            msg = 'Ok, how about endpt now? Enter=yes, Backspace=no' ;
                            disp(msg)
                            title(msg)
                            button = waitforbuttonpress() ;
                            if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                                disp('Great, all done.')
                            elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                                msg = 'Find the endpoint manually by clicking';
                                disp(msg)
                                title(msg)
                                endpt = ginput(1) ;
                                endpt = [endpt(2) endpt(1)] ;
                            end
                        end
                    else
                        error('bad button press')
                    end

                    % get shortest path via fastmarching by Gabriel Peyre
                    options.propagation_type = 'normal';
                    options.Tmax = sum(size(cij))*1.2;
                    options.start_points = startpt';
                    clf;
                    disp('Performing FM');
                    options.reduc_factor = 1.0 ;
                    options.weight = 1.0;
                    %   'W' is the weight matrix (the highest, the slowest the front will move).
                    %   'start_points' is a 2 x k array, start_points(:,i) is the ith starting point .
                    %   'end_points' is a 2 x 1 array, it is the goal.
                    % Define speed of pixels' movement for marching
                    % Blur the image a bit to straighten lines, then find curve
                    % WW = imgaussfilt((max(cij(:)) - cij) - min(cij(:)), 1) ;
                    % WW = imgaussfilt(cij) ;
                    WW = imgaussfilt((cij - min(cij(:))).^ cij_exponent + Woffset) ;
                    % WW = WW - min(WW(:)) + 1e-2 ;
                    % WW = cij - min(cij(:)) + 1e-2;
                    % Compute distance transform
                    % didn't get slow version to work
                    % [DD,S] = perform_front_propagation_2d_slow(WW',...
                    %     [startpt(2);startpt(1)], [endpt(2); endpt(1)], 4000, []);
                    %
                    [DD,S] = perform_fast_marching(WW', [startpt(2)+0.1; startpt(1)+0.1], options) ;

                    % Check DD
                    % imagesc(DD)
                    % button = waitforbuttonpress() ;

                    % perform_fmstar_2d(WW', spt, ept, options);    
                    disp('Extracting Paths');
                    stepsize = 0.1 ;
                    thres_dist = 2 ;
                    str_options = [stepsize 10000];

                    % path extraction
                    options.str_options = str_options ;
                    options.trim_path = true ;
                    options.startpt = [startpt(2), startpt(1)] ;
                    options.thres_dist = thres_dist ;
                    str_options = options.str_options ;

                    % grad = compute_grad(DD);
                    % grad = -perform_vf_normalization(grad);
                    % Dx = squeeze(grad(:, :, 1)) ;
                    % Dy = squeeze(grad(:, :, 2)) ;

                    % figure;
                    % subplot(1, 2, 1)
                    % imagesc(squeeze(grad(:, :, 1))); title('grad(1)')
                    % colormap(bwr)
                    % caxis([-1,1])
                    % subplot(1, 2, 2)
                    % imagesc(Dy); title('Dy')
                    % colormap(bwr)
                    % caxis([-0.01,0.01])
                    % 
                    % figure;
                    % subplot(1, 2, 1)
                    % imagesc(squeeze(grad(:, :, 2))); title('grad(2)')
                    % colormap(bwr)
                    % caxis([-1,1])
                    % subplot(1, 2, 2)
                    % imagesc(Dx); title('Dx')
                    % colormap(bwr)
                    % caxis([-1,1])


                    % Compute gradient
                    [Dy, Dx] = gradient(DD) ;
                    Dx = - Dx ;
                    Dy = - Dy ;
                    normalization = (Dx.^2 + Dy.^2);
                    stationary = find(normalization < 1e-6) ;
                    normalization(stationary) = 1;
                    Dx = Dx .* (1./sqrt(normalization)) ;
                    Dy = Dy .* (1./sqrt(normalization)) ;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    disp('Show gradient')
                    imagesc(DD'); hold on;
                    quiver(Dx', Dy', 1)
                    title('DT with gradient')
                    ccaxis = caxis ; % get current caxis
                    caxis([0, min(3, ccaxis(2))])
                    cb = colorbar() ;
                    ylabel(cb, 'DT') ;
                    axis equal 
                    axis tight
                    pause(0.001)
                    % saveas(gcf, [outfnBase '_DD.png'])

                    % Plot with streamlines
                    disp('Showing cijstream')
                    [xx, yy ] = meshgrid(1:size(Dx, 1), 1:size(Dx, 2)) ;
                    xx2=xx(1:10:size(xx, 1), 1:10:size(xx, 2));
                    yy2=yy(1:10:size(yy, 1), 1:10:size(yy, 2));
                    cxyz = stream2(xx, yy, Dx', Dy', xx2', yy2') ;
                    %displays log plot if indicated
                    if (log_display == 1)
                        imagesc(log(cij_orig))
                    else
                        imagesc(cij_orig)
                    end
                    hold on;
                    quiver(xx, yy, Dy', Dx', 0)
                    title('cij with gradient')
                    plot(startpt(2), startpt(1), 'ko')
                    plot(endpt(2), endpt(1), 'ko')
                    title(['Cross-correlation for datasets ' exptIDs{ii} ' and ' exptIDs{jj}])
                    streamline(cxyz)
                    axis equal 
                    axis tight
                    pause(0.001)
                    % saveas(gcf, [outfnBase '_cijstream.png'])

                    % path extraction
                    % grad = cat(3, Dx, Dy) ;
                    % grad = perform_vf_normalization(grad) ;
                    % Dx = squeeze(grad(:, :, 1)) ;
                    % Dy = squeeze(grad(:, :, 2)) ;
                    % works but is uphill
                    % path = stream2(Dx, Dy, endpt(1),endpt(2), str_options);
                    % messing around here
                    path = stream2(xx, yy, Dx', Dy', endpt(2)-1, endpt(1)-1, str_options);
                    % works but does not connect to startpt
                    % path = stream2(Dy, Dx, endpt(1), endpt(2), str_options) ;
                    % path = stream2(-Dx', -Dy', startpt(1), startpt(2), str_options) ;
                    path = path{1} ;
                    path = [path(:, 2), path(:, 1)] ;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    disp('Showing DT')
                    figure(1)
                    imagesc(DD') ;
                    title('DT')
                    colorbar()
                    hold on;
                    plot(path(:, 2), path(:, 1), '.-', 'color', yellow)
                    plot(startpt(:, 2), startpt(:, 1), 'ko')
                    plot(endpt(:, 2), endpt(:, 1), 'ko')
                    xlabel(['time, dataset ', exptIDs{jj}])
                    ylabel(['time, dataset ', exptIDs{ii}])
                    title(['DT for datasets ' exptIDs{ii} ' and ' exptIDs{jj}])
                    axis equal 
                    axis tight
                    pause(0.001)
                    % saveas(gcf, [outfnBase '_DT.png'])

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    disp('Showing Dx')
                    figure(2);
                    imagesc(Dx') ;
                    colormap(bwr)
                    title('Dx')
                    colorbar()
                    %caxis([-1,1])
                    hold on;
                    plot(path(:, 2), path(:, 1), '.-', 'color', yellow)
                    plot(startpt(:, 2), startpt(:, 1), 'ko')
                    plot(endpt(:, 2), endpt(:, 1), 'ko')
                    xlabel(['time, dataset ', exptIDs{jj}])
                    ylabel(['time, dataset ', exptIDs{ii}])
                    title(['\partial_xW for datasets ' exptIDs{ii} ' and ' exptIDs{jj}])
                    axis equal 
                    axis tight
                    pause(0.001)
                    % saveas(gcf, [outfnBase '_Dx.png'])

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    disp('Showing Dy')
                    figure(3);
                    imagesc(Dy') ;
                    title('Dy')
                    colormap(bwr)
                    colorbar()
                    %caxis([-1,1])
                    hold on; 
                    plot(path(:, 2), path(:, 1), '.-', 'color', yellow)
                    plot(startpt(:, 2), startpt(:, 1), 'ko')
                    plot(endpt(:, 2), endpt(:, 1), 'ko')
                    xlabel(['time, dataset ', exptIDs{jj}])
                    ylabel(['time, dataset ', exptIDs{ii}])
                    title(['\partial_yW for datasets ' exptIDs{ii} ' and ' exptIDs{jj}])
                    axis equal 
                    axis tight
                    pause(0.001)
                    % saveas(gcf, [outfnBase '_Dy.png'])

                    disp('Showing cij')
                    figure(4);
                    %displays log plot if indicated
                    if (log_display == 1)
                        imagesc(log(cij_orig))
                    else
                        imagesc(cij_orig)
                    end
                    title('cij')
                    colormap(bwr)
                    colorbar()
                    caxis([-1,1])
                    hold on; 
                    plot(path(:, 2), path(:, 1), '.-', 'color', yellow)
                    plot(startpt(:, 2), startpt(:, 1), 'ko')
                    plot(endpt(:, 2), endpt(:, 1), 'ko')
                    xlabel(['time, dataset ', exptIDs{ii}])
                    ylabel(['time, dataset ', exptIDs{jj}])
                    title(['Correlation between dataset ' exptIDs{ii} ' and ' exptIDs{jj}])
                    axis equal 
                    axis tight
                    pause(0.001)
                    set(gca,'fontsize', 12);
                    % saveas(gcf, [outfnBase '_cij.png'])

                    disp('Showing streamlines')
                    xx2=xx(1:10:size(xx, 1), 1:10:size(xx, 2));
                    yy2=yy(1:10:size(yy, 1), 1:10:size(yy, 2));
                    cxyz = stream2(xx, yy, Dx', Dy', xx2', yy2') ;
                    plot(startpt(:, 2), startpt(:, 1), 'ko')
                    plot(endpt(:, 2), endpt(:, 1), 'ko')
                    figure(5);
                    %displays log plot if indicated
                    if (log_display == 1)
                        imagesc(log(cij_orig))
                    else
                        imagesc(cij_orig)
                    end 
                    hold on;
                    quiver(Dx', Dy', 0)
                    streamline(cxyz)
                    xlabel(['time, dataset ', exptIDs{ii}])
                    ylabel(['time, dataset ', exptIDs{jj}])
                    title(['Correlation between dataset ' exptIDs{ii} ' and ' exptIDs{jj}])
                    pause(1)
                    % saveas(gcf, [outfnBase '_streamlines.png'])

                    % cpath = compute_geodesic(DD, flipud(ept), opt);
                    % Clip the path near the start point
                    d2start = vecnorm(path-startpt,2,2) ;
                    path = path(d2start > thres_dist, :) ;
                    [tmp, ind] = unique(path(:, 1)) ;
                    inds = sort(ind) ;
                    path = path(inds, :) ;
                    cpath = [endpt; path; startpt] ;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Visualize resulting path
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    close all
                    %displays log plot if indicated
                    if (log_display == 1)
                        imagesc(log(cij_orig))
                    else
                        imagesc(cij_orig)
                    end
                    hold on;
                    plot(startpt(2), startpt(1), 'ro') 
                    plot(endpt(2), endpt(1), 'ks')
                    plot(cpath(:, 2), cpath(:, 1), '-') ;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Convert cpath to tpath: x timeline is integer, y timeline
                    % is float
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ntimeptsA = round(endpt(1) - startpt(1) + 1) ; 
                    xxx = round(startpt(1)):round(endpt(1)) ;
                    yyy = interp1(flipud(cpath(:, 1)), ...
                        flipud(cpath(:, 2)), xxx, 'pchip') ;
                    tpath = [xxx; yyy]' ;

                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % Continue visualization
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%
                    plot(tpath(:, 2), tpath(:, 1), 'o') ;
                    xlabel(['time, dataset ', exptIDs{jj}])
                    ylabel(['time, dataset ', exptIDs{ii}])
                    cb = colorbar() ;
                    ylabel(cb, 'correlation')            
                    title('Path ok? Enter=yes, backspace=no/redo')
                    button = waitforbuttonpress() ;
                    if button && strcmp(get(gcf, 'CurrentKey'), 'return')
                        move_on = true ;
                    elseif button && strcmp(get(gcf, 'CurrentKey'), 'backspace')
                        move_on = false;
                        disp('starting over with path detection')
                        Woffset = input('Set weight offset = ') ;
                        disp(' --> new Woffset')
                        if isempty(Woffset)
                            Woffset = input('Set weight offset = ') ;
                        end

                        cij_exponent = input('Set cij potential exponent = ') ;
                        disp(' --> new cij exponent')
                        if isempty(cij_exponent)
                            cij_exponent = input('Set cij potential exponent = ') ;
                        end
                    end
                end

                % If the resulting correspondence path is good, save path and
                % move on. 
                if ~abort

                    % Save Woffset and cij_exponent        
                    cpath_param_fn = fullfile(corrPathOutDir, ['cpath_params' ijstr '.mat']) ;
                    save(cpath_param_fn, 'Woffset', 'cij_exponent')

                    % Save image
                    figfn = fullfile(corrImOutDir, ['correspondence' ijstr '_' corr_method '.png']) ;
                    disp(['Saving ' figfn])
                    title(['Correspondences between ' exptIDs{ii} ' and ' exptIDs{jj} ', ' corr_method ])
                    axis equal
                    axis tight
                    set(gca,'fontsize', 12);
                    saveas(gcf, figfn) ;

                    % Save the path
                    time_correspondences = tpath;
                    save(cpathfn, 'time_correspondences') ;
                end
                
                
                
                
                
                
            end
        end
        
    end
end
disp('Done with correspondence curves')


%% only continue if not only drawing curves but doing timeline as well
if (curves_only == 0)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% HERE expts is rebuilt from disk if network of timeline correspondences 
    % was previously saved
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% VII. Minimize difference between timing points
    % {tj} = {dtj} + frame
    % Fit line of best fit to each ridge extraction using time-time-corr cell

    timelineDir = fullfile(outdir, sprintf('timeline_ss%02d_realspacecorr', ssfactor)) ;
    if ~exist(timelineDir, 'dir')
        mkdir(timelineDir)
    end

    % Build ttc: time-time correlation cell
    ttcfn = fullfile(timelineDir, 'ttc.mat') ;
    if exist(ttcfn, 'file') && ~overwrite
        tmp = load(ttcfn, 'ttc', 'exptIDs', 'hard', 'hard_reference_ID') ;
        ttc = tmp.ttc ;

        % Ensure that the exptIDs corresponding to the indices i of i_tau0j are
        % loaded from the mat file itau0jfn
        if ~isfield(tmp, 'exptIDs') || ~isfield(tmp, 'hard')
            % this is bad: the exptIDs were not saved with itau0jfn, so we
            % cannot for certain know that the indices i in i_tau0j correctly
            % index exptIDs. The solution is to erase the file on disk which is
            % missing exptIDs so we can recreate it.
            disp('exptIDs or hard reference index were not saved with itau0jfn!')
            msg = ['We must recreate it!', ...
                'Please erase ', ttcfn, ' from disk now and rerun'];
            error(msg)
        end

        exptIDs = tmp.exptIDs ;
        hard = tmp.hard ;
        hard_reference_ID = tmp.hard_reference_ID ;

        % Rebuild expts in RAM in case it is different from the saved one
        [filepath,~] = fileparts(dynamic_embryos.folders{ii}) ;
        for ii = 1:length(exptIDs)
            expts{ii} = fullfile(filepath, exptIDs{ii}) ;
        end
    else
        ttc = cell(length(expts), 1) ;
        for ii = 1:length(expts)
            % create the cell
            ttc{ii} = cell(length(expts), 1) ;
        end
        % Populate the correspondences into cell from mat files
        for ii = 1:length(expts)
            disp(['dataset ii = ', exptIDs{ii}])
            for jj = 1:length(expts)
                % Define the correlation matrix filename
                ijstr = [ '_' exptIDs{ii} '_' exptIDs{jj} extn ] ;
                cfn = fullfile(corrOutDir, ['corr' ijstr '.mat']) ;
                disp(['Seeking cfn = ' cfn])
                cpathfn = fullfile(corrPathOutDir, ['cpath' ijstr '.mat']) ;
                if exist(cpathfn, 'file')
                    load(cpathfn, 'time_correspondences');
                    ttc{ii}{jj} = time_correspondences ;
                end
            end
        end

        % Save the time-time corr cell, along with the exptIDs it indexes
        hard_reference_ID = exptIDs{hard} ;
        save(ttcfn, 'ttc', 'exptIDs', 'hard', 'hard_reference_ID')
    end

    % Plot the dynamics
    close all
    for use_offset = [true false]
        offset = 0 ;
        for ii = 1:length(ttc)
            subplot(round(length(ttc)*0.5)+1, 2, ii)
            % Do earlier dataset correspondences
            for jj = 1:ii-1
                if ~isempty(ttc{jj}{ii})
                    % displace vertically to match
                    if use_offset
                        i0 = ttc{jj}{ii}(1, 2) ;
                        j0 = ttc{jj}{ii}(1, 1) ;
                        offset = i0 - j0 ;
                    end
                    plot(ttc{jj}{ii}(:, 2), ...
                        ttc{jj}{ii}(:, 1) + offset,...
                        '.-', 'color', colorset(jj, :))
                    hold on;
                end
            end

            % Do this dataset to itself
            for jj = ii 
                plot(ttc{ii}{jj}(:, 1), ttc{ii}{jj}(:, 2), '.', ...
                    'color', colorset(jj, :))
                hold on;
            end

            % Do later dataset correspondences
            for jj = ii+1:length(ttc)
                if ~isempty(ttc{ii}{jj})
                    % displace vertically to match
                    if use_offset
                        i0 = ttc{ii}{jj}(1, 1) ;
                        j0 = ttc{ii}{jj}(1, 2) ;
                        offset = i0 - j0 ;
                    end
                    plot(ttc{ii}{jj}(:, 1), ...
                        ttc{ii}{jj}(:, 2) + offset,...
                        '.-', 'color', colorset(jj, :))
                    hold on;
                end
            end

            % labels
            ylabel(['$\tau(t_' exptIDs{ii} ')$' ], 'Interpreter', 'Latex')
            if ii == 5 || ii == 6
                xlabel('time, $t_i$', 'Interpreter', 'latex')
            end
            xlim([1, 150])
            ylim([1, 150]) 
        end

        % Legend
        subplot(round(length(ttc)*0.5)+1, 2, [length(ttc) + 1, length(ttc)+2])
        for ii = 1:length(ttc)
            plot([-1],[-1],'.-', 'color', colorset(ii, :))
            hold on
        end
        legend({'1','2','3','4','5','6'}, 'orientation', 'horizontal')
        xlim([0, 1])
        ylim([0, 1])
        axis off
        set(gcf, 'Units', 'centimeters');
        set(gcf, 'Position', [0 0 16 30]) ;
        set(gca,'fontsize', 12);
        version_of_matlab = version ;
        if use_offset
            figurefn = fullfile(timelineDir, 'time_correspondences_offset.png') ;
        else
            figurefn = fullfile(timelineDir, 'time_correspondences.png') ;
        end
        if contains(version_of_matlab, '2020')
            exportgraphics(gcf, figurefn)
        else     
            saveas(gcf, figurefn)
        end

        close all
    end

    %% VIII. Relax timepoints to reference curve (time of dataset #hard)
    % This is assigning an intial guess to the the entire timeline for tau0.
    % This will be updated in the next sections so its OK that its just an
    % intial guess and only includes for example TTC{1}{hard} and not TTC
    % {hard}{1}
    % so-called 'hard reference' is the master timeline
    to_relax = 1:length(ttc) ;
    to_relax = setdiff(to_relax, hard) ;

    % Fit other curves to the hard reference timeline
    % Tau0 is a map from timeline ti to reference timeline
    % to use, do yy = ppval(pp, xq);
    for ii = 1:length(ttc)
        if ii <= hard
            pp = spline(ttc{ii}{hard}(:, 1), ttc{ii}{hard}(:, 2)) ;
        elseif ii > hard
            pp = spline(ttc{hard}{ii}(:, 2), ttc{hard}{ii}(:, 1)) ;
        end
        if ii == 1
            tau0 = [pp] ;
        else
            tau0 = [tau0, pp] ;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Build i_tau0j master timeline lookup table
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Note that a map from indices to exptIDs must exist on disk too
    % Check for existing result
    itau0jfn = fullfile(timelineDir, 'i_tau0j.mat') ;
    if exist(itau0jfn, 'file') && ~overwrite
        % NOTE exptIDs in RAM is overwritten here to ensure that the indices of
        % the saved network correspond to the correct datasets
        tmp = load(itau0jfn, 'i_tau0j', 'ntp_tot', 'addt', 'ntps', ...
            'exptIDs', 'hard', 'hard_reference_ID') ;
        i_tau0j = tmp.i_tau0j ;
        ntp_tot = tmp.ntp_tot ;
        addt = tmp.addt ;
        ntps = tmp.ntps ;
        % Ensure that the exptIDs corresponding to the indices i of i_tau0j are
        % loaded from the mat file itau0jfn
        if ~isfield(tmp, 'exptIDs')
            % this is bad: the exptIDs were not saved with itau0jfn, so we
            % cannot for certain know that the indices i in i_tau0j correctly
            % index exptIDs. The solution is to erase the file on disk which is
            % missing exptIDs so we can recreate it.
            msg = ['exptIDs were not saved with itau0jfn! We must recreate it!', ...
                'Please erase ', itau0jfn, ' from disk now and rerun'];
            error(msg)
        end

        % IF exptIDs were loaded, we must ensure that they are the same as the
        % ones we already have on RAM. This would only be false if we partially
        % ran this method and overwrote ttc but not i_tau0j AND happened to
        % create exptIDs that were scrambled. This is unlikely, but let's check
        % that everything is ok. 
        if length(exptIDs) == length(tmp.exptIDs)
            for qq = 1:length(exptIDs) 
                if ~strcmp(exptIDs{qq}, tmp.exptIDs{qq})
                    error(['exptIDs in ', itau0jfn, ' does not match variable in RAM. ', ...
                        'You must erase ' itau0jfn ' on disk and try again'])
                end
            end
        else
            error(['exptIDs in ', itau0jfn, ' does not match variable in RAM. ', ...
                'You must erase ' itau0jfn ' on disk and try again'])
        end
    else
        % First get #tps in each dataset, which is stored in diagonals of ttc
        ntps = zeros(length(ttc), 1) ;
        addt = zeros(length(ttc), 1) ;
        tpid = [] ;
        for ii = 1:length(ttc)
            addt(ii) = sum(ntps) ;
            ntps(ii) = length(ttc{ii}{ii}) ;
            % Now label each timepoint in linear indexing
            if ii == 1
                tpid = cat(2, tpid, 1:ntps(ii)) ;
            else
                tpid = cat(2, tpid, addt(ii) + (1:ntps(ii))) ;
            end

            % Store all indices in one array
            tau0extra = fnxtr(tau0(ii), 2) ;
            tauv = ppval(tau0extra, (1:ntps(ii))) ;
            i_tau0j(addt(ii) + (1:ntps(ii)), :) = [ii * ones(ntps(ii), 1), tauv'] ;
        end
        ntp_tot = sum(ntps) ;

        % save data
        save(itau0jfn, 'i_tau0j', 'ntp_tot', 'addt', 'ntps', 'exptIDs', ...
            'hard', 'hard_reference_ID') ;

        % de-clutter
        clearvars tauv tau0extra tpid ii
    end
    % assigning each frame a linear index.  So ordering here is like
    % (1,...N1,N1+1, ... N1+N2, .... N1+N2+N3, etc).  This indexes these points
    % so that they can be treated like balls on a string which can then be
    % relaxed.

    %% IX. Next build NL, KL, BL of correspondences for master timeline
    % NL = neighbor list meaning which frames are connect to whom, KL = spring
    % constant list, and BL = bond list. i.e if frame 1 is associated with
    % frame 10000 where 10000 would be in another dataset, one element of bl
    % would be 1 10000
    NLKLBLfn = fullfile(timelineDir, 'NL_KL_BL.mat') ;

    if exist(NLKLBLfn, 'file') && ~overwrite
        tmp = load(NLKLBLfn, 'NL', 'KL', 'BL', 'exptIDs','hard_reference_ID') ;
        NL = tmp.NL ;
        KL = tmp.KL ;
        BL = tmp.BL ;

        if length(exptIDs) == length(tmp.exptIDs)
            for qq = 1:length(exptIDs) 
                if ~strcmp(exptIDs{qq}, tmp.exptIDs{qq})
                    error(['exptIDs in ', NLKLBLfn, ' does not match variable in RAM. ', ...
                        'You must erase ' NLKLBLfn ' on disk and try again'])
                end
            end
        else
            error(['exptIDs in ', NLKLBLfn, ' does not match variable in RAM. ', ...
                'You must erase ' NLKLBLfn ' on disk and try again'])
        end
        disp('loaded bond list BL, neighbor list NL, k list KL from disk')

        % Assert that hard reference is unchanged
        if ~strcmp(exptIDs{hard}, tmp.hard_reference_ID)
            error(['Hard reference stored with NLKLBL does not match current. ', 
                'Remove ' NLKLBLfn '  from disk and rerun'])
        end
    else
        % preallocate Nxlarge array for NL, KL, do not preallocate BL
        NL = zeros(ntp_tot, 300) ;
        KL = zeros(ntp_tot, 300) ;
        first = true ;
        for ii = 1:length(ttc)
            for jj = 1:length(ttc)
                if ii ~= jj
                    % Ensure that there are some correspondences
                    if ~isempty(ttc{ii}{jj})
                        nodei = round(ttc{ii}{jj}(:, 1)) + addt(ii);
                        nodej = round(ttc{ii}{jj}(:, 2)) + addt(jj);
                        bladd = [nodei, nodej] ;

                        % Add to BL
                        if first
                            BL = bladd ;
                            first = false ;
                        else
                            BL = cat(1, BL, bladd) ;
                        end

                        % build NL, KL with redundancy of correspondences built into KL
                        for id = 1:length(bladd)
                            pair = bladd(id, :) ;
                            nodei = pair(1) ; 
                            nodej = pair(2) ;
                            % ij
                            if ismember(nodej, NL(nodei, :))
                                ind = find(NL(nodei, :) == nodej) ;
                                assert(NL(nodei, ind) == nodej) ;
                                KL(nodei, ind) = KL(nodei, ind) + 1 ;
                            else
                                firstzero = find(NL(nodei, :)==0, 1) ; 
                                assert(~isempty(firstzero)) 
                                NL(nodei, firstzero) = nodej ; 
                                KL(nodei, firstzero) = 1 ;
                            end

                            disp(['ii ', num2str(ii)])
                            disp(['jj: ', num2str(jj)])
                            disp(['id: ', num2str(id)])
                            disp(['nodei: ', num2str(nodei)])
                            disp(['nodej: ', num2str(nodej)])
                            % ji
                            if ismember(nodei, NL(nodej, :))
                                ind = find(NL(nodej, :) == nodei) ;
                                assert(NL(nodej, ind) == nodei) ;
                                KL(nodej, ind) = KL(nodej, ind) + 1 ;
                            else
                                firstzero = find(NL(nodej, :)==0, 1) ; 
                                if isempty(firstzero)
                                   tester = 'hi'
                                end
                                assert(~isempty(firstzero)) 
                                NL(nodej, firstzero) = nodei ; 
                                KL(nodej, firstzero) = 1 ;
                            end
                        end
                    end
                end
            end
        end

        % Find first of the columns that are all zero
        colsZero = find(all(NL==0), 1) ;  
        NL = NL(:, 1:colsZero) ;
        KL = KL(:, 1:colsZero) ;

        disp('done building bond list BL, neighbor list NL, k list KL ')

        % SAVE NLKLBLfn
        save(NLKLBLfn, 'NL', 'KL', 'BL', 'exptIDs', 'hard', 'hard_reference_ID')
    end

    %% X. Build timeline network visualization (pairwise_corr_timeline_XXX.png)

    for use_BL = [true false]
        if use_BL
            disp('Visualizing network using BL')
        else
            disp('VIsualizing network using KL, NL')
        end
        for ii = 1:1:length(ttc)    
            close all
            for qq = 1:size(i_tau0j, 1)
                plot(i_tau0j(qq, 2), i_tau0j(qq, 1), '.', 'color', colorset(i_tau0j(qq, 1), :))
                hold on
            end
            % labels
            ylabel('dataset $i$', 'Interpreter', 'Latex')
            xlabel('time, $t_0$', 'Interpreter', 'latex')
            xlim([0, max(i_tau0j(:, 2)) + 1])

            if use_BL
                % Add bonds to plot using BL (non-reciprocal)
                for qq = 1:size(BL, 1)
                    if i_tau0j(BL(qq, 1), 1) == ii 
                        plot([i_tau0j(BL(qq, 1), 2), i_tau0j(BL(qq, 2), 2)], ...
                             [i_tau0j(BL(qq, 1), 1), i_tau0j(BL(qq, 2), 1)], '-', ...
                             'color', colorset(i_tau0j(BL(qq, 2), 1), :), 'linewidth', 0.01)
                    elseif i_tau0j(BL(qq, 2), 1) == ii
                        plot([i_tau0j(BL(qq, 1), 2), i_tau0j(BL(qq, 2), 2)], ...
                             [i_tau0j(BL(qq, 1), 1), i_tau0j(BL(qq, 2), 1)], '-', ...
                             'color', colorset(i_tau0j(BL(qq, 1), 1), :), 'linewidth', 0.01)
                    end
                end

                title(['Time relaxation network: dataset ' exptIDs{ii}])
                ylim([0, max(i_tau0j(:, 1)) + 1])
                set(gca,'fontsize', 12);
                saveas(gcf, fullfile(timelineDir, sprintf('pairwise_corr_timeline_BL_%03d.png', ii)))
                % disp('press any button on figure')
                pause(0.1)
            else
                % Add bonds to plot using NL,KL 
                % cycle only through the linear point indices that belong to ii
                for qq = (addt(ii)+1):(addt(ii) + ntps(ii))
                    for nn = NL(qq, :)
                        % plot only if real neighbor
                        if nn > 0
                            plot([i_tau0j(qq, 2), i_tau0j(nn, 2)], ...
                                 [i_tau0j(qq, 1), i_tau0j(nn, 1)], '-', ...
                                 'color', colorset(i_tau0j(nn, 1), :), 'linewidth', 0.01)
                        end
                    end
                end
                title(['Time relaxation network: dataset ' exptIDs{ii}])
                ylim([0, max(i_tau0j(:, 1)) + 1])
                set(gca,'fontsize', 12);
                saveas(gcf, fullfile(timelineDir, sprintf('pairwise_corr_timeline_%03d.png', ii)))
                % disp('press any button')
                pause(0.1)
            end
        end
    end
    disp('done with network visualization')
    close all

    %% XII. Relax the timeline network & visualize (i_tau0j_tau0jrelaxed.mat)

    % check for existing result
    fn = fullfile(timelineDir, 'i_tau0j_tau0jrelaxed.mat') ;
    if exist(fn, 'file') && ~overwrite
        load(fn, 'i_tau0j_tau0jrelaxed')
    else
        % Relax: fix the time nodes of the hard reference dataset, move all others
        options = optimset('PlotFcns',@optimplotfval, 'TolX', 1e-4, 'TolFun', 1e-6);
        x0 = i_tau0j(:, 2) ;
        % pop the indices of the fixed times from the array x0
        fixed_ind = find(i_tau0j(:, 1) == hard) ;
        fixed_xx = i_tau0j(fixed_ind, 2) ;
        x0(fixed_ind) = [] ;
        fun = @(x)springEnergy1D(x, BL, fixed_ind, fixed_xx);
        xf = fminsearch(fun, x0, options) ;
        % reinsert indices of fixed times
        xnew1 = xf(1:fixed_ind(1)-1) ;
        xnew2 = xf(fixed_ind(1):end) ;
        xnew = [xnew1; fixed_xx; xnew2 ] ;
        disp('done')

        % Plot the relaxed network result 
        close all
        offset = -0.25 ;
        for qq = 1:size(i_tau0j, 1)
            plot(i_tau0j(qq, 2), i_tau0j(qq, 1) + offset, ...
                'o', 'color', colorset(i_tau0j(qq, 1), :))
            hold on
        end
        for qq = 1:size(i_tau0j, 1)
            plot(xnew(qq), i_tau0j(qq, 1), '.', 'color', colorset(i_tau0j(qq, 1), :))
            hold on
        end

        % Show sliding movement 
        for qq = 1:size(i_tau0j, 1)
            plot([i_tau0j(qq, 2), xnew(qq)], i_tau0j(qq, 1) + [offset, 0], ...
                '-', 'color', colorset(i_tau0j(qq, 1), :))
            hold on
        end

        % Figure labels
        ylabel('dataset $i$', 'Interpreter', 'Latex')
        xlabel('time, $t_0$', 'Interpreter', 'latex')
        xlim([0, max(i_tau0j(:, 2)) + 1])

        % Title and formatting
        title(['Time relaxation network'])
        ylim([0, max(i_tau0j(:, 1)) + 1])
        set(gcf, 'Units', 'centimeters');
        set(gcf, 'Position', [0 0 36 16]) ;
        set(gca,'fontsize', 12);
        saveas(gcf, fullfile(timelineDir, sprintf('relaxation_results.png')))

        % Save the result
        i_tau0j_tau0jrelaxed = cat(2, i_tau0j, xnew) ;
        save(fn, 'i_tau0j_tau0jrelaxed')
    end
    disp('done')

    disp('Done with makeMasterTimeLinePIV.m')


end
