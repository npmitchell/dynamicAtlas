%% MASTER TIMING PIPELINE
% Master pipeline for computing lagrangian-advected patterns.
% Here we demonstrate by advecting Runt pattern along average PIV field
% built solely from sqh-mCherry and comparing to a later timepoint.
%
% NPMitchell 2021

%% First mount the server onto your machine 
% for example, on a Mac: Apple+K afp://flydrive.synology.me
% mount minimalData/ 

%% Let's clear our environment
clc
clear
close all


%% Add paths (this part can be slow)
% Add time_align_embryos directory to path so that dynamicAtlas package is
% available to use.
% tlaDir = '/Volumes/minimalData/code/';
% tlaDir = '/Users/mattlefebvre/Desktop/Code/code';
tlaDir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/dynamicAtlas/code';
cd(fullfile(tlaDir)) ;
addpath(genpath('dynamicAtlas')) ;
cd('dynamicAtlas')
addpath(genpath('+dynamicAtlas'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define atlasPath to be where the dynamicAtlas resides (the parent
% dynamicAtlas directory, not the project directory '+dynamicAtlas')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%atlasPath = '/Volumes/minimalData/Atlas_Data' ;
%atlasPath = '/Users/mattlefebvre/Desktop/WT_data_server/'
atlasPath = '/Volumes/Elements/Atlas_Data' ;

%% Build the dynamicAtlas
% Build dynamic atlas with all genotypes in the atlasPath
% da = dynamicAtlas.dynamicAtlas(atlasPath) ;
% Or choose which genotypes to include in atlas (default=all of them)
options = struct() ;
options.labels = {'sqh-mCherry'} ;
da = dynamicAtlas.dynamicAtlas(atlasPath, {'WT'}, options) ;

%% Search for dynamic data 
genotype = 'WT' ;
label = 'sqh-mCherry' ;
qs = da.findDynamicGenotypeLabel(genotype, label) ;

%% Build average flow field at each point in time
pivStack = qs.buildPIVStack(da, genotype, label) ;


function applyOpticalFlow(pivStack, embryoID, Options)

    npivx = pivStack.npivx ;
    npivy = pivStack.npivy ;
    % Note on convention: 
    % The number of unique X0 = size(PIV, 1)
    % The number of unique Y0 = size(PIV, 2)
    options = struct() ;
    timePoints = mint:maxt ;
    options.timePoints = timePoints ;
    options.Lx = szY ;
    options.Ly = szX ;
    [YY, XX] = pullbackPathlines(pivStack, Y0, X0, 0, options) ;

    % Visualize the result
    if plot_scatterpaths
        clf
        tidx2do = 1:10:length(timePoints) ;
        tidx2do = [tidx2do, setdiff(1:length(timePoints), tidx2do)] ;
        for tidx = tidx2do
            xx = squeeze(XX(tidx, :, :)) ;
            yy = squeeze(YY(tidx, :, :)) ;
            scatter(xx(:), yy(:), 10, 'filled')
            axis equal
            axis off
            title(['$t = $' num2str(timePoints(tidx)) ' min'], 'interpreter', 'latex')

            xlim([0, szX])
            ylim([0, szY])
            outfn = fullfile(outDir, sprintf('scat_%06d.png', tidx)) ;
            saveas(gcf, outfn)
            pause(0.1)
        end
    end

    %% Texture patch
    if plot_texturepaths
        gutDir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/' ;
        addpath(fullfile(gutDir, 'gut_matlab', 'mesh_handling'))
        addpath(fullfile(gutDir, 'gut_matlab', 'mesh_handling', 'rectilinearMesh'))
        addpath(fullfile(gutDir, 'TexturePatch'))
        % Reference image
        imfn = '/Volumes/Elements/Atlas_Data/WT/Runt/202001150004/MAX_Cyl1_2_000000_c1_rot_scaled_view1_ss04.tif' ;
        t0_furrow = 12 ;
        im0 = imread(imfn, t0_furrow);
        % Define a triangulation on the image
        faces = defineFacesRectilinearGrid([], npivx, npivy) ;
        % [pbX, pbY] = meshgrid(linspace(min(X0(:)), max(X0(:)), npivy)', ...
        %     linspace(min(Y0(:)), max(Y0(:)), npivx)') ;
        pbV = [X0(:), Y0(:)] ;
        [imX, imY] = meshgrid(linspace(1, size(im0, 1), npivy)', ...
            linspace(1, size(im0, 2), npivx)') ;
        imV = [imX(:), imY(:)] ;

        % check the triangulation -- periodic dimension is in x
        % trisurf(faces, pbX(:), pbY(:), 0*pbX(:))
        % view(2) ; axis equal ;
        % waitfor(gcf)

        % Test texturepatch
        % Options = struct() ;
        % imout = texture_patch_to_image(faces, pbV, faces, imV, im0, Options) ;
        % imshow(imout)
        % waitfor(gcf) 

        %% Make texture patch vertices into pathlines
        % outDir = fullfile(atlasPath, 'pullbackPathlines') ;
        outDir = fullfile('/Users/npmitchell/Desktop/tmp/WT_pullbackPathlines') ;
        if ~exist(outDir, 'dir')
            mkdir(outDir)
        end
        tidx2do = 1:10:length(timePoints) ;
        tidx2do = [tidx2do, setdiff(1:length(timePoints), tidx2do)] ;
        for tidx = tidx2do
            disp(['tidx = ' num2str(tidx)])
            % [imx, imy] = ndgrid(1:size(im0, 1), 1:size(im0, 2)) ;
            % griddedInterpolant(imx, imy, double(im0), 'natural', 'nearest')

            xtmp = squeeze(XX(tidx,:,:))' ;
            ytmp = squeeze(YY(tidx,:,:))' ;
            pbV = [xtmp(:), ytmp(:)] ;

            % CHECK 1
            check1 = false ;
            if check1
                % Show coord vertices
                subplot(2, 1, 1)
                scatter(pbV(:, 1), pbV(:, 2), 10, 'filled')
                axis equal 

                subplot(2, 1, 2)
                scatter(imV(:, 1), imV(:, 2), 10, 'filled')
                axis equal 
                waitfor(gcf)
            end

            % CHECK 2
            check2 = false ;
            if check2
                % Show coord vertices with faces
                subplot(2, 1, 1)
                trisurf(faces, pbV(:, 1), pbV(:, 2), 0*pbV(:, 1))
                view(2)
                axis equal 

                subplot(2, 1, 2)
                trisurf(faces,imV(:, 1), imV(:, 2), 0*imV(:, 2))
                view(2)
                axis equal 
                waitfor(gcf)
            end
            Options = struct() ;
            Options.baseSize = 100 ;
            imout = texture_patch_to_image(faces, pbV, faces, imV, im0, Options) ;

            outfn = fullfile(outDir, 'texture', sprintf('texture_%06d.png', tidx)) ;
            title(['$t = $' num2str(timePoints(tidx)) ' min'], 'interpreter', 'latex')

            if tidx == 1
                outsize = size(imout) ;
            else
                if any(size(imout) ~= outsize)
                    imout = imresize(imout, outsize) ;
                end
            end

            imwrite(imout, outfn)
            close all

            % MONTAGE
            % tiffpage = max(1, round(tidx + t0_furrow + mint)) ;
            % disp(['reading page ' num2str(tiffpage)])
            % imref = imread(imfn, tiffpage);
            % imshowpair(imref, imresize(imout, size(im0)), 'montage')
            % 
            % outfn = fullfile(outDir, 'montage', sprintf('montage_%06d.png', tidx)) ;
            % title(['$t = $' num2str(timePoints(tidx)) ' min'], 'interpreter', 'latex')
            % saveas(gcf, outfn)
            % close all
        end
    end
end







%% Search for dynamic data with PIV
% qs = findFlowGenotypeLabel(da, genotype, label) ;
% da.makeMasterFlowField('WT', 'Runt', Options)


