function applyOpticalFlow(pivStack, im0, Options)
%applyOpticalFlow(pivStack, im0, Options)
% Given an input image im0, advect the  intensity along the pivStack using
% Euler integration and output the deformed image
%
% Parameters
% ----------
% Options : struct with fields
%   mandatory fields:
%       outDir
%       minT
%       maxT
%   optional fields:
%       plot_texturepaths : bool (default=true)
%           make movie of optical flow
%       plot_scatterpaths : bool (default=false)
%           make movie of tracer beads
%       clipImage : bool (default=true)
%           clip the texturepatch image to a rectangle 
%       szX : size of resampled images in which PIV is performed
%       szY : size of resampled images in which PIV is performed
%       X0 : evaluation grid points, x coord
%       Y0 : evaluation grid points, y coord
%
% NPMitchell 2021

% Default options
plot_texturepaths = true ;
plot_scatterpaths = false ;
clipImage = true ;

outDir = Options.outDir ;
mint = Options.minT ;
maxt = Options.maxT ;

% resized dimensions of piv grid --> nearly (szX_orig, szY_orig) * isf
szX = 696 ;
szY = 820 ;

EdgeLength = 15;
% rescaleFactor = 0.4;
% szX_orig = 1738 ;
% szY_orig = 2050 ;

% in resized pixels
[X0,Y0] = meshgrid(EdgeLength/2:EdgeLength:(szX-EdgeLength/2), ...
    EdgeLength/2:EdgeLength:(szY-EdgeLength/2)); 


% Unpack options
if isfield(Options, 'plot_texturepaths')
    plot_texturepaths = Options.plot_texturepaths ;
end
if isfield(Options, 'plot_scatterpaths')
    plot_scatterpaths = Options.plot_scatterpaths ;
end
if isfield(Options, 'szX')
    szX = Options.szX ;
end
if isfield(Options, 'szY')
    szX = Options.szY ;
end
if isfield(Options, 'X0')
    X0 = Options.X0 ;
end
if isfield(Options, 'Y0')
    Y0 = Options.Y0 ;
end
if isfield(Options, 'clipImage')
    clipImage = Options.clipImage ;
end

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
options.clipY = false ;
[YY, XX] = pullbackPathlines(pivStack, Y0, X0, 0, options) ;

% Visualize the result
if plot_scatterpaths
    if ~exist(fullfile(outDir, 'scatter'), 'dir')
        mkdir(fullfile(outDir, 'scatter'))
    end
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
        outfn = fullfile(outDir, 'scatter', sprintf('scat_%06d.png', tidx)) ;
        saveas(gcf, outfn)
        pause(0.1)
    end
end

%% Texture patch
if plot_texturepaths
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
        Options.scaleData = false ;
        if clipImage
            pbV(:, 1) = min(pbV(:, 1), szX) ;
            pbV(:, 2) = min(pbV(:, 2), szY) ;
        end
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

        disp(['Writing image: ' outfn])
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

