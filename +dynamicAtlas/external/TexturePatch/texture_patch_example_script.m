%% TEXTURE_PATCH EXAMPLE SCRIPT ===========================================
% In this script we show how to use TEXTURE_PATCH_2D and TEXTURE_PATCH_3D
% to create faithfully textured mesh triangulations from 2D/3D texture
% images/image volumes.
%
% By Dillon Cislo
%==========================================================================

%% TEXTURE_PATCH_2D Example ===============================================
% This example shows how to texture a mesh triangulation using a 2D texture
% image. See the documentation of 'texture_patch_2d.m' for a complete
% description of the input arguments

clear; close all; clc;

% Load the data (Grecian Urn)
load('testdata_2d.mat');

% Show the texture-space triangulation.  Notice that there are fewer
% texture vertices than real-space vertices, but the number of texture
% space faces equals the number of real space faces.  Also notice that TV
% is given in (row,column) format
figure

imshow( I, [] );
hold on

% Re-scale texture vertex coordinates to proper pixel indicies
TVS = [ (size(I,2)-1) .* TV(:,2) + 1, (size(I,1)-1) .* TV(:,1) + 1 ];
triplot( triangulation( TF, TVS ), 'r', 'LineWidth', 2 );

hold off

clear TVS

%% Show the Textured Real-Space Mesh --------------------------------------

texture_patch_2d( FF, VV, TF, TV, I );
axis equal

%% TEXTURE_PATCH_3D Example ===============================================
% This example shows how to texture a mesh triangulation using a 3D texture
% image volume. See the documentation of 'texture_patch_3d.m' for a
% complete description of the input arguments

clear; close all; clc;

% Load the data (Gryllus Bimaculatus Limb)
load('testdata_3d.mat');

%% Simple Interpolation ---------------------------------------------------
% For comparison we first show the results of a simplistic interpolation
% scheme using the MATLAB built-in patch functionality. We find a set of
% vertex texture coordinates by indexing into the image volume.  Notice
% that TV is given in (row,column,page) format:
% i.e., TV = [ VV(:,2) VV(:,1) VV(:,3) ].
CDataIDx = sub2ind( size(IV), ...
    round(TV(:,1)), round(TV(:,2)), round(TV(:,3)) );
CData = IV( CDataIDx );

patch( 'Faces', FF, 'Vertices', VV, 'FaceVertexCData', CData, ...
    'FaceColor', 'interp', 'EdgeColor', 'none' );
colormap bone
axis equal

clear CData CDataIDx

%% Texture-Mapped 3D Patch ------------------------------------------------

% Reduce the size of the face image patch.  See documentation for details
Options.PSize = 6;
Options.EdgeColor = 'none';
tic
texture_patch_3d( FF, VV, TF, TV, IV, Options );
colormap bone
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
toc

%% Texture-Mapped 3D Patch with Rotation + Translation --------------------

% Reduce the size of the face image patch.  See documentation for details
for i = 1
    Options.PSize = 16;
    Options.EdgeColor = 'none';
    alph = i * 1.0 ;
    Options.Rotation = [[cos(alph), -sin(alph), 0]; ...
        [sin(alph), cos(alph), 0 ]; ...
        [0, 0, 1]] ;
    Options.Translation = [300 * i, 0, 0] ;

    texture_patch_3d( FF, VV, TF, TV, IV, Options );
    colormap bone
    hold on
    ttext = ['Rotated by \alpha=' num2str(alph) ] ;
    ttext = [ttext ', translated by ' num2str(Options.Translation(1))] ;
    title(ttext)
    pause(1)
end
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

