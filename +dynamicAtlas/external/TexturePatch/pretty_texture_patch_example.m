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
VV = VV(:, [1 3 2]);


%% Show Texture Space Triangulation ---------------------------------------
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

%% Show the Pretty Textured Real-Space Mesh -------------------------------

Options = struct();
Options.PSize = 32;
Options.EdgeColor = 'none';
Options.FaceLighting = 'gouraud';
Options.AddLights = false;
Options.ApplyAmbientOcclusion = true;
Options.AmbientOcclusionFactor = 0.8;

texture_patch_2d( FF, VV, TF, TV, I, Options);
axis equal

% Add a non-visible patch object for shadows
hold on
ts = trisurf(triangulation(FF,VV));
ts.Visible = 'off';
hold off

axis equal
view([52 42]);

set(gca, 'Visible', 'off')
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], ...
    'Color', [254 194 194]/255 );

% Add shadows to the non-visible patch
add_shadow(ts, []);

%% TEXTURE_PATCH_3D Example ===============================================
% This example shows how to texture a mesh triangulation using a 3D texture
% image volume. See the documentation of 'texture_patch_3d.m' for a
% complete description of the input arguments

clear; close all; clc;

% Load the data (Gryllus Bimaculatus Limb)
load('testdata_3d.mat');

%% Texture-Mapped 3D Patch ------------------------------------------------
clc; close all;

% Reduce the size of the face image patch.  See documentation for details
Options = struct();
Options.PSize = 32;
Options.EdgeColor = 'none';
Options.FaceLighting = 'gouraud';
Options.AddLights = false;
Options.ApplyAmbientOcclusion = true;
Options.ScalarField = VV(:,3);
Options.ScalarColorMap = 'cool';
Options.ScalarAlpha = 0.5;
Options.SpecularStrength = 0.1;
Options.DiffuseStrength = 0.8;
Option.AmbientStrength = 0.9;

tic
texture_patch_3d( FF, -VV, TF, TV, IV, Options );
colormap bone
axis equal
view([-130 60]);
toc

% Add a non-visible patch object for shadows
hold on
ts = trisurf(triangulation(FF,-VV));
ts.Visible = 'off';
hold off

set(gca, 'Visible', 'off')
set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1], ...
    'Color', [254 194 194]/255 );

% Add shadows to the non-visible patch
add_shadow(ts, []);

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


