%test script

clear all
close all
clc

%% Define QueriedSample on which to inspect timing
tlaDir = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/dynamicAtlas/code';
% tlaDir = '/mnt/data/code/';

cd(fullfile(tlaDir)) ;
addpath(genpath('dynamicAtlas')) ;
cd('dynamicAtlas')
addpath(genpath('+dynamicAtlas'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define atlasPath to be where the dynamicAtlas resides (the parent
% dynamicAtlas directory, not the project directory '+dynamicAtlas')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
atlasPath = '/Volumes/minimalData/Atlas_Data' ;
%atlasPath = '/Users/mattlefebvre/Desktop/WT_data_server/'
% atlasPath = '/Volumes/Elements/Atlas_Data' ;
% atlasPath = '/run/user/1001/gvfs/afp-volume:host=flydrive.local,user=npmitchell,volume=minimalData/Atlas_Data/' ;

% Build the dynamicAtlas
% Build dynamic atlas with all genotypes in the atlasPath
% da = dynamicAtlas.dynamicAtlas(atlasPath) ;
% Or choose which genotypes to include in atlas (default=all of them)
options = struct() ;
options.labels = {'endoECad-GFP', 'ECad-GFP'} ;
da = dynamicAtlas.dynamicAtlas(atlasPath, {'WT'}, options) ;

% Search for dynamic data 
genotype = 'WT' ;
% label = 'sqh-mCherry' ;
label = 'ECad-GFP' ;
qs = da.findDynamicGenotypeLabel(genotype, label) ;


%% 

%whether or not to display the movie
dispMov = 0;

%whether or not to save the movie
saveMov = 0;

%path to realspace stripe folder
stripe_dir = '/Volumes/minimalData/Atlas_Data/timing/WT/Runt/realspacecorr_ss04/';

%handle for the stripe itself
stripe_handle = 'curve7stats_collapsed_filtered.mat';

%path to stripe
stripe_path = fullfile(stripe_dir, stripe_handle);

%overall directory with all the pullbacks
master_dir = '/mnt/data/embryo/Marion_Atlas_All_LR-Preserved_Pullbacks/';
%cd into the master directory
cd(master_dir)

%get the directory info: enables one to iterate over ALL the genotypes
%MAKE SURE ALL FOLDERS HAVE AT LEAST ONE DASH! ENABLES EASY WAY TO AVOID . and .. IN THIS!
d_master = dir('*-*');
numGenoFolders = length(d_master);

%Genotype to look at
genoIDX = 4;

%experiment to look at 
exptIDX = 1;

%name of genotype
geno_name = d_master(genoIDX).name;
%path of genotype
geno_path = fullfile(master_dir, geno_name);

%path of genotype
cd(geno_path)
d_geno = dir('*20*');

%name of experiment
expt_name = d_geno(exptIDX).name;

%path to expt
expt_path = fullfile(geno_path, expt_name);

%handle for quantitative computations
comp_handle = 'Quantitative Computations';

%handle for name of structure to load
vort_handle = 'vorts_cc_noSym.mat';

%handle for name of piv data
piv_handle = 'piv-data.mat';

% Finding t0v 
%
% %NOTE: HAS to be mean rather than rms because need SIGNED quantity, need 
% %the time when goes from LOW to HIGH values since that is where the 
% %desired transition happens, at the higher value border of the cross!!
% 
% %maximum time length for this embryo
% max_time = size(vort_corr_map, 1);
% 
% %take the gradient
% [GX, GY] = gradient(vort_corr_map);
% 
% %value of mean of gradient columns
% mean_GX_cols = zeros(max_time, 1);
% 
% %stores values
% for c = 1 : max_time
%     mean_GX_cols(c) = mean(GX(:,c));
% end
% 
% %max index of the columns
% [~, mIdx] = max(mean_GX_cols);
% 
% %represents the old time which BECOMES the new t=1. 
% %
% %this value is one more than the maxIdx to get the start of the 
% %stationary region, the left border of this pixel, placed in the top
% %left hand corner! 
% optimal_new_t1_idx = mIdx+1;

%%

%using t0V to find stripes

%load stripe itself
stripes = load(stripe_path);

%max time
maxT = size(stripes.LX, 2);

sampleT = 21;
%start and end times
tStart = sampleT;%1;
tEnd  = sampleT;%maxT;

for t_sample = tStart : tEnd

    % %sample time to use
    % t_sample = 35;

    %parametrization coords
    xx = stripes.xx;

    %structure of all leading edge stripes
    LX = stripes.LX;

    %structure of smoothed standard deviations
    LSs = stripes.LSs;

    %structure of smoothed standard errors
    LVs = stripes.LVs;

    %structure of all trailing edge stripes
    TX = stripes.TX;

    %stripe AP values
    stripe_t = LX(:,t_sample);
    %stripe AP standard deviations
    stdev_t = LSs(:,t_sample);
    %stripe AP standard errors
    stderr_t = LVs(:,t_sample);


    %scale for AP and DV used in resized pullbacks
    scale_AP_resized = 696;
    scale_DV_resized = 820;

    %dimensions of original images
    orig_AP = 1738;
    orig_DV = 2050;

    %scale factor between resized pixels and micrometers
    %isf = 0.4;
    isf = scale_AP_resized / orig_AP;
    px2um = 0.2619;
    resized_pix_to_um = px2um / isf;

    if dispMov

        %plot the stripe

        %resized coords according to usual dimensions
        resized_xx_pix = scale_DV_resized * xx;
        resized_xx_um = resized_pix_to_um * resized_xx_pix;

        %resized stripe according to usual dimensions, multiplies by 1-stripe_t
        %to flip the orientation properly so A is on top and P is on bottom
        resized_stripe_t_pix =  scale_AP_resized * stripe_t;
        resized_stripe_t_um = resized_pix_to_um * resized_stripe_t_pix;

        %resized standard deviation
        resized_stdev_t_pix = scale_AP_resized * stdev_t;
        %resized standard error
        resized_stderr_t_pix = scale_AP_resized * stderr_t;


        %plot the stripe data
        plot(resized_xx_um, resized_stripe_t_um, '-b', 'LineWidth', 1.5)
        xlabel('Position along DV axis (um)')
        ylabel('Position along AP axis (um)')

        %reverse the active axes
        h = gca;
        set(h, 'YDir', 'reverse')

        xlim([0 resized_pix_to_um * scale_DV_resized])
        ylim([0 resized_pix_to_um * scale_AP_resized])

        title(['Runt Stripe 7 Leading Edge at time ', num2str(t_sample), ' min'])

        pause(.1)

        M(t_sample) = getframe(gcf);

    end

end

if (saveMov == 1)

    cd(geno_path)
    v = VideoWriter(strcat('runt_stripe7_leading_edge_mov_program'), 'Uncompressed AVI');
    open(v)
    writeVideo(v, M)
    close(v)

end


%%
% % 
% figure(2)
% 
% %plot(resized_xx_pix, resized_stripe_t_pix, '-b')
% 
% errorbar(resized_xx_pix, resized_stripe_t_pix, resized_stdev_t_pix, '-b')
% 
% xlabel('Position along DV axis (resized pix)')
% ylabel('Position along AP axis (resized pix)')
% 
% %reverse the active axes
% h = gca;
% set(h, 'YDir', 'reverse')
% 
% xlim([0 scale_DV_resized])
% ylim([0 scale_AP_resized])

%%

clear M_sample

%path to t0V stripe structure
path_to_t0V_stripe = '/mnt/data/embryo/Marion_Atlas_All_LR-Preserved_Pullbacks/Runt-nanobody/stripe_t0V.mat';

%loads in the t0V stripe (will be in resized pixels, same as piv
s = load(path_to_t0V_stripe);
stripe_t0V = s.stripe_t0V;

%The standard grid used for resized piv
EdgeLength = 20;
[X_grid,Y_grid] = meshgrid(EdgeLength/2:EdgeLength:scale_DV_resized-EdgeLength/2,EdgeLength/2:EdgeLength:scale_AP_resized-EdgeLength/2); 

%IMPORTANT: MODIFICATION SO THAT GRID IS ORIENTED AT BOTTOM LEFT, NOT TOP LEFT!
Y_grid = flipud(Y_grid);

%Genotype to look at
genoIDX_sample = 2;

%experiment to look at 
exptIDX_sample = 1;


%name of genotype
geno_name_sample = d_master(genoIDX_sample).name;
%path of genotype
geno_path_sample = fullfile(master_dir, geno_name_sample);

%path of genotype
cd(geno_path_sample)
d_geno_sample = dir('*20*');

%name of experiment
expt_name_sample = d_geno_sample(exptIDX_sample).name;

%path to expt
expt_path_sample = fullfile(geno_path_sample, expt_name_sample);

%handle for stripe structure to save
stripeSaveHandle = 'PIV-Advected_Runt_stripe7curve.mat';

%handle for text file containing t0V, ASSUME JUST CONTAINS ONE NUMBER!
t0V_handle = 't0V.txt';

%IMPORTANT: the t0V for this PARTICULAR EMBRYO
t0V_SAMPLE = load(fullfile(expt_path_sample, 't0V.txt'));

%loads in the piv data
p = load(fullfile(expt_path_sample, comp_handle, piv_handle));
piv_data_sample = p.piv_data;

%number of timepoints in the sample (will be one MORE than the PIV!)
numT_sample = length(piv_data_sample) + 1;

%gets the initial stripe pointcloud at t0V from the stripe
pc_t0V_x = stripe_t0V.xx;
pc_t0V_y = stripe_t0V.stripe_pix_t0V;
pc_t0V = [pc_t0V_x, pc_t0V_y];

%to store the temporary advected pointclouds
pc_t = pc_t0V;

%number of adjusted timepoints
num_adj_T = numT_sample-(t0V_SAMPLE-1);

%cell ary of all pointclouds, stored for examining advected stripe
pc_ary = cell(num_adj_T, 1);

%IMPORTANT: cell ary of all pointclouds, stored TO SAVE
pc_ary_TOSAVE = cell(num_adj_T, 1);

%IMPORANT: scale factor to use to save
%e.g. isf was 0.4 so divided orig dimensions by 2.5 for piv-advected stripes
%AND runt stripes subsampled by 4, so divided orig dimensions by 4
%need to save the piv-advected as runt stripes format, so multiply by 2.5/4
scalefac_TOSAVE = 2.5/4;

%uses the piv to advect the original pointcloud over time
for t = t0V_SAMPLE : numT_sample
    
    %stores the current pointcloud in the appropriate time for advectio
    pc_ary{t-(t0V_SAMPLE-1)} = pc_t;
    
    %IMPORTANT: switches the order, and scales by appropriate factor,
    %to the structure to save. Also reverses AP direction because of the
    %conventions used for saving
    pc_ary_TOSAVE{t-(t0V_SAMPLE-1)} = scalefac_TOSAVE * [scale_AP_resized - pc_t(:,2), pc_t(:,1)];
    
    %gets the stripe coordinates at this time
    pc_t_X = pc_t(:,1);
    pc_t_Y = pc_t(:,2);
    
    
    %plots the current stripe
    
    figure(3)

    plot(pc_t_X, pc_t_Y, '-b')

    %errorbar(resized_xx_pix, resized_stripe_t_pix, resized_stderr_t_pix, '-b')

    xlabel('Position along DV axis (resized pix)')
    ylabel('Position along AP axis (resized pix)')
    xlim([0 scale_DV_resized])
    ylim([0 scale_AP_resized])
    title({['Advected stripe7 pointcloud for ', geno_name_sample, '-', expt_name_sample], [' at time ', num2str(t), ' min (t0V = ', num2str(t0V_SAMPLE), 'min)']})
    
    hold on
    
%     q_scale = 1;
%     quiver(X_grid, Y_grid, q_scale * Interp_t_VX(X_grid, Y_grid), q_scale * Interp_t_VY(X_grid, Y_grid), 0, 'k');
%     
    hold off
    
    pause(.1)
    M_sample(t-(t0V_SAMPLE-1)) = getframe(gcf);
    
    %DOES NOTHING OTHER THAN ADD THE POINTCLOUD AND PLOT IF REACH THE LAST TIME
    if t == numT_sample
        continue;
    end
    
    %gets the velocity field at this time
    piv_t = piv_data_sample{t};
    VX_t = piv_t.VX;
    VY_t = piv_t.VY;
    
    %IMPORTANT: MODIFICATION SO THAT GRID IS ORIENTED AT BOTTOM LEFT, NOT TOP LEFT!
    VY_t = -VY_t;
    
%     %computes the velocities of the pointcloud points using the
%     %interpolants of the gridded VX and VY data 
%     pc_t_VX = interp2(X_grid, Y_grid, VX_t, pc_t_X, pc_t_Y, 'nearest');
%     pc_t_VY = interp2(X_grid, Y_grid, VY_t, pc_t_X, pc_t_Y, 'nearest');
%     pc_t_V = [pc_t_VX , pc_t_VY];
    
    %interpolants of VX and VY given their coordinates
    Interp_t_VX = scatteredInterpolant(X_grid(:), Y_grid(:), VX_t(:), 'natural');
    Interp_t_VY = scatteredInterpolant(X_grid(:), Y_grid(:), VY_t(:), 'natural');
    
    %evaluating the interpolants at the pointcloud coordinates
    pc_t_VX = Interp_t_VX(pc_t_X, pc_t_Y);
    pc_t_VY = Interp_t_VY(pc_t_X, pc_t_Y);
    pc_t_V = [pc_t_VX , pc_t_VY];
    
    %advects the pointcloud with these interpolated velocities
    pc_advect = pc_t + pc_t_V;  
        
    %MODIFIED: AFTER advected, makes DV coordinate modulo DV length, because
    %connects to other side!
    pc_advect(:,1) = mod(pc_advect(:,1), scale_DV_resized);
    
    %sorts by the first coordinate so coords in order of least to greatest!
    c1 = pc_advect(:,1);
    c2 = pc_advect(:,2);
    [~,I_sort] = sort(c1, 'ascend');
    
    %this will have the stripe ordered by DV coordinate (1st column)
    pc_advect = [c1(I_sort), c2(I_sort)];
    
    %makes sure none of the coordinates are out of AP bounds
    %does so by selecting all the coordinates which are in bounds
    num_pts = size(pc_advect, 1);
    %temporary structure to store in bounds coordinates
    pc_advect_temp = [];
    %counter to use to add in bounds coordinates
    counter = 1;
    for p = 1 : num_pts
        
        %the AP coordinate of this point
        testCoord = pc_advect(p,2);
        
        %if (only if) the AP coordinate is in bounds, picks it
        if (testCoord >= 0 && testCoord <= scale_AP_resized)
            pc_advect_temp(counter,:) = pc_advect(p,:);
        end
        
        counter = counter+1;
        
    end
    
    %sets the advected point cloud to the in-bounds one
    pc_advect = pc_advect_temp;
        
    %sets the temporary pointcloud to the advected one 
    %(now in bounds, ordered)
    pc_t = pc_advect;


end




%saves the advected stripe 7 structure in the expt directory
stripe7curves = pc_ary_TOSAVE;
save(fullfile(expt_path_sample, stripeSaveHandle), 'stripe7curves');



saveMov = 1;
if (saveMov == 1)

    cd(expt_path_sample)
    v = VideoWriter(strcat('advected_stripe7_mov_program'), 'Uncompressed AVI');
    open(v)
    writeVideo(v, M_sample)
    close(v)

end


disp('Advected Stripe 7 Saved!')

%%

%looking at the runt stripes

clear M_advect

%sample time to look at
t_sample = 25;

%the t0 to use here for the sample
t0V_sample = 20;

%the original t0 for the overall master timeline
t0V_orig = 21;

%path to realspace stripe folder
stripe_dir = '/mnt/data/embryo/Runt Stripe Info/realspacecorr_ss04';

%handle for the stripe itseslf
stripe_handle = 'curve7stats_collapsed_filtered.mat';

%path to stripe
stripe_path = fullfile(stripe_dir, stripe_handle);
%load stripe itself
stripes = load(stripe_path);

%parametrization coords
xx = stripes.xx;

%structure of all leading edge stripes at all times
LX = stripes.LX;

%structure of all trailing edge stripes at all times
TX = stripes.TX;

for t_sample = t0V_orig : 150

%scale for AP and DV used in resized pullbacks
scale_AP_resized = 696;
scale_DV_resized = 820;

%figure(1)

L_stripe_t = LX(:,t_sample);
T_stripe_t = TX(:,t_sample);

%resized coords according to usual dimensions
resized_xx_pix = scale_DV_resized * xx;
%resized_xx_um = resized_pix_to_um * resized_xx_pix;

%resized stripe according to usual dimensions, multiplies by 1-stripe_t
%to flip the orientation properly so A is on top and P is on bottom
resized_L_stripe_t_pix =  scale_AP_resized * L_stripe_t;
%resized_stripe_t_um = resized_pix_to_um * resized_stripe_t_pix;

resized_T_stripe_t_pix = scale_AP_resized * T_stripe_t;

%plot the leading edge stripe data
plot(resized_xx_pix, resized_L_stripe_t_pix, '-b', 'LineWidth', 1.5)
%reverse the active axes
h = gca;
set(h, 'YDir', 'reverse')

xlim([0 scale_DV_resized])
ylim([0 scale_AP_resized])


hold on

% %plot the trailing edge stripe data
% hold on
% plot(resized_xx_pix, resized_T_stripe_t_pix, '-r', 'LineWidth', 1.5);
% 
% xlabel('Position along DV axis (um)')
% ylabel('Position along AP axis (um)')
% 
% %reverse the active axes
% h = gca;
% set(h, 'YDir', 'reverse')
% 
% xlim([0 scale_DV_resized])
% ylim([0 scale_AP_resized])
% 
% title(['Runt Stripe 7 Leading and Trailing Edges at time ', num2str(t_sample), ' min'])
% 
% hold off

%

%figure(2)
% 
% stripe_t = stripe7curves_orig{t_sample};
% 
% scale_fac = 1.6;
% 
% %resized coords according to usual dimensions
% resized_xx_pix = scale_fac * stripe_t(:,2);
% %resized_xx_um = resized_pix_to_um * resized_xx_pix;
% 
% %resized stripe according to usual dimensions, multiplies by 1-stripe_t
% %to flip the orientation properly so A is on top and P is on bottom
% resized_stripe_t_pix = scale_fac * stripe_t(:,1);
% 
% 
% %plot the leading edge stripe data
% plot(resized_xx_pix, resized_stripe_t_pix, '-b', 'LineWidth', 1.5)
%
% xlabel('Position along DV axis (um)')
% ylabel('Position along AP axis (um)')
% 
% %reverse the active axes
% h = gca;
% set(h, 'YDir', 'reverse')
% 
% xlim([0 scale_DV_resized])
% ylim([0 scale_AP_resized])
% 
% title(['Orig and Advected Runt Stripe 7 Edges at time ', num2str(t_sample), ' min'])


%figure(3)

stripe_t = stripe7curves_advected{t_sample-(t0V_sample-1)};

scale_fac = 1.6;

%resized coords according to usual dimensions
resized_xx_pix = scale_fac * stripe_t(:,2);
%resized_xx_um = resized_pix_to_um * resized_xx_pix;

%resized stripe according to usual dimensions, multiplies by 1-stripe_t
%to flip the orientation properly so A is on top and P is on bottom
resized_stripe_t_pix = scale_fac * stripe_t(:,1);

hold on

%plot the leading edge stripe data
plot(resized_xx_pix, resized_stripe_t_pix, '-r', 'LineWidth', 1.5)

xlabel('Position along DV axis (um)')
ylabel('Position along AP axis (um)')

%reverse the active axes
h = gca;
set(h, 'YDir', 'reverse')

xlim([0 scale_DV_resized])
ylim([0 scale_AP_resized])

title(['Master and Advected Runt Stripe 7 Leading Edge at time ', num2str(t_sample), ' min'])

legend('Master', 'Advected')

hold off

pause(.1)
M_advect(t_sample-(t0V_orig-1)) = getframe(gcf);
    

end


saveMov = 1;
if (saveMov == 1)

    cd('/mnt/data/embryo/Marion_Atlas_All_LR-Preserved_Pullbacks/Runt-nanobody/202001150004')
    v = VideoWriter(strcat('orig_and_advected_stripe7_mov_program'), 'Uncompressed AVI');
    open(v)
    writeVideo(v, M_advect)
    close(v)

end
