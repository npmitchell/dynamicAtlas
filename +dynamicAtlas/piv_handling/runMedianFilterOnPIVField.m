function runMedianFilterOnPIVField(input_PIV_path, output_path, median_order)

% runMedianFilterOnPIVField(input_PIV_path, output_path)
%   Runs a median filter on the PIV data located in input_PIV_path and 
%   saves the data to output_path
%   
% Parameters
% ----------
% input_PIV_dir : path to directory containing all the PIV .mat structures
% output_dir : path to directory in which to save median filtered PIV
% median_order : order of median filter--total range to take median over
%
% Outputs
%
% new directory with name output_directory that containts saved PIV fields 
% which have been run through median filter
% 
%
% Vishank Jain-Sharma 2020

%navigate to path with all the original PIV mat structures
cd(input_PIV_path)

%finds all .mat files within the directory
d = dir('*.mat');

%length of d will be the number of PIV timepoints
numT = length(d(:));

%finds spatial dimensions of the PIV matrices, all will be identical
test_piv = load(fullfile(input_PIV_path, d(1).name));
dim1 = size(test_piv.VX, 1);
dim2 = size(test_piv.VX, 2);

%3d arrays to store the PIV structures, 3rd dim is number of timepoints
piv_master_x = zeros(dim1, dim2, numT);
piv_master_y = zeros(dim1, dim2, numT);

%loads in all the PIV structures for X and Y components separately
for t = 1 : numT
    piv_temp = load(fullfile(input_PIV_path, d(t).name));
    piv_master_x(:,:,t) = piv_temp.VX;
    piv_master_y(:,:,t) = piv_temp.VY;
end

%3d arrays to store the filtered PIV structures
piv_filtered_x = zeros(dim1, dim2, numT);
piv_filtered_y = zeros(dim1, dim2, numT);

%performs median filter on each vector in time in the array, will be 
%dim1*dim2 medfilt1 operations performed. Uses median_order for order of
%medfilt1, representing the number of entries which median is taken over
for i = 1 : dim1
    for j = 1 : dim2
        piv_filtered_x(i,j,:) = medfilt1(piv_master_x(i,j,:), median_order);
        piv_filtered_y(i,j,:) = medfilt1(piv_master_y(i,j,:), median_order);
    end
end

%makes output directory
mkdir(output_path);

%extracts components from each timepoint, puts into struct, and saves
%to output_path folder
for t = 1 : numT
    
    %x and y components
    VX = piv_filtered_x(:,:,t);
    VY = piv_filtered_y(:,:,t);
    
    %file name used to save the filtered struct
    save_filename = fullfile(output_path, sprintf('/VeloT_medfilt_%06d.mat', t));
    
    %saves the filtered piv for this timepoint to the given file name
    save(save_filename, 'VX', 'VY');
    
    clear VX VY
end


end