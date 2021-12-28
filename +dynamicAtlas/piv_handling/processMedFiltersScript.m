%script to run med filters on PIV data

clear all
close all
clc

%%

%directory with desired pullback data
genotypeDirectory = '/Users/Vishank/Documents/Pullback_data/WT/CAAX-mCherry';

cd(genotypeDirectory)

%extension BEFORE adding filtered
piv_extension = 'PTV';
%piv_extension = 'PIV_equalized';

%all dataset folders within dir
d = dir('20*');

%order of median filters
median_order = 3;

%loops over all folders and runs median filters on all of them
for i = 1 : length(d(:))
    
    dataset_name = d(i).name;
    input_PIV_path = fullfile(genotypeDirectory, dataset_name, piv_extension);
    output_path = fullfile(genotypeDirectory, dataset_name, [piv_extension, '_filtered']);
    
    runMedianFilterOnPIVField(input_PIV_path, output_path, median_order);
    
end

%navigates back into genotype directory
cd(genotypeDirectory)

disp('Done with processMedFiltersScript.m')