% script_example_relaxPairwiseCorrespondenceNetwork

addpath('/mnt/data/code/dynamicAtlas/+dynamicAtlas/data_handling/')
addpath('/mnt/data/code/dynamicAtlas/timeline_handling/')

%% A network to stretch
ttc = {} ; ttc{1} = {} ; ttc{2} = {}; ttc{3} = {} ;
ttc{1}{1} = cat(2, linspace())';
ttc{1}{2} = [1,2.5,3;1,2,3]'; 
ttc{1}{3} = [1,2.5,3;1,2,3]'; 
ttc{2}{2} = [1,2,3;1,2,3]';
ttc{2}{3} = [1,2,3;1,2,3]';
ttc{3}{3} = [1,2,3;1,2,3]';
options = struct('save', false, 'use_offdiagonals_only', true) ;
featureMatchedStr = 'test' ;
timelineDir = '' ;
expts = {'testa','testb','testc'} ;
exptIDs = {'testa','testb', 'testc'} ;
hard = 1;
relaxPairwiseCorrespondenceNetwork(ttc, hard, ...
expts, exptIDs, timelineDir, featureMatchedStr, options) 


%% A network with no energy and no net force
ttc = {} ; ttc{1} = {} ; ttc{2} = {}; ttc{3} = {} ;
ttc{1}{1} = [1,2,3;1,2,3]';
ttc{1}{2} = [1,2.5,3;1,2,3]'; 
ttc{1}{3} = [1,2.5,3;1,2,3]'; 
ttc{2}{2} = [1,2,3;1,2,3]';
ttc{2}{3} = [1,2,3;1,2,3]';
ttc{3}{3} = [1,2,3;1,2,3]';
options = struct('save', false, 'use_offdiagonals_only', true) ;
featureMatchedStr = 'test' ;
timelineDir = '' ;
expts = {'testa','testb','testc'} ;
exptIDs = {'testa','testb', 'testc'} ;
hard = 1;
relaxPairwiseCorrespondenceNetwork(ttc, hard, ...
expts, exptIDs, timelineDir, featureMatchedStr, options) 


%% A network with energy and no net force except for dt within expts
ttc = {} ; ttc{1} = {} ; ttc{2} = {}; ttc{3} = {} ;
ttc{1}{1} = [1,2,3; 1,2,3]';
ttc{1}{2} = [1,2.5,3; 1,2,3]';  % timeline d2 maps 2 to 2.5 in d1
ttc{1}{3} = [1,2.5,3; 1,2,3]';  % timeline d3 maps 2 to 2.5 in d1
ttc{2}{1} = [1,1.66,3; 1,2,3]';  % timeline d1 maps 2 to 1.5 in d2
ttc{2}{2} = [1,2,3; 1,2,3]';
ttc{2}{3} = [1,2,3; 1,2,3]';
ttc{3}{1} = [1,1.66,3; 1,2,3]';  % timeline d1 maps 2 to 1.5 in d3
ttc{3}{2} = [1,2,3; 1,2,3]'; 
ttc{3}{3} = [1,2,3; 1,2,3]';
options = struct('save', false, 'use_offdiagonals_only', true) ;
featureMatchedStr = 'test' ;
timelineDir = '' ;
expts = {'testa','testb','testc'} ;
exptIDs = {'testa','testb', 'testc'} ;
hard = 1;
relaxPairwiseCorrespondenceNetwork(ttc, hard, ...
expts, exptIDs, timelineDir, featureMatchedStr, options) 


%% A network with energy and net force
ttc = {} ; ttc{1} = {} ; ttc{2} = {}; ttc{3} = {} ;
ttc{1}{1} = [1,2,3; 1,2,3]';
ttc{1}{2} = [1,1.5,2.5; 1,2,3]';  % timeline d2 maps 2 to 2.5 in d1
ttc{1}{3} = [1,1.5,2.5; 1,2,3]';  % timeline d3 maps 2 to 2.5 in d1
% ttc{2}{1} = [1,2,3; 1,2,3]';  % timeline d1 maps 2 to 1.5 in d2
ttc{2}{2} = [1,2,3; 1,2,3]';
ttc{2}{3} = [1,2,3; 1,2,3]';
% ttc{3}{1} = [1,2,3; 1,2,3]';  % timeline d1 maps 2 to 1.5 in d3
% ttc{3}{2} = [1,2,3; 1,2,3]'; 
ttc{3}{3} = [1,2,3; 1,2,3]';
options = struct('save', false, 'use_offdiagonals_only', true) ;
featureMatchedStr = 'test' ;
timelineDir = '' ;
expts = {'testa','testb','testc'} ;
exptIDs = {'testa','testb', 'testc'} ;
hard = 1;
i_tau0j_tau0jrelaxed = relaxPairwiseCorrespondenceNetwork(ttc, hard, ...
    expts, exptIDs, timelineDir, featureMatchedStr, options) 

%% A network with energy and net force
ttc = {} ; ttc{1} = {} ; ttc{2} = {}; ttc{3} = {} ;
ttc{1}{1} = [1,2,3; 1,2,3]';
ttc{1}{2} = [1,1.5,2.5; 1,2,3]';  % timeline d2 maps 2 to 2.5 in d1
ttc{1}{3} = [1,1.5,2.5; 1,2,3]';  % timeline d3 maps 2 to 2.5 in d1
ttc{2}{1} = [1,2,3; 1,2,3]';  % timeline d1 maps 2 to 1.5 in d2
ttc{2}{2} = [1,2,3; 1,2,3]';
ttc{2}{3} = [1,2,3; 1,2,3]';
ttc{3}{1} = [1,2,3; 1,2,3]';  % timeline d1 maps 2 to 1.5 in d3
ttc{3}{2} = [1,2,3; 1,2,3]'; 
ttc{3}{3} = [1,2,3; 1,2,3]';
options = struct('save', false, ...
    'use_offdiagonals_only', true, ...
    'timelineStiffness', 0.25) ;
featureMatchedStr = 'test' ;
timelineDir = '' ;
expts = {'testa','testb','testc'} ;
exptIDs = {'testa','testb', 'testc'} ;
hard = 1;
i_tau0j_tau0jrelaxed = relaxPairwiseCorrespondenceNetwork(ttc, hard, ...
    expts, exptIDs, timelineDir, featureMatchedStr, options) 



