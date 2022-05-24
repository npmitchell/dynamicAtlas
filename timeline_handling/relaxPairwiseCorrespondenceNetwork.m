function [i_tau0j_tau0jrelaxed, itau0jInfo, NLKLBLInfo] =...
    relaxPairwiseCorrespondenceNetwork(ttc, hard, ...
    expts, exptIDs, timelineDir, featureMatchedStr, options) 
% relaxPairwiseCorrespondenceNetwork(ttc, hard, expts, exptIDs, timelineDir, featureMatchedStr, overwrite, roundTimepoints)
%
% Parameters
% ----------
% ttc: cell of cell arrays
%   ttc{ii}{jj} is an Nx2 array where 
%   ttc{ii}{jj}(:, 1) is timeline indices for ii, as integers >=1, and 
%   ttc{ii}{jj}(:, 2) is the corresponding placement into timeline jj, 
%   as floats >=1.
% hard : int
%   the index of the dataset for which to make 
% expts : cell of string paths
%   paths to the experiments used in the timelines, only used if
%   options.saveResults is true
% exptIDs : cell of string identifiers
%   name identifiers of each experiment, such as '202001011200', only used
%   if options.saveResults is true
% timelineDir : str path
%   where the timeline will be stored on disk, if options.saveResults is true
% featureMatchedStr : string (ex 'ssR', '_curve7_chisq', 'pbPathlines')
%    description of what feature is matched to build this network 
% options : struct with fields
%   overwrite : bool
%       overwrite existing results on disk
%   roundTimepoints : bool
%       round the relaxed timepoints to nearest integer values (for use as
%       indices)
%   use_offdiagonals_only : bool
%       compare only 1->1, 1->2, 1->3, ... 
%                          2->2, 2->3, ...
%                                3->3, ... etc
%       instead of   1->1, 1->2, 1->3, ... 
%                    2->1, 2->2, 2->3, ...
%                    3->1, 3->2, 3->3, ... etc
%   saveResults : bool (default=true)
%       write results to disk
%   ensureReciprocalInteractions : bool (default=true)
%       if there is a bond linking ii to jj, also send bond with same
%       strength from jj to ii. Useful if use_offdiagonals_only==true
%   dts : #expts x 1 cell of numeric or of #timepoints x 1 numeric arrays
%       timestep between integer-indexed frames for each experiment
%
% Returns
% -------
%
% Saves to disk
% -------------
% itau0jfn = fullfile(timelineDir, 'i_tau0j.mat') ;
%   saved variables: 'i_tau0j', 'ntp_tot', 'addt', 'ntps', ...
%                    'exptIDs', 'hard', 'hard_reference_ID'
%
%
% Example usage
% -------------
% ttc = {} ; ttc{1} = {} ;
% ttc{1}{1} = [1,2,3;1,2,3]';
% ttc{2}{2} = [1,2,3;1,2,3]';
% ttc{3}{3} = [1,2,3;1,2,3]';
% ttc{1}{2} = [1,2,3;1,3,3]'; ttc{1}{3} = [1,2,3;1,3,3]'; 
% test{2}{3} = [1,2,3;1,2,3]';
% options = struct('saveResults', false, 'use_offdiagonals_only', true) ;
% featureMatchedStr = 'test' ;
% timelineDir = '' ;
% expts = {'testa','testb','testc'} ;
% exptIDs = {'testa','testb', 'testc'} ;
% hard = 1;
% relaxPairwiseCorrespondenceNetwork(ttc, hard, ...
%    expts, exptIDs, timelineDir, featureMatchedStr, options) 
%
% See also 
% --------
% script_example_relaxPairwiseCorrespondences.m
%
% NPMitchell 2020-2022


%% Default options
maxDistReference = 30 ;
roundTimepoints = false ;
overwrite = false ;
ensureReciprocalInteractions = false ;
use_offdiagonals_only = false ;

% Parse options
if nargin < 7
    options = struct() ;
end

if isfield(options, 'overwrite')
    overwrite = options.overwrite ;
end
if isfield(options, 'saveResults')
    saveResults = options.saveResults ;
else
    saveResults = ~isempty('timelineDir') ;
end
if isfield(options, 'roundTimepoints')
    roundTimepoints = options.roundTimepoints;
end

if isfield(options, 'ensureReciprocalInteractions')
    ensureReciprocalInteractions = options.ensureReciprocalInteractions ;
end
if isfield(options, 'use_offdiagonals_only')
    use_offdiagonals_only = options.use_offdiagonals_only ;
end
if isfield(options, 'maxDistReference')
    maxDistReference = options.maxDistReference ;
end

% load any non-standard timesteps
if isfield(options, 'dts')
    dts = options.dts ;
else
    % default timestep is just 1 for each adjacent timepoint pair
    dts = cell(length(ttc), 1) ;
    for qq = 1:length(ttc)
        dts{qq} = ones(size(ttc{qq}{qq},1) - 1, 1) ;
    end
end

% Determine relative stiffness of self-timeline and cross-timeline
if isfield(options, 'timelineStiffness')
    timelineStiffness = options.timelineStiffness ;
else
    timelineStiffness = 1.0 ;
end

if isempty(featureMatchedStr)
    featureMatchedStr = '_curve7_chisq';
end
hard_reference_ID = exptIDs{hard} ;

%% VIII. Relax timepoints to reference curve (time of dataset #hard)
% This is assigning an intial guess to the the entire timeline for tau0.
% This will be updated in the next sections so its OK that its just an
% intial guess and only includes for example TTC{1}{hard} and not TTC
% {hard}{1}.
% so-called 'hard reference' is the master timeline

% Fit other curves to the hard reference timeline
% Tau0 is a map from timeline ti to reference timeline
% to use, do yy = ppval(pp, xq);
for ii = 1:length(ttc)
    % decide if we use all-to-all (1->2 and 2->1) or half-all-to-all (1->2 only)
    if use_offdiagonals_only
        if ii <= hard
            pp = spline(ttc{ii}{hard}(:, 1), ttc{ii}{hard}(:, 2)) ;
        elseif ii > hard
            pp = spline(ttc{hard}{ii}(:, 2), ttc{hard}{ii}(:, 1)) ;
        end
    else
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
    
    itau0jInfo = struct('i_tau0j', i_tau0j, ...
        'ntp_tot', ntp_tot, 'addt', addt, 'ntps', ntps, ...
        'exptIDs', exptIDs, ...
        'hard', hard, ...
        'hard_reference_ID', hard_reference_ID) ;
else
    % First get #tps in each dataset, which is stored in diagonals of ttc
    ntps = zeros(length(ttc), 1) ;
    addt = zeros(length(ttc), 1) ;
    tpid = [] ;
    for ii = 1:length(ttc)
        addt(ii) = sum(ntps) ;
        if ~isempty(ttc{ii}{ii})    
            ntps(ii) = length(ttc{ii}{ii}) ;  
        end
        % Now label each timepoint in linear indexing
        tpid = cat(2, tpid, addt(ii) + (1:ntps(ii))) ;
        
        % Store all indices in one array, evaluating polynomial interpolant
        % of tau0 here
        tau0extra = fnxtr(tau0(ii), 2) ;
        tauv = ppval(tau0extra, (1:ntps(ii))) ;
        i_tau0j(addt(ii) + (1:ntps(ii)), :) = [ii * ones(ntps(ii), 1), tauv'] ;
    end
    ntp_tot = sum(ntps) ;
    
    % save data
    if nargout > 0
        itau0jInfo = struct('i_tau0j', i_tau0j, ...
            'ntp_tot', ntp_tot, 'addt', addt, 'ntps', ntps, ...
            'exptIDs', exptIDs, ...
            'hard', hard, ...
            'hard_reference_ID', hard_reference_ID) ;
    end
    if saveResults
        save(itau0jfn, 'i_tau0j', 'ntp_tot', 'addt', 'ntps', 'exptIDs', ...
            'hard', 'hard_reference_ID') ;
    end
    
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
    tmp = load(NLKLBLfn, 'NL', 'KL', 'BL', 'BL0', 'kL', ...
        'exptIDs','hard_reference_ID') ;
    NL = tmp.NL ;
    KL = tmp.KL ;
    BL = tmp.BL ;
    BL0 = tmp.BL0 ;
    kL = tmp.kL ;
    
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
        error(['Hard reference stored with NLKLBL does not match current. ', ...
            'Remove ' NLKLBLfn '  from disk and rerun'])
    end
    
    NLKLBLInfo = struct('NL', NL, ...
        'KL', KL, 'BL', BL, 'BL0', BL0, 'kL', kL, 'exptIDs', exptIDs, ...
        'hard', hard, 'hard_reference_ID', hard_reference_ID) ;
else
    % preallocate Nxlarge array for NL, KL, do not preallocate BL
    NL = zeros(ntp_tot, 300) ;
    KL = zeros(ntp_tot, 300) ;
    BM0 = zeros(ntp_tot, 300) ;
    first = true ;
    for ii = 1:length(ttc)
        for jj = 1:length(ttc)
            if ii ~= jj
                % Ensure that there are some correspondences
                if ~isempty(ttc{ii}{jj})
                    if ~roundTimepoints
                        msg = ['Two bonds ascribed to each node ', ...
                            'since timestamping is non-integer, ', ...
                            'each bond stiffness weighted by distance from partner'] ;
                        disp(msg)
                        disp(['Asserting that partners are integer timestamps of timeline ' num2str(ii)])
                        assert(all(mod(ttc{ii}{jj}(:, 1), 1) == 0))
                        
                        % if the match is a perfect integer, add a single
                        % bond
                        for qq = 1:length(ttc{ii}{jj}(:, 1))
                            if ceil(ttc{ii}{jj}(qq, 2)) == floor(ttc{ii}{jj}(qq, 2)) && ttc{ii}{jj}(qq, 2) >= 1
                                % the match is a perfect integer
                                nodei = ttc{ii}{jj}(qq, 1) + addt(ii);
                                nodej = ttc{ii}{jj}(qq, 2) + addt(jj);
                                bladd = [nodei, nodej] ;
                                kladd = [1]  ;
                                bl0add = 0 ;
                            elseif ttc{ii}{jj}(qq, :) >= 1
                                % Add two bonds since the match is non-integer
                                % ceil
                                if ceil(ttc{ii}{jj}(qq, 2)) <= ntps(jj) && floor(ttc{ii}{jj}(qq, 2)) >= 1
                                    if ceil(ttc{ii}{jj}(qq, 2)) <= ntps(jj)
                                        nodei = ttc{ii}{jj}(qq, 1) + addt(ii);
                                        nodej = ceil(ttc{ii}{jj}(qq, 2)) + addt(jj);
                                        bladd1 = [nodei, nodej] ;
                                        % Make bond strongest when timepoints
                                        % very nearly match an integer given by
                                        % ceil(ttc{ii}{jj}(qq,1)), and
                                        % weakest when nearly match floor().
                                        kladd1 = abs(1 - (ceil(ttc{ii}{jj}(qq, 2)) - ttc{ii}{jj}(qq, 2)) ) ;
                                        % bond length is position minus node
                                        % position (so it is negative)
                                        bl0add1 = ttc{ii}{jj}(qq, 2) - ceil(ttc{ii}{jj}(qq, 2))  ;
                                    else
                                        error('handle edge point here -- sample extends beyond reference')
                                    end

                                    % floor
                                    if floor(ttc{ii}{jj}(qq, 2)) >= 1
                                        nodei = ttc{ii}{jj}(qq, 1) + addt(ii);
                                        nodej = floor(ttc{ii}{jj}(qq, 2)) + addt(jj);
                                        bladd2 = [nodei, nodej] ;
                                        % Make bond strongest when timepoints
                                        % very nearly match an integer given by
                                        % floor(ttc{ii}{jj}(qq,1)), and
                                        % weakest when nearly match ceil().
                                        kladd2 = abs(1 - (ttc{ii}{jj}(qq, 2) - floor(ttc{ii}{jj}(qq, 2))) ) ;
                                        % bond length is position minus node
                                        % position (so it is positive)
                                        bl0add2 = ttc{ii}{jj}(qq, 2) - floor(ttc{ii}{jj}(qq, 2))  ;
                                    else
                                        error('handle edge point here -- sample extends beyond reference')
                                    end

                                    % cat them
                                    bladd = [bladd1; bladd2] ;
                                    kladd = [kladd1; kladd2] ;
                                    bl0add = [bl0add1; bl0add2];
                                    
                                    try
                                        assert(sum(kladd) == 1)
                                        assert(sum(abs(bl0add)) == 1)
                                    catch
                                        error('bond strengths did not sum to 1')
                                    end
                                else
                                    disp('index out of range, ignoring')
                                    bladd = [] ;
                                    kladd = [] ;
                                    bl0add = [] ;
                                end
                            else
                                disp('index out of range, ignoring')
                                bladd = [] ;
                                kladd = [] ;
                                bl0add = [] ;
                            end

                            % Check that these are not too far apart in the
                            % hard reference timeline
                            if ~isempty(bladd) 
                                try
                                    distInRefTimeline = abs(i_tau0j(bladd(:, 1), 2) - i_tau0j(bladd(:, 2), 2)) ;
                                    assert(all(distInRefTimeline < maxDistReference))
                                catch
                                    error('This bond would be too long in the reference dataset timeline')
                                end
                            end
                            
                            % Add to BL
                            if first
                                BL = bladd ;
                                first = false ;
                                BL0 = bl0add ;
                                kL = kladd ;
                            else
                                BL = cat(1, BL, bladd)  ;
                                BL0 = cat(1, BL0, bl0add) ; 
                                kL = cat(1, kL, kladd)  ;
                                assert(size(BL0, 2) == 1)
                            end
                            % build NL, KL with redundancy and non-reciprocities
                            %   of correspondences built into KL
                            for id = 1:size(bladd,1)
                                pair = bladd(id, :) ;
                                kv = kladd(id) ;
                                bl0 = bl0add(id) ;
                                nodei = pair(1) ; 
                                nodej = pair(2) ;

                                % ij. Note that we treat each pairing anew 
                                % to allow for non-reciprocal bonds
                                % if ismember(nodej, NL(nodei, :))
                                %     ind = find(NL(nodei, :) == nodej) ;
                                %     assert(NL(nodei, ind) == nodej) ;
                                %     KL(nodei, ind) = KL(nodei, ind) + kv ;
                                % else
                                firstzero = find(NL(nodei, :)==0, 1) ; 
                                assert(~isempty(firstzero)) 
                                NL(nodei, firstzero) = nodej ; 
                                KL(nodei, firstzero) = kv ;
                                BM0(nodei, firstzero) = bl0 ;

                                % % ji
                                % if ismember(nodei, NL(nodej, :))
                                %     ind = find(NL(nodej, :) == nodei) ;
                                %     assert(NL(nodej, ind) == nodei) ;
                                %     KL(nodej, ind) = KL(nodej, ind) + kv ;
                                % else
                                %     firstzero = find(NL(nodej, :)==0, 1) ; 
                                %     if isempty(firstzero)
                                %        tester = 'testing here' ;
                                %        disp(tester)
                                %     end
                                %     assert(~isempty(firstzero)) 
                                %     NL(nodej, firstzero) = nodei ; 
                                %     KL(nodej, firstzero) = kv ;
                                % end
                            end
                        end
                    else
                        error('check here')
                        nodei = round(ttc{ii}{jj}(:, 1)) + addt(ii);
                        nodej = round(ttc{ii}{jj}(:, 2)) + addt(jj);
                        bladd = [nodei, nodej] ;
                        bl0add = 0 ;

                        % Add to BL
                        if first
                            BL = bladd ;
                            BL0 = bl0add ;
                            first = false ;
                            kL = kladd ;
                        else
                            BL = cat(1, BL, bladd) ;
                            BL0 = cat(1, BL0, bl0add) ;
                            kL = cat(1, kL, kl0add) ;
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

                            % ji
                            if ismember(nodei, NL(nodej, :))
                                ind = find(NL(nodej, :) == nodei) ;
                                assert(NL(nodej, ind) == nodei) ;
                                KL(nodej, ind) = KL(nodej, ind) + 1 ;
                            else
                                firstzero = find(NL(nodej, :)==0, 1) ; 
                                if isempty(firstzero)
                                   tester = 'testing here' ;
                                   disp(tester)
                                end
                                assert(~isempty(firstzero)) 
                                NL(nodej, firstzero) = nodei ; 
                                KL(nodej, firstzero) = 1 ;
                            end
                        end
                    end
                end
            else
                % Add bonds between nodes for self-timeline (bond length dt)
                for qq = 1:length(ttc{ii}{jj}(:, 1))-1
                    nodei = qq + addt(ii);
                    nodej = qq + 1 + addt(jj);
                    bladd = [nodei, nodej] ;
                    kladd = timelineStiffness  ;
                    bl0add = dts{ii}(qq) ;

                    % Add to BL
                    if first
                        BL = bladd ;
                        first = false ;
                        BL0 = bl0add ;
                        kL = kladd ;
                    else
                        BL = cat(1, BL, bladd) ; 
                        BL0 = cat(1, BL0, bl0add) ;
                        kL = cat(1, kL, kladd) ;
                        assert(size(BL0, 2) == 1)
                    end
                    % build NL, KL with any redundancy and non-reciprocities
                    %   of correspondences built into KL
                    for id = 1:size(bladd,1)
                        pair = bladd(id, :) ;
                        kv = kladd(id) ;
                        bl0 = bl0add(id) ;
                        nodei = pair(1) ; 
                        nodej = pair(2) ;

                        if(any([nodei, nodej] == 113))
                            disp('pause')
                        end

                        % ij. Note that we treat each pairing anew 
                        % to allow for non-reciprocal bonds
                        firstzero = find(NL(nodei, :)==0, 1) ; 
                        assert(~isempty(firstzero)) 
                        NL(nodei, firstzero) = nodej ; 
                        KL(nodei, firstzero) = kv ;
                        BM0(nodei, firstzero) = bl0 ;
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
    NLKLBLInfo = struct('NL', NL, ...
        'KL', KL, 'BL', BL, 'BL0', BL0, 'kL', kL, 'exptIDs', exptIDs, ...
        'hard', hard, 'hard_reference_ID', hard_reference_ID) ;
    if saveResults
        save(NLKLBLfn, 'NL', 'KL', 'BL', 'BL0', 'kL', ...
            'exptIDs', 'hard', 'hard_reference_ID')
    end
end

%% X. Build timeline network visualization (pairwise_corr_timeline_XXX.png)
colorset = define_colors() ;
close all
for use_BL = [true false]
    if use_BL
        disp('Visualizing network using BL')
    else
        disp('Visualizing network using KL, NL')
    end
    for ii = 1:1:length(ttc)  
        if saveResults
            close all
            subplot(1,2,1)
        else
            subplot(1,length(ttc),ii)
        end
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
                kv = kL(qq) ;
                % col = find(NL(BL(qq, 1), :) == BL(qq, 2)) ;
                % kv = KL(BL(qq, 1), col) ;
                % if we only use off-diagonals, then everything should be
                % symmetric
                if ensureReciprocalInteractions
                    col2 = find(NL(BL(qq, 2), :) == BL(qq, 1)) ;
                    try
                       assert(KL(BL(qq, 2), col2) == kv)
                    catch
                       error('interactions are non-reciprocal')
                    end
                end
                if kv > 0
                    if i_tau0j(BL(qq, 1), 1) == ii 
                        if roundTimepoints
                            plot(round([i_tau0j(BL(qq, 1), 2), i_tau0j(BL(qq, 2), 2)]), ...
                               round([i_tau0j(BL(qq, 1), 1), i_tau0j(BL(qq, 2), 1)]), '-', ...
                               'color', colorset(i_tau0j(BL(qq, 2), 1), :), 'linewidth', 2*kv)
                        else
                            % non-rounded
                            plot([i_tau0j(BL(qq, 1), 2), i_tau0j(BL(qq, 2), 2)], ...
                                 [i_tau0j(BL(qq, 1), 1), i_tau0j(BL(qq, 2), 1)], '-', ...
                                 'color', colorset(i_tau0j(BL(qq, 2), 1), :), 'linewidth', 2*kv)
                        end
                    end
                    if ensureReciprocalInteractions
                        if i_tau0j(BL(qq, 2), 1) == ii
                            if roundTimepoints
                                plot(round([i_tau0j(BL(qq, 1), 2), i_tau0j(BL(qq, 2), 2)]), ...
                                     round([i_tau0j(BL(qq, 1), 1), i_tau0j(BL(qq, 2), 1)]), '-', ...
                                     'color', colorset(i_tau0j(BL(qq, 1), 1), :), 'linewidth', 2*kv)
                            else
                                % non-rounded                     
                                plot([i_tau0j(BL(qq, 1), 2), i_tau0j(BL(qq, 2), 2)], ...
                                     [i_tau0j(BL(qq, 1), 1), i_tau0j(BL(qq, 2), 1)], '-', ...
                                     'color', colorset(i_tau0j(BL(qq, 1), 1), :), 'linewidth', 2*kv)
                            end
                        end
                    
                    end
                
                end
            end

            if isnumeric(exptIDs{ii})
                title(['Time relaxation network: dataset ' num2str( exptIDs{ii} ) ])
            else
                title(['Time relaxation network: dataset ' exptIDs{ii}])
            end
            ylim([0, max(i_tau0j(:, 1)) + 1])
            if saveResults
                subplot(1,2,2)
                for pairID = 1:length(ttc{ii})
                    pair = ttc{ii}{pairID} ;
                    plot(pair(:, 1), pair(:, 2), '.-', 'color', colorset(pairID, :)) 
                    hold on;
                end
                axis equal
                axis square
                xlabel(['timeline i, ' exptIDs{ii}])
                ylabel('timeline j')
                saveas(gcf, fullfile(timelineDir, sprintf('pairwise_corr_timeline_BL_%03d.pdf', ii)))
            end
            % disp('press any button on figure')
            pause(0.1)
        else
            % Add bonds to plot using NL,KL 
            % cycle only through the linear point indices that belong to ii
            for qq = (addt(ii)+1):(addt(ii) + ntps(ii))
                for dmyk = 1:length(NL(qq,:))
                    nn = NL(qq, dmyk) ;
                    kv = KL(qq, dmyk) ;
                    % plot only if real neighbor
                    if nn > 0
                        plot([i_tau0j(qq, 2), i_tau0j(nn, 2)], ...
                             [i_tau0j(qq, 1), i_tau0j(nn, 1)], '-', ...
                             'color', colorset(i_tau0j(nn, 1), :), 'linewidth', 2*kv)
                    end
                end
            end
            title(['Time relaxation network: dataset ' exptIDs{ii}])
            ylim([0, max(i_tau0j(:, 1)) + 1])
            if saveResults
                subplot(1,2,2)
                for pairID = 1:length(ttc{ii})
                    pair = ttc{ii}{pairID} ;
                    plot(pair(:, 1), pair(:, 2), '.-', 'color', colorset(pairID, :)) 
                    hold on;
                end
                axis equal
                xlabel(['timeline i, ' exptIDs{ii}])
                ylabel('timeline j')
                saveas(gcf, fullfile(timelineDir, sprintf('pairwise_corr_timeline_%03d.pdf', ii)))
            end
            % disp('press any button')
            pause(0.1)
        end
    end
    if ~saveResults
        msg = 'Inspect figures and press Continue button on figure to continue' ;
        continueButtonOnFigure(gcf, msg) ;
        clf
    end
end
disp('done with network visualization')
close all

%% XII. Relax the timeline network & visualize (i_tau0j_tau0jrelaxed.mat)

% check for existing result
close all
fn = fullfile(timelineDir, 'i_tau0j_tau0jrelaxed.mat') ;
if exist(fn, 'file') && ~overwrite
    load(fn, 'i_tau0j_tau0jrelaxed')
else
    % Relax: fix the time nodes of the hard reference dataset, move all others
    options = optimset('PlotFcns',@optimplotfval, 'TolX', 1e-4, ...
        'TolFun', 1e-6, 'MaxIter', 1e4);
    x0 = i_tau0j(:, 2) ;
    % pop the indices of the fixed times from the array x0.
    % Here, fix only one single timepoint to avoid zero mode translation!
    fixed_ind = find(i_tau0j(:, 1) == hard, 1) ;
    fixed_xx = i_tau0j(fixed_ind, 2) ;
    x0(fixed_ind) = [] ;
    fun = @(x)springEnergy1D(x, BL, fixed_ind, fixed_xx, kL, BL0);
    xf = fminsearch(fun, x0, options) ;
    % reinsert indices of fixed times
    xnew1 = xf(1:fixed_ind(1)-1) ;
    xnew2 = xf(fixed_ind(1):end) ;
    xnew = [xnew1; fixed_xx; xnew2 ] ;
    disp('done')

    % Plot the relaxed network result 
    if saveResults
        close all
    else
        figure()
    end
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
    h1 = plot([NaN], [NaN], 'ko') ;
    h2 = plot([NaN], [NaN], 'k.') ;
    legend([h1, h2], {'timestamp in reference timeline', 'annealed timestamp'})
    
    % Figure labels
    ylabel('dataset $i$', 'Interpreter', 'Latex')
    xlabel('time, $t_0$', 'Interpreter', 'latex')
    xlim([0, max(i_tau0j(:, 2)) + 1])
    
    % Title and formatting
    title('Time relaxation network')
    ylim([0, max(i_tau0j(:, 1)) + 1])
    set(gcf, 'Units', 'centimeters');
    set(gcf, 'Position', [0 0 36 16]) ;
    if saveResults
        saveas(gcf, fullfile(timelineDir, sprintf('relaxation_results.pdf')))
    end
    % Save the result
    i_tau0j_tau0jrelaxed = cat(2, i_tau0j, xnew) ;
    if saveResults
        save(fn, 'i_tau0j_tau0jrelaxed')
    else
        % msg = 'Time relaxation network: press Enter to continue' ;
        % continueButtonOnFigure(gcf, msg) ;
    end
end
disp('done')

%% Show relaxed timelines
% % Plot the dynamics
% close all
% colorset = define_colors() ;
% for use_offset = [true false]
%     offset = 0 ;
%     for ii = 1:length(ttc)
%         subplot(round(length(ttc)*0.5)+1, 2, ii)
%         
%         if useOffDiagonalsOnly
%             % Do earlier dataset correspondences
%             for jj = 1:ii-1
%                 if ~isempty(ttc{jj}{ii})
%                     % displace vertically to match
%                     if use_offset
%                         i0 = ttc{jj}{ii}(1, 2) ;
%                         j0 = ttc{jj}{ii}(1, 1) ;
%                         offset = i0 - j0 ;
%                     end
%                     plot(ttc{jj}{ii}(:, 2), ...
%                         ttc{jj}{ii}(:, 1) + offset,...
%                         '.-', 'color', colorset(jj, :))
%                     hold on;
%                 end
%             end
% 
%             % Do this dataset to itself
%             for jj = ii 
%                 if ~isempty(ttc{ii}{jj})
%                     plot(ttc{ii}{jj}(:, 1), ttc{ii}{jj}(:, 2), '.', ...
%                         'color', colorset(jj, :))
%                     hold on;
%                 end
%             end
% 
%             % Do later dataset correspondences
%             for jj = ii+1:length(ttc)
%                 if ~isempty(ttc{ii}{jj})
%                     % displace vertically to match
%                     if use_offset
%                         i0 = ttc{ii}{jj}(1, 1) ;
%                         j0 = ttc{ii}{jj}(1, 2) ;
%                         offset = i0 - j0 ;
%                     end
%                     plot(ttc{ii}{jj}(:, 1), ...
%                         ttc{ii}{jj}(:, 2) + offset,...
%                         '.-', 'color', colorset(jj, :))
%                     hold on;
%                 end
%             end
%         else
%             
%             % Do all correspondences
%             for jj = 1:length(ttc)
%                 if ~isempty(ttc{ii}{jj})
%                     % displace vertically to match
%                     if use_offset
%                         i0 = ttc{ii}{jj}(:, 1) ;
%                         j0 = ttc{ii}{jj}(:, 2) ;
%                         offset = mean(i0 - j0) ;
%                     end
%                     plot(ttc{ii}{jj}(:, 1), ...
%                         ttc{ii}{jj}(:, 2) + offset,...
%                         '.-', 'color', colorset(jj, :))
%                     hold on;
%                 end
%             end            
%         end
% 
%         % labels
%         if isnumeric(exptIDs{ii})
%             ylabel(['$\tau(t_{' num2str(exptIDs{ii}) '})$' ], ...
%                 'Interpreter', 'Latex')
%         else
%             ylabel(['$\tau(t_{' exptIDs{ii} '})$' ], 'Interpreter', 'Latex')
%         end
%         if ii == 5 || ii == 6
%             xlabel('time, $t_i$', 'Interpreter', 'latex')
%         end
%         xlim([0, Inf])
%         ylim([0, Inf]) 
%     end
% 
%     % Legend
%     subplot(round(length(ttc)*0.5)+1, 2, length(ttc) + 1)
%     labels = {} ;
%     for ii = 1:length(ttc)
%         plot([-1],[-1],'.-', 'color', colorset(ii, :))
%         labels{ii} = num2str(ii) ;
%         hold on
%     end
%     legend(labels, 'orientation', 'horizontal')
%     xlim([0, 1])
%     ylim([0, 1])
%     axis off
%     set(gcf, 'Units', 'centimeters');
%     set(gcf, 'Position', [0 0 16 30]) ;
%     set(gca,'fontsize', 12);
%     
%     if use_offset
%         figurefn = fullfile(timelineDir, 'time_correspondences_offset.pdf') ;
%     else
%         figurefn = fullfile(timelineDir, 'time_correspondences.pdf') ;
%     end
%     
%     saveas(gcf, figurefn)
%     
%     close all
% end

%% XIIb. Save timestamps in the original embryo folders
if saveResults
    maxv = -Inf ;
    clf
    figure('Units', 'centimeters', 'position', [0,0,6,6])
    for ii = 1:length(expts)
        tjr = i_tau0j_tau0jrelaxed(i_tau0j_tau0jrelaxed(:, 1) == ii, 3) ;
        tjr0 = i_tau0j_tau0jrelaxed(i_tau0j_tau0jrelaxed(:, 1) == ii, 2) ;

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
    saveas(gcf, fullfile(timelineDir, sprintf('relaxation_results_linearplot.pdf')))

end


