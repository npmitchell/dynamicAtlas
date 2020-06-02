classdef queriedSample < handle
    %QUERIEDSAMPLE Object for manipulating subset of atlas
    %   Class for manipulating a queried subset of the dynamicAtlas
    %
    % Example
    % -------
    % a = dynamicAtlas('./', {'WT'}) 
    % qs = a.findLabel('Eve')
    % qs.getData()
    %
    % NPMitchell 2020
    
    properties
        meta        % struct with fields
        data        
    end
    
    
    methods
        function obj = queriedSample(sampleStruct)
            obj.meta = sampleStruct ;
        end
        
        function dataCell = getData(obj)
            %GETDATA(obj) Load the TIFFs for all data in this queriedSample
            %
            if isempty(obj.data)
                dataCell = cell(1, length(obj.meta.folders)) ;
                for qq = 1:length(obj.meta.folders)
                    dataCell{qq} = loadtiff(obj.meta.folders{qq}, obj.meta.names{qq}) ;
                end
                obj.data = dataCell ;
            else
                % If we have requested the data to be returned as a cell,
                % grab that cell 
                if nargout > 0
                    dataCell = obj.data ;
                end
            end
        end
        
    end
    
end

