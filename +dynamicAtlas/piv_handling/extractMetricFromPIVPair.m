function [MetricMat, t_min, t_max] = extractMetricFromPIVPair(VX, VY, WX, WY, piv_computation_method, truncation_dim_img, truncation_factor, V_RMS, W_RMS)
% extractMetricFromPIVPair(VX, VY, WX, WY, computation_method, truncation_factor)
%   Extract a comparison metric from computation on two 2D PIV fields V,W
%   returning a matrix of metric values calculated at all points in 
%   the field. VX, VY, WX, WY are assumed to have the same dimensions.
%   The field will be spatially truncated if truncation_factor is not 0.
%   The method of computation used is specified by piv_computation_method.
%   The dimension to truncate is given by truncation_dim (1 or 2)
%   The factor to truncate the image by, on both sides, is given by 
%   truncation_factor. The default dimension is 2, the horizontal 
%   axis of the image i.e. the vertical axis of the PIV field (the 
%   PIV program transposes the image when computing PIV to be symmetric
%   about the vertical axis). Also returns the minimum and maximum
%   time indices used along the truncation dimension, the matrix values
%   outside of that range should all be zeros.
%
%   Modified on 8/31 to input V_RMS and W_RMS as the RMS values for the 
%   given fields, where smoothed and truncated RMS values can be taken 
%   into account rather than simply the standard calculation.
%   
% Parameters
% ----------
% VX : X component of first PIV field (V), 2D matrix
% VY : Y component of first PIV field (V), 2D matrix
% WX : X component of second PIV field (W), 2D matrix
% WY : Y component of second PIV field (W), 2D matrix
%
% computation_method : string with potential values
%   normalized residual
%   unnormalized residual
%   inner product
%
% truncation_dim_img: dimension to truncate over, int with value 1 or 2,
% default dimension 2              
%
% truncation_factor: double factor to truncate matrix by on both sides, 
% must be between 0 (inclusive) and 0.5 (exclusive)
%
% Returns
%   MetricMat, a 2D matrix of correlation values calculated based on the 
%   chosen method.
%
% Outputs
% -------
%
% Vishank Jain-Sharma 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Verifies that the matrix width and length dimensions are all the same size

%verifies width
if size(VX,1) == size(VY,1) && size(VY,1) == size(WX,1) && size(WX,1) == size(WY,1)
    matrixHeight = size(VX,1);
else
    disp('Error: inputted matrix heights not all the same size!')
    matrixHeight = NaN;
end

%verifies length
if size(VX,2) == size(VY,2) && size(VY,2) == size(WX,2) && size(WX,2) == size(WY,2)
    matrixLength = size(VX,2);
else
    disp('Error: inputted matrix lengths not all the same size!')
    matrixLength = NaN;
end


%verifies truncation factor within range, defaults to no truncation if not
if (truncation_factor < 0 || truncation_factor >= 0.5)
    disp('Truncation factor must be between 0 inclusive and 0.5 exclusive!')
    disp('Defaulting to 0.')
    truncation_factor = 0;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%sets truncation indices based on dimension, default dimension 2
%this is horizontal axis of the image, which is vertical axis of PIV field 
%in convention with horizontal pullback AP-axis, this truncates the AP axis

%non-default dim of 1, the vertical axis of image, horizontal axis of PIV
if (truncation_dim_img == 1)
    
    %determine truncation indices, make sure within array index boundaries
    t_min = max(1, ceil(matrixLength*truncation_factor));
    t_max = min(matrixLength, floor(matrixLength*(1-truncation_factor)));
    
else
    
    %default truncation dim is 2, horiz axis of image, vert axis of PIV
    if (truncation_dim_img ~= 2)
        disp('Truncation dimension not set to 1 or 2! Defaulting to 2 (horizontal axis).')
    end
    
    %determine truncation indices, make sure within array index boundaries
    t_min = max(1, ceil(matrixHeight*truncation_factor));
    t_max = min(matrixHeight, floor(matrixHeight*(1-truncation_factor)));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%uses the residual described in Streichan 2018 E-Life paper, normalized 
%using the global field root-mean-square values
if strcmpi(piv_computation_method, 'normalized_residual')
    
    %RMS values of V and W matrices (when their X,Y components combined)
    %V_RMS = sqrt(mean(VX(:).^2) + mean(VY(:).^2));
    %W_RMS = sqrt(mean(WX(:).^2) + mean(WY(:).^2));  
    
    %residual computed and returns a matrix,
    %terms broken into pieces for easier readability
    %suffix indicates normalizing by global RMS value
    quadraticTerm_V_Norm = (VX.^2 + VY.^2)   / (V_RMS^2);
    quadraticTerm_W_Norm = (WX.^2 + WY.^2)   / (W_RMS^2);
    crossTerm_Norm       = (VX.*WX + VY.*WY) / (V_RMS * W_RMS);
    
    %returned matrix of correlation values
    MetricMatUntruncated = (quadraticTerm_V_Norm - 2 * crossTerm_Norm + quadraticTerm_W_Norm) / 2;
     
%uses the same residual as above but without normalizing for the global 
%RMS values of the vector fields
elseif strcmpi(piv_computation_method, 'unnormalized_residual')
    
    %residual computed and returns a matrix,
    %terms broken into pieces for easier readability
    %no global RMS normalization
    quadraticTerm_V = (VX.^2 + VY.^2);
    quadraticTerm_W = (WX.^2 + WY.^2);
    crossTerm       = (VX.*WX + VY.*WY);
    
    %returned matrix of correlation values
    MetricMatUntruncated = (quadraticTerm_V - 2 * crossTerm + quadraticTerm_W) / 2;
    
%simply takes inner product of V and W as the returned matrix
elseif strcmpi(piv_computation_method, 'inner_product')
    
    %inner product between V and W
    MetricMatUntruncated = VX.*WX + VY.*WY;
    
%if method not one of the above, defaults to inner product method
else
    disp('None of the ''normalized_residual'', ''unnormalized_residual'', or ''inner_product'' methods selected for comparison')
    disp('Defaulting to ''inner_product'' method')
    
    %inner product between V and W
    MetricMatUntruncated = VX.*WX + VY.*WY;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%truncates final output matrix according to the indices required, by
%setting all the values outside the range to 0

%starts from baseline of all zeros
MetricMatTemp = zeros(matrixHeight, matrixLength);

%truncation dim 1 means horizontal PIV axis to be truncated 
if (truncation_dim_img == 1) 
    
    %sets equal to computed values only in non-truncated range
    MetricMatTemp(:, t_min:t_max) = MetricMatUntruncated(:, t_min:t_max);
    
    MetricMat = MetricMatTemp;
    
%truncation dim 2 means vertical PIV axis to be truncated
else
    
    %sets equal to computed values only in non-truncated range
    MetricMatTemp(t_min:t_max, :) = MetricMatUntruncated(t_min:t_max, :);
    
    MetricMat = MetricMatTemp;    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
