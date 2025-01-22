function [bootstrapStatistics] = bootstrapData(rate, fieldName, convert)

    % Number of bootstrap samples
    numBootstrapSamples = 1000;
    
    % Preallocate array for storing statistics for each bootstrap sample
    bootstrapStatistics = zeros(numBootstrapSamples, size(rate.(fieldName), 2));
    
    gausWindowSize = 100;
    
    if convert
        % Bootstrap resampling: takes raw trial data and computes rate
        % iteratively 
        
        parfor i = 1:numBootstrapSamples
        
            % Generate a bootstrap sample by sampling rows with replacement
            bootstrapSampleIndices = randi(size(rate.(fieldName), 1), size(rate.(fieldName), 1), 1);
            bootstrapSample = rate.(fieldName)(bootstrapSampleIndices, :);
            
            bootstrapStatistic = raster2rate(bootstrapSample, gausWindowSize);
        
            % Store the statistic for this bootstrap sample
            bootstrapStatistics(i, :) = bootstrapStatistic;
        end
    else

        % assumes units are in rate (hz) already (e.g., bootstrapping
        % across subjects)
        
        % Number of rows to sample
        numSamples = numBootstrapSamples;

        % Generate random row indices with replacement
        randomIndices = randi(size(rate.(fieldName), 1), numSamples, 1);

        % Sample rows from the matrix
        bootstrapStatistics = rate.(fieldName)(randomIndices, :);

    end
    
    
end