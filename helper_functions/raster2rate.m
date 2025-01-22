function smoothedRateHz = raster2rate(raster, gausWindowSize)

    % Number of trials (rows)
    numTrials = size(raster, 1);

    % Initialize the count of first occurrences
    firstOccurrenceCount = zeros(1, size(raster, 2));

    % Iterate through each trial to find the first occurrences
    for trial = 1:numTrials
        for timepoint = 2:size(raster, 2)
            % Check if the current value is 1 and the previous value is not 1
            if raster(trial, timepoint) == 1 && (isnan(raster(trial, timepoint - 1)) || raster(trial, timepoint - 1) ~= 1)
                firstOccurrenceCount(timepoint) = firstOccurrenceCount(timepoint) + 1;
            end
        end
        % Check the first column separately
        if raster(trial, 1) == 1
            firstOccurrenceCount(1) = firstOccurrenceCount(1) + 1;
        end
    end

    % Compute the rate of first occurrences per millisecond (proportion of trials with first occurrence)
    % instead of # of trials, consider only the trials with valid time bins
    firstOccurrenceRatePerMs = firstOccurrenceCount ./ sum(~isnan(raster),1); % numTrials;

    % OR TO COUNT EVERY BIN
    %firstOccurrenceRatePerMs = sum(raster,1) ./ sum(~isnan(raster),1);  %

    % Convert rate to Hertz (events per second)
    firstOccurrenceRateHz = firstOccurrenceRatePerMs * 1000;

    % Smooth the rate using a Gaussian sliding window of 100 ms
    smoothedRateHz = smoothdata(firstOccurrenceRateHz, 'gaussian', gausWindowSize);

end