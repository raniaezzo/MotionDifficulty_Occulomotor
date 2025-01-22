function [lowerBound, upperBound] = findCI(bootstrapStatistics, percentileRange)

    if percentileRange == 68
        lPerc = 16; uPerc = 84;
    elseif percentileRange == 95
        lPerc = 2.5; uPerc = 97.5;
    end

    % Calculate confidence intervals for each time point
    lowerBound = prctile(bootstrapStatistics, lPerc, 1);
    upperBound = prctile(bootstrapStatistics, uPerc, 1);

end