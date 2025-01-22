function [x_fit, y_fit] = fitLine(RT_outputCond, postOnsetMS_cond)

    % Fit a line to the data
    valid_indices = ~isnan(RT_outputCond) & ~isnan(postOnsetMS_cond'); % Filter out NaN values
    x_valid = RT_outputCond(valid_indices);
    y_valid = postOnsetMS_cond(valid_indices)';
    
    % Fit a line to the valid data
    p = polyfit(x_valid, y_valid, 1);
    x_fit = linspace(min(x_valid), max(x_valid), 100);
    y_fit = polyval(p, x_fit);

end