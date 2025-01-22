function [sig_cluster, p_values_for_clusters] = permtestBinary(condition1, condition2, threshold, n_permutations)

    % INITIALIZE
    % Set optional arguments:
    if nargin < 4 || isempty(n_permutations)
        n_permutations = 1000; %10^4;
    end
    if nargin < 3 || isempty(threshold)
        threshold = 0.05;
    end
    if nargin < 2
        err
    end
    
    % Step 1: Compute the observed proportion of microsaccade starts
    rate_condition1 = (sum(condition1, 1) / size(condition1, 1))*1000; % Proportion for condition1
    rate_condition2 = (sum(condition2, 1) / size(condition2, 1))*1000; % Proportion for condition2
    observed_diff = rate_condition1 - rate_condition2; % Observed difference in proportions
    
    % Step 2: Initialize permutation testing
    n_trials1 = size(condition1, 1); % Number of trials in condition1
    n_trials2 = size(condition2, 1); % Number of trials in condition2
    combined_trials = [condition1; condition2]; % Combine trials from both conditions
    
    % Create a matrix to hold all permuted differences
    permuted_diffs = zeros(n_permutations, size(condition1, 2));
    
    % Step 3: Generate null distribution by shuffling trials
    for i = 1:n_permutations
        % Step 3.1: Shuffle trial labels
        permuted_indices = randperm(n_trials1 + n_trials2);
        permuted_condition1 = combined_trials(permuted_indices(1:n_trials1), :);
        permuted_condition2 = combined_trials(permuted_indices(n_trials1+1:end), :);
        
        % Step 3.2: Compute the permuted proportion difference for each time bin
        permuted_rate1 = (sum(permuted_condition1, 1) / n_trials1)*1000;
        permuted_rate2 = (sum(permuted_condition2, 1) / n_trials2)*1000;
        permuted_diffs(i, :) = permuted_rate1 - permuted_rate2;
    end
    
    % Step 4: Compare observed differences to the permuted distribution
    % For each time bin, calculate the proportion of permuted differences that are more extreme than the observed difference
    p_values = zeros(1, size(condition1, 2)); % Preallocate p-values vector
    for t = 1:size(condition1, 2)
        % Compute two-tailed p-values (this is the proportion that permuted
        % differences is larger than observed differences (small is good)
        p_values(t) = mean(abs(permuted_diffs(:, t)) >= abs(observed_diff(t)));
    end
    
    % Step 5: Form clusters of adjacent significant time bins based on the p-values
    % Here we use a significance threshold (e.g., p < 0.05) to define significant time bins
    significant_bins = p_values < threshold; % Identify significant time bins
    
    % Step 6: Form clusters from adjacent significant bins
    clusters = bwlabel(significant_bins); % Label contiguous significant bins as clusters
    
    % Step 7: Compute observed cluster statistics (sum of differences within each cluster)
    n_clusters = max(clusters); % Number of clusters
    observed_cluster_stats = zeros(1, n_clusters); % Preallocate
    for c = 1:n_clusters
        observed_cluster_stats(c) = sum(observed_diff(clusters == c)); % Sum of differences within each cluster
    end
    
    % Step 8: Compute permuted cluster statistics for each permutation
    permuted_cluster_stats = []; % Initialize empty array for storing permuted cluster stats
    for i = 1:n_permutations
        perm_significant_bins = abs(permuted_diffs(i, :)) >= quantile(abs(permuted_diffs(:)), 0.95); % Use 95th percentile of null distribution
        perm_clusters = bwlabel(perm_significant_bins); % Identify clusters in permuted data
        
        n_perm_clusters = max(perm_clusters); % Number of clusters in permuted data
        for c = 1:n_perm_clusters
            permuted_cluster_stats = [permuted_cluster_stats, sum(permuted_diffs(i, perm_clusters == c))]; %#ok<AGROW>
        end
    end
    
    % Step 9: Compare observed cluster statistics to the null distribution
    p_values_for_clusters = zeros(1, n_clusters);
    for c = 1:n_clusters
        p_values_for_clusters(c) = mean(permuted_cluster_stats >= observed_cluster_stats(c));
    end
    
    % Step 10: Identify significant clusters (e.g., p-value < 0.05)
    significant_clusters = find(p_values_for_clusters < 0.05);
    
    sigLabels = unique(significant_clusters);
    sig_cluster = {[], []};
    
    for cc=1:length(sigLabels)
        sig_cluster{cc} = find(clusters==significant_clusters(cc));
    end

end