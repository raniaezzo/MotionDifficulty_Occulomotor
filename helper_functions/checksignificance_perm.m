function [sig_cluster,sig_vals] = checksignificance_perm(rateSummary, dependent_samples)
        windowSize = 10; 
        
        % for permutation testing
        p_threshold = 0.05;
        num_permutations = 256; %1000; 
        two_sided = false;
        num_clusters =[];

        fields = fieldnames(rateSummary);
        fields{1}
        fields{2}
        data1 = rateSummary.(fields{1});
        data2 = rateSummary.(fields{2});
    
        % Apply the moving average filter to each row separately
        smoothedData = zeros(size(data1));
        smoothedData2 = zeros(size(data2));

        for i = 1:size(data1, 1)
            smoothedData1(i, :) = movmean(data1(i, :), windowSize);
        end
        for i = 1:size(data2, 1)
            smoothedData2(i, :) = movmean(data2(i, :), windowSize);
        end
    
        [clusters, p_values, ~, ~ ] = permutest( smoothedData1', smoothedData2', dependent_samples, ...
        p_threshold, num_permutations, two_sided, num_clusters );
        sig_cluster = clusters(p_values<p_threshold); % find significant cluster(s)
        sig_vals = p_values(p_values<p_threshold);

end