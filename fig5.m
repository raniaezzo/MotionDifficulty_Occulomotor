% combining timeseries plots
clc; clear all; %close all;

jsonparamfile = 'setup.json';
jsonParams = jsondecode(fileread(jsonparamfile));
subjects = jsonParams.subjects.Subjectids; 
protocols = jsonParams.protocols.Protocolids; 
datadir = jsonParams.datadir.Path; 
savedir = jsonParams.figsavedir.Path; 
nSampleCutOff = 4000; 

gausWindowSize = 100;
x = (1 : nSampleCutOff) - 1300;

trialwise = 1;
percentileRange = 68; % 68% CIs
subjectids = 1:8;
analyses = {'direction', 'tilt'};

for aa=1:length(analyses)
    analysis_type = analyses{aa}; %'direction'; %'direction'; % 'tilt'; %'location'
    
    if strcmp(analysis_type, 'direction')
        fieldNames = {'cardinal', 'oblique'};
        color = {[17, 119, 51],[51, 34, 136]};
        subfig = 'A';
    elseif strcmp(analysis_type, 'tilt')
        fieldNames = {'largeoffset','smalloffset'};
        color = {[0, 0, 0],[175, 175, 175]}; 
        subfig = 'B';
    elseif strcmp(analysis_type, 'location')
        fieldNames = {'horizontalLoc','verticalLoc'};
        color = {[0, 0, 0],[175, 175, 175]}; 
        subfig = 'C';
    elseif strcmp(analysis_type, 'dirtilt')
        fieldNames = {'easycardinal','hardoblique'};
        color = {[0, 0, 0],[175, 175, 175]}; 
        subfig = 'D';
    end
    
    combinedTimeseries = nan(length(subjects), length(x(1000:end)), length(fieldNames));
    combinedTimeseriesMeta = nan(length(x(1000:end)), length(fieldNames));
    metaRate = struct();
    
    for ii = subjectids
    
        summaryMSPath = fullfile(datadir, subjects{ii}, 'ProcessedData', 'Summary', 'microsaccades');
        fileName = sprintf('%s_%s_%s_%s_%i_%i', subjects{ii}, analysis_type, fieldNames{1}, fieldNames{2}, trialwise, percentileRange);
        
        summaryName = fullfile(summaryMSPath, fileName);
        
        load(strcat(summaryName, '.mat'))
        
        for pp=1:length(fieldNames)
            fieldName = fieldNames{pp};
            currRate = raster2rate(rate.(fieldName), gausWindowSize); % convert all trials for a condition to RATE (per subject)
            combinedTimeseries(ii, 1:length(currRate), pp) = currRate;   
            
            if ~isfield(metaRate, fieldName)
                metaRate.(fieldName) = rate.(fieldName);
            end
            
            % for meta subject
            metaRate.(fieldName) = [metaRate.(fieldName); rate.(fieldName)];
        end
        
    end
    
    for pp=1:length(fieldNames)
        fieldName = fieldNames{pp};
        % convert metasubject to RATE
        currRate = raster2rate(metaRate.(fieldName), gausWindowSize); % convert all trials for a condition to RATE (per subject)
        combinedTimeseriesMeta(1:length(currRate), pp) = currRate;  
    end
    
    %%
    rng(1)
    
    CI=1; % bootstrap
    metasubject = 0; %1;
    
    if metasubject
        meanTimeseries = combinedTimeseriesMeta;
    else
        meanTimeseries = squeeze(median(combinedTimeseries(subjectids,:,:), 1));
    end
    
    % plot a timeseries to check
    figure
    
    summaryFigPath = fullfile(datadir, 'ALLSUBJECTS', 'figures');

    if ~isfolder(summaryFigPath)
        mkdir(summaryFigPath)
    end

    figName = sprintf('%s_%s_%s_%s_%i_%i', 'allsubjects', analysis_type, fieldNames{1}, fieldNames{2}, trialwise, percentileRange);
    figPath = fullfile(summaryFigPath, figName);
    groupSummary = struct(); pairwiseSummary = struct();
    
    
    summaryMSPath = fullfile(datadir, 'ALLSUBJECTS','microsaccades');

    if ~isfolder(summaryMSPath)
        mkdir(summaryMSPath)
    end

    fileName = sprintf('%s_%s_%s_%s_%i_%i', 'ALLSUBJECTS', analysis_type, fieldNames{1}, fieldNames{2}, trialwise, percentileRange);
    summaryName = fullfile(summaryMSPath, fileName);
    
    for pp=1:length(fieldNames)
        fieldName = fieldNames{pp};
        meanTimeseries_cond = meanTimeseries(:,pp);
        a = plot(x(1001:1950), meanTimeseries_cond(1:950), '-','LineWidth',2, 'Color', color{pp}/255);
        hold on
        
        if metasubject
            groupSummary.(fieldName) = metaRate.(fieldName);
            CI = 1; convert2rate = 1;
        else
            groupSummary.(fieldName) = combinedTimeseries(subjectids,1:950,pp);
            convert2rate = 0; % already rate
        end
        
        pairwiseSummary.(fieldName) = combinedTimeseries(subjectids,1:950,pp);
        
        % confidence intervals
        if CI
            bootFile = strcat(summaryName, sprintf('%sbootstraps.mat', fieldName));
            if isfile(bootFile)
                load(bootFile)
            else
                % bootrap across subjects
                disp('Bootstrapping Data .. ')
                bootstrapStatistics = bootstrapData(rate, fieldName, 1); % 1 to convert to rate
                save(strcat(summaryName, sprintf('%sbootstraps.mat', fieldName)), 'bootstrapStatistics');
            end
            
            [lowerBound, upperBound] = findCI(bootstrapStatistics, percentileRange);
            upperError = upperBound - meanTimeseries_cond';
            lowerError = meanTimeseries_cond'-lowerBound;
            %a = shadedErrorBar(x(1001:1950), meanTimeseries_cond(1:950), [upperError(301:1250);lowerError(301:1250)], 'lineprops', {'-','Color',color{pp}/255},'transparent',1,'patchSaturation',0.2);
            hold on
        else
            % SEM
            semVal = std(groupSummary.(fieldName),0,1)/sqrt(size(groupSummary.(fieldName),1));
            %a = shadedErrorBar(x(1001:1950), meanTimeseries_cond(1:950),semVal(1:950), 'lineprops', {'-','Color',color{pp}/255},'transparent',1,'patchSaturation',0.2);
        end
       
        hold on
    end
    
    x_patch = [0, 500, 500, 0]; % x-coordinates
    y_patch = [0, 0, 2, 2]; % y-coordinates
    fill(x_patch, y_patch, [173, 216, 230]/255, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    hold on
    
    ylim([0,2])
    xlim([-300,650]); %2000])
    ylabel('microsaccade rate (hz)'); %, 'FontSize', 20)
    xlabel('time (ms)')
    f = gcf;
    f.Position = [5 996 1431 311];
    
    yLimits = ylim; % get current ylim
    
    
    % plot significant clusters (this method only when doing trialwise)
    if length(fieldNames) == 2 %&& ~trialwise % do cluster permutation
        % false because we are not directly comparing sessions (sometimes different # of sessions) 
        [sig_cluster,sig_vals] = checksignificance_perm(pairwiseSummary, 'true');
    end
    
    for si=1:length(sig_cluster)
        segLength = length(sig_cluster{si}(1:end));
        plot(sig_cluster{si}(1:end)-300, ones(1,segLength).*(yLimits(2)*.85), 'LineWidth',2, 'Color', 'k'); %color{1}/255);
        hold on % plot multiple clusters
    
        % compute cohen's d:
        % 1) take the pairwise difference
        effectDiff = pairwiseSummary.(fieldNames{1}) - pairwiseSummary.(fieldNames{2});
        deltas = effectDiff(:, sig_cluster{si});
        meanDeltas = mean(deltas,1);
        sdDeltas = nanstd(deltas);  % standard deviation of the diffferences
        s        = sdDeltas;        % re-name
        d_timepoint        =  meanDeltas ./ s;   % Cohen's d (paired version) per timepoint
        d = mean(d_timepoint)
    end
    

    f1 = gcf;
    f1.Position = [5 916 1469 391];
    dpi = get(0, 'ScreenPixelsPerInch');
    width_inch = 1469 / dpi;
    height_inch = 391 / dpi;
    
    set(gca, 'LineWidth', 1.5);
    dpi = get(0, 'ScreenPixelsPerInch');
    set(gca, 'FontName', 'Arial', 'FontSize', (391/dpi)*4);
    
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
    set(gcf, 'PaperPosition', [0 0 11.69 (11.69/width_inch)*height_inch]);
    print(gcf, '-dpdf', fullfile(savedir,sprintf('fig5%s.pdf', subfig)));  % for PDF
    close all;
end