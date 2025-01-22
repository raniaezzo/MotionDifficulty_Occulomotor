clc; clear all; %close all;

jsonparamfile = 'setup.json';
jsonParams = jsondecode(fileread(jsonparamfile));
subjects = jsonParams.subjects.Subjectids; 
protocols = jsonParams.protocols.Protocolids; 
datadir = jsonParams.datadir.Path; 
savedir = jsonParams.figsavedir.Path; 

analysis_type = 'direction';

stimonset = 1300; % this data has 1000 ms cutoff from beginning
stimoffset = 1800; 

analyses = {'direction', 'tilt', 'outcome'};
x_values = linspace(-5, 5, 100);
z_score_init = 60000; % arbitrary, but never goes beyond

msOnset = nan(length(subjects), 2, length(analyses)); % subjects, easy/hard, analysis
zscored_msOnset = nan(length(subjects), 2, length(analyses));
plotOn = 0;

% Initialize the cell array for probability density
probDensity = cell(1, length(analyses));
grand_z_scores_MSonset = cell(1, length(analyses)); 
grand_possible = cell(1, length(analyses)); 

amp = cell(1, length(analyses));
boxamp = cell(1, length(analyses));
latvals = cell(1, length(analyses));
eventtimes = cell(1, length(analyses));

for cc = 1:length(analyses)
    matrix = nan(length(subjects), length(x_values), 2);
    matrix2 = nan(length(subjects), z_score_init, 2);
    matrix3 = nan(length(subjects), 1, 2);
    matrix4 = nan(length(subjects),3,2);
    % Assign the matrix to the cell
    probDensity{cc} = matrix;
    grand_z_scores_MSonset{cc} = matrix2; % mat2 is for zscores
    grand_possible{cc} = matrix3; % mat2 is for zscores
    amp{cc} = matrix4;
    boxamp{cc} = matrix3; % only 1 val
    latvals{cc} = matrix4;
    eventtimes{cc} = matrix4;
end

totalQualifiedTrials = nan(length(subjects), length(analyses));

grand_z_scores_MSonset_easy = [];
grand_z_scores_MSonset_hard = [];
grand_possible_easy = [];
grand_possible_hard = [];

rtConds = nan(length(subjects), 3, length(analyses)); % 3 time measures saved

if plotOn ~= 0
    %f1 = figure(1);
    f2 = figure(2);
    f3 = figure(3);
    f4 = figure(4);
    f5 = figure(5);
    f6 = figure(6);
end

f6.Position = [135 739 1848 400];

% get the scale factor for PSC
psc_scaleFactor = nan(length(subjects), 1);
    
%%
for si=1:length(subjects)

    subj=subjects{si};

    savePath = sprintf('%s/%s/ProcessedData/Summary/pupil', datadir, subj);
    load(fullfile(savePath,sprintf('%s_allpupilFits.mat', subj)), 'output');
    load(fullfile(savePath,sprintf('%s_allpupilData.mat', subj))); % to load in tab file

    psc_scaleFactor(si,1) = (nanmean(alltrialSignalSummary(:,2))/nanmax(alltrialSignalSummary(:,4)));


    rate = load(fullfile(strrep(savePath, 'pupil', 'microsaccades'), sprintf('%s_allmsData.mat', subj)));

    % S07 has VL file missing
    if si==7
         alltrialTimeStamps = alltrialTimeStamps(alltab(:,10)~=7,:);
         alltrialSignalSummary = alltrialSignalSummary(alltab(:,10)~=7,:);
         allfilteredPupilData = allfilteredPupilData(alltab(:,10)~=7,:);
         alltab = alltab(alltab(:,10)~=7,:);
         rate.allMSData = rate.allMSData(alltab(:,10)~=7,:);
    end

    % Correction for some trials with <500 RT: some early sessions counted
    % ms relative to stimOff rather than stimOn
    temp = alltab(:,1); % get the trial number
    next_session = [false; diff(temp) < 0];     % Find the points where the numbers start to decrease
    indices = cumsum(next_session) + 1;
    rtTrialwise.allMSData = nan(length(alltab),1);
    suppressionTimes = nan(length(alltab),1);
    postsupPeakTimes = nan(length(alltab),1);

    % to find the post stim inflection point
    window_size = 250; % Window size for the moving average

    for sn=1:max(indices) % this is by session
        [selectTrials, ~] = find(indices==sn);
    
%         % find the suppression for the mean rate over session
        meanMS = nanmean(rate.allMSData(selectTrials,:),1);

        smoothed_data = smoothdata(meanMS, 'movmean', window_size);
            
        decreasing_index = find(diff(smoothed_data(1800:end)) < 0, 1);

        if isempty(decreasing_index)
            decreasing_index = find(diff(smoothed_data(1800:end)) == 0, 1); % if the MS rate does not go back down, select when it plateaus
        end

        postsupPeakTimes(selectTrials) = ones(length(selectTrials),1)*(decreasing_index+1800);

        [~,suppressionTime] = min(meanMS(stimonset:stimoffset)); 
        suppressionTimes(selectTrials) = ones(length(selectTrials),1)*suppressionTime; % suppression computed within session (but this doesn't matter)

        RT_unfiltered = alltab(selectTrials,13)*1000; % all RTs, no filter
        if any(RT_unfiltered<500)
            rtTrialwise.allMSData(selectTrials) = RT_unfiltered +500;
        else
            rtTrialwise.allMSData(selectTrials) = RT_unfiltered;
        end
    end

    rtTrialwise.allMSData(rtTrialwise.allMSData>=2000) = nan;

    trialData.allMSData = alltab;

    [postOnsetMS_cond, latencyCond, RT_outputCond, accuracy_Cond, qualifyingTrials, dataCheck, supTimes, postOffsetTimes, selectedTrials] = ...
        compute_latency(rate, 'allMSData', suppressionTimes, postsupPeakTimes, rtTrialwise, trialData, stimonset);

    tabnew = qualifyingTrials{1}; % 1 b/c all conditions are in cell 1

    pupilFitsMatch = output(selectedTrials);
    pupilData = allfilteredPupilData(selectedTrials,:);
    
    for ai=1:length(analyses)

        analysis_type = analyses{ai};

        % Correalte MS latency with Pupil Latency
        if strcmp(analysis_type, 'direction')
            fieldNames = {'cardinal', 'oblique'};
            color = {[17, 119, 51],[51, 34, 136]}; 
        elseif strcmp(analysis_type, 'tilt')
            fieldNames = {'largeoffset','smalloffset'};
            color = {[0, 0, 0],[175, 175, 175]}; 
        elseif strcmp(analysis_type, 'location')
            fieldNames = {'horizontalLoc','verticalLoc'};
            color = {[0, 0, 0],[175, 175, 175]}; 
        elseif strcmp(analysis_type, 'dirtilt')
            fieldNames = {'easycardinal','hardoblique'};
            color = {[0, 0, 0],[175, 175, 175]}; 
        elseif strcmp(analysis_type, 'outcome')
            fieldNames = {'correct','incorrect'};
            color = {[0, 255, 0],[255, 0, 0]}; 
        end


        if strcmp(analysis_type, 'direction')
            easyTrials = ismember(tabnew(:,10), 5:8); % & tabnew(:,14)==1; 
            hardTrials = ismember(tabnew(:,10), 1:4); % & tabnew(:,14)==1;
            possibleEasy = sum(ismember(trialData.allMSData(:,10), 5:8));
            possibleHard = sum(ismember(trialData.allMSData(:,10), 1:4));
        elseif strcmp(analysis_type, 'tilt')
            easyTrials = tabnew(:,11)==8; %  & tabnew(:,14)==1;
            hardTrials = tabnew(:,11)<=4; %  & tabnew(:,14)==1;
            possibleEasy = sum(trialData.allMSData(:,11)==8); 
            possibleHard = sum(trialData.allMSData(:,11)<=4); %2);
        elseif strcmp(analysis_type, 'outcome')
            easyTrials = tabnew(:,14)==1;
            hardTrials = tabnew(:,14)==0;
            possibleEasy = sum(trialData.allMSData(:,14)==1);
            possibleHard = sum(trialData.allMSData(:,14)==0);
        end

        msOnset(si,1, ai) = nanmean(postOnsetMS_cond(easyTrials));
        msOnset(si,2, ai) = nanmean(postOnsetMS_cond(hardTrials));

        if plotOn ~= 0
            [x_fit, y_fit] = fitLine(RT_outputCond(easyTrials), postOnsetMS_cond(easyTrials));
            plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', color{1}/255);
        
            hold on
            [x_fit2, y_fit2] = fitLine(RT_outputCond(hardTrials), postOnsetMS_cond(hardTrials));
            plot(x_fit2, y_fit2, '-', 'LineWidth', 2, 'Color', color{2}/255);
    
            text(1800, 1200, num2str(correlation_coefficient1)); % Coordinates for the first line of text (2, 5)
            text(1800, 1100, num2str(correlation_coefficient2)); 
        end

        % recalculate the suppression for the hard vs. easy (based on data output)
        
        [postsupPeakTimes_easy, supTime_easy] = calculateSuppPeak(rate.allMSData, easyTrials, window_size);
        [postsupPeakTimes_hard, supTime_hard] = calculateSuppPeak(rate.allMSData, hardTrials, window_size);
        
        rtConds(si,1,ai) = nanmedian(RT_outputCond(easyTrials));
        rtConds(si,2,ai) = nanmedian(RT_outputCond(hardTrials));
        rtConds(si,3,ai) = nanmedian(RT_outputCond(hardTrials))-nanmedian(RT_outputCond(easyTrials));

        % pupil data
        pupilFit_easy = pupilFitsMatch(easyTrials);
        pupilFit_hard = pupilFitsMatch(hardTrials);

        varExp_easy = cell2mat({pupilFit_easy.('R2')});
        varExp_hard = cell2mat({pupilFit_hard.('R2')});

        % which values do I need to plot the signal (reconstruction)?
        pupilFit_easy_filt = pupilFit_easy(varExp_easy>0.50);
        pupilFit_hard_filt = pupilFit_hard(varExp_hard>0.50);

        amp{ai}(si,1:3,1) = median(cell2mat({pupilFit_easy_filt.('ampvals')}'), 1);
        boxamp{ai}(si,1,1) = median(cell2mat({pupilFit_easy_filt.('boxampvals')}'), 1);
        latvals{ai}(si,1:3,1) = median(cell2mat({pupilFit_easy_filt.('latvals')}'), 1);
        eventtimes{ai}(si,1:3,1) = median(cell2mat({pupilFit_easy_filt.('eventtimes')}'), 1);
        amp{ai}(si,1:3,2) = median(cell2mat({pupilFit_hard_filt.('ampvals')}'), 1);
        boxamp{ai}(si,1,2) = median(cell2mat({pupilFit_hard_filt.('boxampvals')}'), 1);
        latvals{ai}(si,1:3,2) = median(cell2mat({pupilFit_hard_filt.('latvals')}'), 1);
        eventtimes{ai}(si,1:3,2) = median(cell2mat({pupilFit_hard_filt.('eventtimes')}'), 1);

        avTotal = mean(postOnsetMS_cond);
        stdevTotal = std(postOnsetMS_cond);
        z_scored_data = (postOnsetMS_cond-1300);

        z_scores_MSonset_easy = z_scored_data(easyTrials);
        z_scores_MSonset_hard = z_scored_data(hardTrials);

        grand_z_scores_MSonset{ai}(si,1:length(z_scores_MSonset_easy),1) = z_scores_MSonset_easy;
        grand_z_scores_MSonset{ai}(si,1:length(z_scores_MSonset_hard),2) =  z_scores_MSonset_hard;
        grand_possible{ai}(si,1:length(possibleEasy),1) =  possibleEasy;
        grand_possible{ai}(si,1:length(possibleHard),2) =  possibleHard;

        zscored_msOnset(si,1, ai) = nanmedian(z_scores_MSonset_easy);
        zscored_msOnset(si,2, ai) = nanmedian(z_scores_MSonset_hard);


        [fnorm_easy, ~] = ksdensity(z_scores_MSonset_easy, x_values, 'function','cdf');
        [fnorm_hard, xinorm] = ksdensity(z_scores_MSonset_hard, x_values, 'function','cdf');

        f1_normalized = fnorm_easy * (numel(z_scores_MSonset_easy)/possibleEasy);
        f2_normalized = fnorm_hard * (numel(z_scores_MSonset_hard)/possibleHard);

        probDensity{ai}(si,1:length(x_values),1) = f1_normalized; 
        probDensity{ai}(si,1:length(x_values),2) = f2_normalized;

        totalQualifiedTrials(si, ai) = length(postOnsetMS_cond(easyTrials))+length(postOnsetMS_cond(hardTrials));

    end

end


%% Flatten grant array for easy/hard

normalizeData = 0; % this will make the probability distributions to mean "how probable, out of all trials of this trial type, would MS occur at or before X

probDelay = nan(length(subjects), 3, length(analyses)); % 3 time measures saved

% for shaded patch
x_patch = [0, 500, 500, 0]; % x-coordinates
y_patch = [0, 0, 1, 1]; % y-coordinates

%%
for si=1:length(subjects)
    subjectID = si
    subj = subjects{si};
    figure
    figName = sprintf('%s/%s/ProcessedData/Summary/figures/%s', datadir, subj, sprintf('%sprobDensityFunction', subj));
    for ai=1:2 %length(analyses)
    
        analysis_type = analyses{ai};
    
        if strcmp(analysis_type, 'direction')
            fieldNames = {'cardinal', 'oblique'};     
            color = {[17, 119, 51],[51, 34, 136]}; 
        elseif strcmp(analysis_type, 'tilt')
            fieldNames = {'largeoffset','smalloffset'};
            color = {[0, 0, 0],[175, 175, 175]}; 
        elseif strcmp(analysis_type, 'location')
            fieldNames = {'horizontalLoc','verticalLoc'};
            color = {[0, 0, 0],[175, 175, 175]}; 
        elseif strcmp(analysis_type, 'dirtilt')
            fieldNames = {'easycardinal','hardoblique'};
            color = {[0, 0, 0],[175, 175, 175]}; 
        elseif strcmp(analysis_type, 'outcome')
            fieldNames = {'correct','incorrect'};
            color = {[0, 255, 0],[255, 0, 0]}; 
        end
    
        easy_flat = grand_z_scores_MSonset{ai}(subjectID,:,1);
        hard_flat = grand_z_scores_MSonset{ai}(subjectID,:,2);
        
        easy_flat = easy_flat(:); easy_flat = easy_flat(~isnan(easy_flat));
        hard_flat = hard_flat(:); hard_flat = hard_flat(~isnan(hard_flat));
        
        possibleEasy = sum(grand_possible{ai}(subjectID,:,1));
        possibleHard = sum(grand_possible{ai}(subjectID,:,2));
        
        % randomly draw the # of values in each array with replacement
        % fit a kdensity function 1000 times
        x_values = linspace(0, 1500, 500); %linspace(-5, 5, 100);
        numBootstraps = 1000; %500;
        btSamples_easy = nan(numBootstraps, numel(x_values));
        btSamples_hard = nan(numBootstraps, numel(x_values));
        
        for bi=1:numBootstraps
        
        if normalizeData
            ratio2norm_easy = (numel(easy_flat))/possibleEasy;
            ratio2norm_hard = (numel(hard_flat))/possibleHard;
        else
            ratio2norm_easy = 1;
            ratio2norm_hard = 1;
        end
    
            if bi==1
                % mean values (only compute once)
                [fnorm_easy, ~] = ksdensity(easy_flat, x_values, 'function','cdf');
                [fnorm_hard, xinorm] = ksdensity(hard_flat, x_values, 'function','cdf');
                
                 f1_normalized = fnorm_easy * ratio2norm_easy; % this would tell how probably it is across exp
                 f2_normalized = fnorm_hard * ratio2norm_hard;
            end
        
            % bootstrap
            easy_flat_bt = datasample(easy_flat, numel(easy_flat), 'Replace', true);
            hard_flat_bt = datasample(hard_flat, numel(hard_flat), 'Replace', true);
            
            [fnorm_easy, ~] = ksdensity(easy_flat_bt, x_values, 'function','cdf');
            [fnorm_hard, xinorm] = ksdensity(hard_flat_bt, x_values, 'function','cdf');
            
             %f1_normalized 
             btSamples_easy(bi,1:numel(x_values)) = fnorm_easy * ratio2norm_easy; 
             %f2_normalized
             btSamples_hard(bi,1:numel(x_values)) = fnorm_hard * ratio2norm_hard;
        
        end
        
        % Compute the lower and upper percentiles for the 95% or 68% confidence interval
        lbd_easy = prctile(btSamples_easy, 2.5,1); %16, 1);
        ubd_easy = prctile(btSamples_easy, 97.5,1); %84, 1);
        lbd_hard = prctile(btSamples_hard, 2.5,1); %16, 1);
        ubd_hard = prctile(btSamples_hard, 97.5,1); %84, 1);
        
         subplot(1,2,ai)
         % Plot the square with a shaded fill
         fill(x_patch, y_patch, [173, 216, 230]/255, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
         hold on
         plot(xinorm, f1_normalized, 'Color', color{1}/255, 'LineWidth',1.5); % Plot the KDE
         hold on
         %errorbar(xinorm, f1_normalized, f1_normalized - lbd_easy, ubd_easy - f1_normalized, 'b', 'LineWidth', 1.5, 'LineStyle', 'none');
         %hold on
         shadedErrorBar(xinorm, f1_normalized, [ubd_easy-f1_normalized; f1_normalized-lbd_easy], 'lineprops', {'-','Color',color{1}/255},'transparent',1,'patchSaturation',0.1);
         hold on
         [half_val_1, half_index_1] = min(abs(f1_normalized - .5));
         time1 = xinorm(half_index_1);
         plot([time1, time1], [0, f1_normalized(half_index_1)], ':', 'Color',color{1}/255, 'LineWidth',4);
         hold on
         %plot([0, time1], [f1_normalized(half_index_1), f1_normalized(half_index_1)], ':', 'Color',color{1}/255, 'LineWidth',2);
         %hold on
         plot(xinorm, f2_normalized, 'Color', color{2}/255, 'LineWidth',1.5); % Plot the KDE
         hold on
         %errorbar(xinorm, f2_normalized, f2_normalized - lbd_hard, ubd_hard - f2_normalized, 'b', 'LineWidth', 1.5, 'LineStyle', 'none');
         %hold on
         shadedErrorBar(xinorm, f2_normalized, [ubd_hard-f2_normalized; f2_normalized-lbd_hard], 'lineprops', {'-','Color',color{2}/255},'transparent',1,'patchSaturation',0.2);
         hold on
         [half_val_2, half_index_2] = min(abs(f2_normalized - .5));
         time2 = xinorm(half_index_2);
         plot([time2, time2], [0, f2_normalized(half_index_2)], ':', 'Color',color{2}/255, 'LineWidth',4);
         hold on
         plot([0, time2], [f2_normalized(half_index_2), f2_normalized(half_index_2)], ':', 'Color','k', 'LineWidth',4);
         %title(analyses{ai})
         %xlim([-3 3])
         probDelay(si,1,ai) = time1;
         probDelay(si,2,ai) = time2;
         probDelay(si,3,ai) = time2-time1;
         %set(gca, 'FontName', 'Arial', 'FontSize', 12);
         
         dpi = get(0, 'ScreenPixelsPerInch');
         set(gca, 'FontName', 'Arial', 'FontSize', (365/dpi)*4.5);
         ax = gca;  % Get the current axis
        ax.XAxis.LineWidth = 1.5;  % Set the X-axis line width
        ax.YAxis.LineWidth = 1.5;  % Set the Y-axis line width

         ylabel('cumulative probability density')
         xlabel('time (ms)')
         axis square
         box off
    end

    f = gcf;
    set(f, 'Position', [94 972 893 365])

    width_inch = 893 / dpi;
    height_inch = 365 / dpi;

    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
    if si==3
        print(gcf, '-dpdf', fullfile(savedir,'fig6.pdf'));  % example shown in manuscript
    end
end

%%

ai= 1;
[nu, val, ci, stats] = ttest(probDelay(:,2,ai), probDelay(:,1,ai))
d = computeCohen_d(probDelay(:,2,ai), probDelay(:,1,ai), 'paired')

ai= 2;
[nu, val, ci, stats] = ttest(probDelay(:,2,ai), probDelay(:,1,ai))
d = computeCohen_d(probDelay(:,2,ai), probDelay(:,1,ai), 'paired')


%%
figure(1)

n = length(subjects);
ones_vector = zeros(n, 1);
jittered_vector = ones_vector + (rand(n, 1) - 0.5) * 0.55;

figName1 = sprintf('%s/ampPupilperCond', datadir);
figName2 = sprintf('%s/peakPupilperCond', datadir);
figName3 = sprintf('%s/ampDiffperCond', datadir);
figName4 = sprintf('%s/peakDiffperCond', datadir);

for ai=1:2 %2 %length(analyses)

    analysis_type = analyses{ai};

    if strcmp(analysis_type, 'direction')
        fieldNames = {'cardinal', 'oblique'};
        color = {[17, 119, 51],[51, 34, 136]}; 
        alphlev1 = 0.5;
        lincol1 = 'k';
        alphlev2 = 0.5;
        lincol2 = 'k';
        condComparison = '(obl - card)';
    elseif strcmp(analysis_type, 'tilt')
        fieldNames = {'largeoffset','smalloffset'};
        color = {[0, 0, 0],[175, 175, 175]}; 
        alphlev1 = 0.5;
        lincol1 = 'k';
        alphlev2 = 0.5;
        lincol2 = 'k';
        condComparison = '(small - large)';
    elseif strcmp(analysis_type, 'location')
        fieldNames = {'horizontalLoc','verticalLoc'};
        color = {[0, 0, 0],[175, 175, 175]}; 
    elseif strcmp(analysis_type, 'dirtilt')
        fieldNames = {'easycardinal','hardoblique'};
        color = {[0, 0, 0],[175, 175, 175]}; 
    elseif strcmp(analysis_type, 'outcome')
        fieldNames = {'correct','incorrect'};
        color = {[0, 255, 0],[255, 0, 0]}; 
    end

    figure(1)
    subplot(1,2,ai)
    [nsubs, ~] = size(amp{ai}(:,1,1));
    psc1 = ((amp{ai}(:,1,2).*psc_scaleFactor)-(amp{ai}(:,1,1).*psc_scaleFactor));
    [nu, val, ci, stats] = ttest(psc1)
    d = computeCohen_d(amp{ai}(:,1,2),amp{ai}(:,1,1), 'paired')
    hold on
    sed1 = nanstd(psc1) / sqrt(nsubs);
    scatter(jittered_vector+(1.5*ones(length(amp{ai}(:,1,1)),1)), psc1, 250, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.2) 
    hold on
    errorbar(1.5, nanmean(psc1), sed1, 'k-', 'LineWidth',5, 'CapSize', 15);
    psc2 = ((amp{ai}(:,2,2).*psc_scaleFactor)-(amp{ai}(:,2,1).*psc_scaleFactor));
    [nu, val, ci, stats] = ttest(psc2)
    d = computeCohen_d(amp{ai}(:,2,2),amp{ai}(:,2,1), 'paired')
    sed2 = nanstd(psc2) / sqrt(nsubs);
    hold on
    scatter(jittered_vector+(5.5*ones(length(amp{ai}(:,2,1)),1)), psc2, 250, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.2) 
    hold on
    errorbar(5.5, nanmean(psc2), sed2, 'k-', 'LineWidth',5, 'CapSize', 15);
    hold on   
    psc3 = ((amp{ai}(:,3,2).*psc_scaleFactor)-(amp{ai}(:,3,1).*psc_scaleFactor));
    [nu, val, ci, stats] = ttest(psc3)
    d = computeCohen_d(amp{ai}(:,3,2),amp{ai}(:,3,1), 'paired')
    sed3 = nanstd(psc3) / sqrt(nsubs);
    scatter(jittered_vector+(9.5*ones(length(amp{ai}(:,3,1)),1)), psc3, 250, 'filled', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.2) 
    hold on
    errorbar(9.5, nanmean(psc3), sed3, 'k-', 'LineWidth',5, 'CapSize', 15);
    psc4 = ((boxamp{ai}(:,1,2).*psc_scaleFactor)-(boxamp{ai}(:,1,1).*psc_scaleFactor));
    [nu, val, ci, stats] = ttest(psc4)
    d = computeCohen_d(boxamp{ai}(:,1,2),boxamp{ai}(:,1,1), 'paired')
    sed4 = nanstd(psc4) / sqrt(nsubs);
    scatter(jittered_vector+(13.5*ones(length(boxamp{ai}(:,1,1)),1)), psc4, 250, 'filled', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'b', 'LineWidth', 2, 'MarkerFaceAlpha', 0.2)
    hold on
    errorbar(13.5, nanmean(psc4), sed4, 'k-', 'LineWidth',5, 'CapSize', 15);
    hold on 
    yline(0, ':k', 'LineWidth',2)
    xlim([0,15])
    ylim([-1,2])
    ax1 = gca;
    ax1.XTick = linspace(1.5, 13.5, 4);
    ylabel(['pupil area % ', condComparison]);
    set(gca, 'XTickLabelRotation', 45);
    ax1.XTickLabel = {"trial start"; "task-evoked"; "response-evoked"; "stimulus boxcar"};
    f3 = gcf;
    set(gca, 'LineWidth', 2);
    set(gca, 'FontName', 'Arial', 'FontSize', (911/dpi)*1.8);
    ax = gca;  % Get the current axis
    ax.XAxis.LineWidth = 1.5;  % Set the X-axis line width
    ax.YAxis.LineWidth = 1.5;  % Set the Y-axis line width
    box off

    figure(2)
    subplot(1,2,ai)
    [nu, val, ci, stats] = ttest((latvals{ai}(:,1,2)+eventtimes{ai}(:,1,2)+930)-(latvals{ai}(:,1,1)+eventtimes{ai}(:,1,1)+930))
    hold on
    scatter(jittered_vector+ones(length(latvals{ai}(:,1,1)),1), latvals{ai}(:,1,1)+eventtimes{ai}(:,1,1)+930, 250, 'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', lincol1, 'LineWidth', 2, 'MarkerFaceAlpha', alphlev1)
    hold on
    scatter(jittered_vector+(3*ones(length(latvals{ai}(:,1,2)),1)), latvals{ai}(:,1,2)+eventtimes{ai}(:,1,2)+930, 250, 'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', lincol2, 'LineWidth', 2, 'MarkerFaceAlpha', alphlev2)
    hold on
    plot([jittered_vector+1 jittered_vector+3]', [latvals{ai}(:,1,1)+eventtimes{ai}(:,1,1)+930, latvals{ai}(:,1,2)+eventtimes{ai}(:,1,2)+930]', 'k', 'LineWidth', 1)
    hold on
    [nu, val, ci, stats] = ttest((latvals{ai}(:,2,2)+eventtimes{ai}(:,2,2)+930)-(latvals{ai}(:,2,1)+eventtimes{ai}(:,2,1)+930))
    hold on
    scatter(jittered_vector+(5*ones(length(latvals{ai}(:,2,1)),1)), latvals{ai}(:,2,1)+eventtimes{ai}(:,2,1)+930, 250, 'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', lincol1, 'LineWidth', 2, 'MarkerFaceAlpha', alphlev1)
    hold on
    scatter(jittered_vector+(7*ones(length(latvals{ai}(:,2,2)),1)), latvals{ai}(:,2,2)+eventtimes{ai}(:,2,2)+930, 250, 'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', lincol2, 'LineWidth', 2, 'MarkerFaceAlpha', alphlev2)
    hold on
    plot([jittered_vector+5 jittered_vector+7]', [latvals{ai}(:,2,1)+eventtimes{ai}(:,2,1)+930, latvals{ai}(:,2,2)+eventtimes{ai}(:,2,2)+930]', 'k', 'LineWidth', 1)
    hold on
    [nu, val, ci, stats] = ttest((latvals{ai}(:,3,2)+eventtimes{ai}(:,3,2)+930)-(latvals{ai}(:,3,1)+eventtimes{ai}(:,3,1)+930))
    scatter(jittered_vector+(9*ones(length(latvals{ai}(:,3,1)),1)), latvals{ai}(:,3,1)+eventtimes{ai}(:,3,1)+930, 250, 'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', lincol1, 'LineWidth', 2, 'MarkerFaceAlpha', alphlev1)
    hold on
    scatter(jittered_vector+(11*ones(length(latvals{ai}(:,3,2)),1)), latvals{ai}(:,3,2)+eventtimes{ai}(:,3,2)+930, 250, 'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', lincol2, 'LineWidth', 2, 'MarkerFaceAlpha', alphlev2)
    hold on
    plot([jittered_vector+9 jittered_vector+11]', [latvals{ai}(:,3,1)+eventtimes{ai}(:,3,1)+930, latvals{ai}(:,3,2)+eventtimes{ai}(:,3,2)+930]', 'k', 'LineWidth', 1)
    ax1 = gca;
    ax1.XTick = linspace(2, 10, 3);
    ylabel('peak dilation (ms)')
    ax1.XTickLabel = {"trial start"; "task-evoked"; "response-evoked"};
    xlim([0,12])
    f2 = gcf;
    set(gca, 'XTickLabelRotation', 45);
    %set(gca, 'FontName', 'Arial', 'FontSize', 20);
    set(gca, 'LineWidth', 2);
    set(gca, 'FontName', 'Arial', 'FontSize', (911/dpi)*1.8);
    ax = gca;  % Get the current axis
    ax.XAxis.LineWidth = 1.5;  % Set the X-axis line width
    ax.YAxis.LineWidth = 1.5;  % Set the Y-axis line width
    box off
    axis square

end

for fi=[f1, f2,f3] %, f4]

    if fi==1
        figName = '9'; %figName3;
        f = figure(2);
        set(f, 'Position', [194 784 956 426]); %[67 547 1971 739])
    elseif fi ==2
        figName = '10'; figName2;
        f = figure(1);
        set(f, 'Position', [334 778 889 513]); %[190 632 1639 648])
    end


    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
    set(gcf, 'PaperPositionMode', 'auto'); 

    tmp = strsplit(figName, '/');

set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);

print(gcf, '-dpdf', fullfile(savedir, sprintf('fig%s.pdf', figName)));  % for PDF

end

%% %% PARTIAL CORRELATIONS

aestheticScalar = 1.3;
dpi = get(0, 'ScreenPixelsPerInch');
analyses = {'direction','tilt'};
observerID = [1:8, 1:8]; observerID = observerID';
condition = [2*ones(8,1); ones(8,1)];
event = 2; % ONLY FOR PUPIL - DOES NOT APPLY TO THIS FIRST SECTION

% REGRESSING OUT CONDITION EFFECT

% Partial correlation (RT - to - MS rebound)

for ai=1:length(analyses)

    analysis_type = analyses{ai};
    
    if ai==1
        fieldNames = {'cardinal', 'oblique'};
        color = {[51, 34, 136],[51, 34, 136],[51, 34, 136],[51, 34, 136],[51, 34, 136],[51, 34, 136],[51, 34, 136],[51, 34, 136], ...
            [17, 119, 51],[17, 119, 51],[17, 119, 51],[17, 119, 51],[17, 119, 51],[17, 119, 51],[17, 119, 51],[17, 119, 51]}; 
        markers = {'o', 's', 'd', '^', '<', 'p','h','>','o', 's', 'd', '^', '<', 'p','h','>'};
    elseif ai==2
        fieldNames = {'largeoffset','smalloffset'};
        color = {[175, 175, 175],[175, 175, 175],[175, 175, 175],[175, 175, 175],[175, 175, 175],[175, 175, 175],[175, 175, 175],[175, 175, 175], ...
            [0, 0, 0],[0, 0, 0],[0, 0, 0],[0, 0, 0],[0, 0, 0],[0, 0, 0],[0, 0, 0],[0, 0, 0]}; 
        markers = {'o', 's', 'd', '^', '<', 'p','h','>','o', 's', 'd', '^', '<', 'p','h','>'};
    end
    
    % retrive analysis-specific values
    RT = [rtConds(:,2,ai); rtConds(:,1,ai)];
    reboundTime = [probDelay(:,2,ai); probDelay(:,1,ai)];
    peak1 = (latvals{ai}(:,event,1)+eventtimes{ai}(:,event,1)+930);
    peak2 = (latvals{ai}(:,event,2)+eventtimes{ai}(:,event,2)+930);
    peakTime = [peak2; peak1];
    
    % this is without regressing out any effects
    [mainCorr, mainPval] = corr(RT, reboundTime)
    
    % regressing out condition
    [RT_residuals_contCond, reboundTime_residuals_contCond, partialCorrelation_contCond, pVal_contCond] = partial_corr(RT, reboundTime, condition);
    partialCorrelation_contCond
    pVal_contCond
    
    figure(1)
    % original data
    subplot(1,3,1)
    for i=1:length(RT_residuals_contCond)
        scatter(RT(i,1), reboundTime(i,1), 250, markers{i}, 'filled', 'MarkerFaceColor', color{i}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.5)
        hold on
    end
    axis equal
    xlim([500 1000])
    ylim([300 800])
    xticks([500 600 700 800 900 1000])
    yticks([300 400 500 600 700 800])
    xlabel('RT (ms)');
    ylabel('rebound time (ms)');
    title('Original Data');
    grid off;
    [x_fit, y_fit] = fitLine(RT, reboundTime');
    plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', 'k');
    
    set(gca, 'LineWidth', 2);
    set(gca, 'LineWidth', 2);
    
    set(gca, 'FontName', 'Arial', 'FontSize', (911/dpi)*aestheticScalar);
         ax = gca;  % Get the current axis
        ax.XAxis.LineWidth = 1.5;  % Set the X-axis line width
        ax.YAxis.LineWidth = 1.5;  % Set the Y-axis line width
    
    subplot(1,3,2)
    for i=1:length(RT_residuals_contCond)
        scatter(RT_residuals_contCond(i,1), reboundTime_residuals_contCond(i,1), 250, markers{i}, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.5)
        hold on
    end
    axis equal
    xlim([-300 300])
    ylim([-300 300])
    xticks([-150 0 150])
    yticks([-300 -150 0 150 300])
    xlabel('RT residuals (ms)');
    ylabel('rebound time residuals (ms)');
    title('Controlled for Difficulty');
    grid off;
    [x_fit, y_fit] = fitLine(RT_residuals_contCond, reboundTime_residuals_contCond');
    plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', 'k');
    
    set(gca, 'LineWidth', 2);
    set(gca, 'LineWidth', 2);
    
    set(gca, 'FontName', 'Arial', 'FontSize', (911/dpi)*aestheticScalar);
         ax = gca;  % Get the current axis
        ax.XAxis.LineWidth = 1.5;  % Set the X-axis line width
        ax.YAxis.LineWidth = 1.5;  % Set the Y-axis line width
    
    % regressing out subject
    [RT_residuals_contObs, reboundTime_residuals_contObs, partialCorrelation_contObs, pVal_contObs] = partial_corr(RT, reboundTime, observerID);
    partialCorrelation_contObs
    pVal_contObs
    
    subplot(1,3,3)
    for i=1:length(RT_residuals_contObs)
        scatter(RT_residuals_contObs(i,1), reboundTime_residuals_contObs(i,1), 250, 'o', 'filled', 'MarkerFaceColor', color{i}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.5)
        hold on
    end
    axis equal
    xlim([-100 100])
    ylim([-100 100])
    xticks([-50 0 50])
    yticks([-100 -50 0 50 100])
    xlabel('RT residuals (ms)');
    ylabel('rebound time residuals (ms)');
    title('Controlled for Observer');
    grid off;
    sgtitle('RT vs. MS REBOUND')
    [x_fit, y_fit] = fitLine(RT_residuals_contObs, reboundTime_residuals_contObs');
    plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', 'k');
    
    set(gca, 'LineWidth', 2);
    set(gca, 'LineWidth', 2);
    
    set(gca, 'FontName', 'Arial', 'FontSize', (911/dpi)*aestheticScalar);
         ax = gca;  % Get the current axis
        ax.XAxis.LineWidth = 1.5;  % Set the X-axis line width
        ax.YAxis.LineWidth = 1.5;  % Set the Y-axis line width
    
    f1 = gcf;
    f1.Position = [-6 927 990 358]; %[334 778 1483 513]; %[334 778 889 513];
    
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
    if ai==1
        print(gcf, '-dpdf', fullfile(savedir,'fig7a.pdf'));  % for PDF
    elseif ai==2
        print(gcf, '-dpdf', fullfile(savedir,'fig7b.pdf'));  % for PDF
    end
    close all


    % Partial correlation (peak - to - MS rebound)
    
    % this is without regressing out any effects
    [mainCorr, mainPval] = corr(peakTime, reboundTime)
    
    % regressing out condition
    [reboundTime_residuals_contCond, peakTime_residuals_contCond, partialCorrelation_contCond, pVal_contCond] = partial_corr(reboundTime, peakTime, condition);
    partialCorrelation_contCond
    pVal_contCond
    
    figure(2)
    subplot(1,3,1)
    for i=1:length(reboundTime_residuals_contCond)
        scatter(reboundTime(i,1), peakTime(i,1), 250, markers{i}, 'filled', 'MarkerFaceColor', color{i}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.5)
        hold on
    end
    axis equal
    xlim([250 850])
    ylim([1150 1750])
    xticks([250,350,450,550,650,750,850])
    yticks([1150,1270,1390,1510,1630,1750])
    xlabel('rebound time (ms)');
    ylabel('peak dilation (ms)');
    title('Original Data');
    grid off;
    [x_fit, y_fit] = fitLine(reboundTime, peakTime');
    plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', 'k');
    
    set(gca, 'LineWidth', 2);
    set(gca, 'LineWidth', 2);
    
    set(gca, 'FontName', 'Arial', 'FontSize', (911/dpi)*aestheticScalar);
         ax = gca;  % Get the current axis
        ax.XAxis.LineWidth = 1.5;  % Set the X-axis line width
        ax.YAxis.LineWidth = 1.5;  % Set the Y-axis line width
    
    
    subplot(1,3,2)
    for i=1:length(reboundTime_residuals_contCond)
        scatter(reboundTime_residuals_contCond(i,1), peakTime_residuals_contCond(i,1), 250, markers{i}, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.5)
        hold on
    end
    [x_fit, y_fit] = fitLine(reboundTime_residuals_contCond, peakTime_residuals_contCond');
    plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', 'k');
    axis equal
    xlim([-350 350])
    ylim([-350 350])
    xticks([-175 0 175])
    yticks([-350 -175 0 175 350])
    xlabel('rebound time residuals (ms)');
    ylabel('peak dilation residuals (ms)');
    title('Controlling for Condition');
    grid off;
    
    set(gca, 'LineWidth', 2);
    set(gca, 'LineWidth', 2);
    
    set(gca, 'FontName', 'Arial', 'FontSize', (911/dpi)*aestheticScalar);
         ax = gca;  % Get the current axis
        ax.XAxis.LineWidth = 1.5;  % Set the X-axis line width
        ax.YAxis.LineWidth = 1.5;  % Set the Y-axis line width
    
    % regressing out observer
    [reboundTime_residuals_contObs, peakTime_residuals_contObs, partialCorrelation_contObs, pVal_contObs] = partial_corr(reboundTime, peakTime, observerID);
    partialCorrelation_contObs
    pVal_contObs
    
    subplot(1,3,3)
    for i=1:length(reboundTime_residuals_contCond)
        scatter(reboundTime_residuals_contObs(i,1), peakTime_residuals_contObs(i,1), 250, 'o', 'filled', 'MarkerFaceColor', color{i}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.5)
        hold on
    end
    [x_fit, y_fit] = fitLine(reboundTime_residuals_contObs, peakTime_residuals_contObs');
    plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', 'k');
    axis equal
    xlim([-100 100])
    ylim([-100 100])
    xticks([-50 0 50])
    yticks([-100 -50 0 50 100])
    xlabel('rebound time residuals (ms)');
    ylabel('peak dilation residuals (ms)');
    title('Controlling for Observer');
    grid off;
    sgtitle('MS REBOUND vs. PUPIL PEAK')
    
    set(gca, 'LineWidth', 2);
    set(gca, 'LineWidth', 2);
    
    set(gca, 'FontName', 'Arial', 'FontSize', (911/dpi)*aestheticScalar);
         ax = gca;  % Get the current axis
        ax.XAxis.LineWidth = 1.5;  % Set the X-axis line width
        ax.YAxis.LineWidth = 1.5;  % Set the Y-axis line width
    
    f1 = gcf;
    f1.Position = [-6 927 990 358]; %[334 778 1483 513]; %[334 778 889 513];
    
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
    
    if ai==1
        print(gcf, '-dpdf', fullfile(savedir,'fig12a.pdf'));  % for PDF
    elseif ai==2
        print(gcf, '-dpdf', fullfile(savedir,'fig12b.pdf'));  % for PDF
    end


    % just plot the remaining
    
    % this is without regressing out any effects
    [mainCorr, mainPval] = corr(peakTime, RT)
    
    
    figure(3)
    subplot(1,3,1)
    for i=1:length(reboundTime_residuals_contCond)
        scatter(RT(i,1), peakTime(i,1), 250, markers{i}, 'filled', 'MarkerFaceColor', color{i}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.5)
        hold on
    end
    axis equal
    xlim([450 1050])
    ylim([1150 1750])
    xticks([450, 570, 690, 810, 930, 1050])
    yticks([1150, 1270, 1390, 1510, 1630, 1750])
    xlabel('RT (ms)');
    ylabel('peak dilation (ms)');
    title('Original Data');
    grid off;
    [x_fit, y_fit] = fitLine(RT, peakTime');
    plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', 'k');
    
    set(gca, 'LineWidth', 2);
    set(gca, 'LineWidth', 2);
    
    set(gca, 'FontName', 'Arial', 'FontSize', (911/dpi)*aestheticScalar);
         ax = gca;  % Get the current axis
        ax.XAxis.LineWidth = 1.5;  % Set the X-axis line width
        ax.YAxis.LineWidth = 1.5;  % Set the Y-axis line width
    
    subplot(1,3,2)
    for i=1:length(RT_residuals_contCond)
        scatter(RT_residuals_contCond(i,1), peakTime_residuals_contCond(i,1), 250, markers{i}, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.5)
        hold on
    end
    [x_fit, y_fit] = fitLine(RT_residuals_contCond, peakTime_residuals_contCond');
    plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', 'k');
    [correl, pVal] = corr(RT_residuals_contCond, peakTime_residuals_contCond)
    axis equal
    xlim([-300 300])
    ylim([-300 300])
    xticks([-150 0 150])
    yticks([-300 -150 0 150 300])
    xlabel('RT residuals (ms)');
    ylabel('peak dilation residuals (ms)');
    title('Controlling for Condition');
    grid off;
    
    set(gca, 'LineWidth', 2);
    set(gca, 'LineWidth', 2);
    
    set(gca, 'FontName', 'Arial', 'FontSize', (911/dpi)*aestheticScalar);
         ax = gca;  % Get the current axis
        ax.XAxis.LineWidth = 1.5;  % Set the X-axis line width
        ax.YAxis.LineWidth = 1.5;  % Set the Y-axis line width
    
    subplot(1,3,3)
    for i=1:length(RT_residuals_contObs)
        scatter(RT_residuals_contObs(i,1), peakTime_residuals_contObs(i,1), 250, 'o', 'filled', 'MarkerFaceColor', color{i}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.5)
        hold on
    end
    [x_fit, y_fit] = fitLine(RT_residuals_contObs, peakTime_residuals_contObs');
    plot(x_fit, y_fit, '-', 'LineWidth', 2, 'Color', 'k');
    [correl, pVal] = corr(RT_residuals_contObs, peakTime_residuals_contObs)
    axis equal
    xlim([-150 150])
    ylim([-150 150])
    xticks([-75 0 75])
    yticks([-150 -75 0 75 150])
    xlabel('RT residuals (ms)');
    ylabel('peak dilation residuals (ms)');
    title('Controlling for Condition');
    grid off;
    sgtitle('RT vs. PUPIL PEAK')
    
    set(gca, 'LineWidth', 2);
    set(gca, 'LineWidth', 2);
    
    set(gca, 'FontName', 'Arial', 'FontSize', (911/dpi)*aestheticScalar);
         ax = gca;  % Get the current axis
        ax.XAxis.LineWidth = 1.5;  % Set the X-axis line width
        ax.YAxis.LineWidth = 1.5;  % Set the Y-axis line width
    
    f1 = gcf;
    f1.Position = [-6 927 990 358]; %[334 778 1483 513]; %[334 778 889 513];
    
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
    
    if ai==1
        print(gcf, '-dpdf', fullfile(savedir,'fig11a.pdf'));  % for PDF
    elseif ai==2
        print(gcf, '-dpdf', fullfile(savedir,'fig11b.pdf'));  % for PDF
    end

end


%%

function [postsupPeakTime, suppressionTime] = calculateSuppPeak(rate, selection, window_size)
    meanMS = nanmean(rate(selection,:),1);

    window_size = 50; % override with 50 for testing

    smoothed_data = smoothdata(meanMS, 'movmean', window_size);
        
    decreasing_index = find(diff(smoothed_data(1800:end)) < 0, 1);
    postsupPeakTime = decreasing_index+1800;
    
    increasing_index = find(diff(smoothed_data(1300:end)) > 0, 1, 'first');
    suppressionTime = increasing_index+1300;

end

function [var1_resid, var2_resid, partialCorrelation, pVal] = partial_corr(var1, var2, regressVar)

    [~, ~, var1_resid] = regress(var1, [dummyvar(regressVar), ones(size(regressVar))]);
    [~, ~, var2_resid] = regress(var2, [dummyvar(regressVar), ones(size(regressVar))]); 

    % Compute partial correlation between residuals
    [partialCorrelation, pVal] = corr(var1_resid, var2_resid);
end

function [var1_resid, var2_resid, partialCorrelation, pVal] = partial_corr2(var1, var2, regressVar1, regressVar2)

    [~, ~, var1_resid] = regress(var1, [dummyvar(regressVar1), dummyvar(regressVar2)]);
    [~, ~, var2_resid] = regress(var2, [dummyvar(regressVar1), dummyvar(regressVar2)]); 

    % Compute partial correlation between residuals
    [partialCorrelation, pVal] = corr(var1_resid, var2_resid);
end