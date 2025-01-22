clc; clear all; %close all;

% replace path to wherever json file lives
jsonparamfile = 'setup.json';
jsonParams = jsondecode(fileread(jsonparamfile));
subjects = jsonParams.subjects.Subjectids; 
protocols = jsonParams.protocols.Protocolids; 
datadir = jsonParams.datadir.Path; 
nSampleCutOff = 4000; 
analysis_type = 'direction'; %'direction' or 'tilt'
savedir = jsonParams.figsavedir.Path; 

dpi = get(0, 'ScreenPixelsPerInch');

 % trialwise will take mean over all trials across sessions and bootstrap
 % with 68% CI
 % if 0, will do SEM across sessions and do permutation test
trialwise = 1; 
percentileRange = 95; % for plotting, do 68% CIs

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
end

direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
time=[0,100];

x = (1 : nSampleCutOff) - 1300;
samplingRateData = 1000; % all 1000 hz after running preprocess_MA_01.m

%%
for ii = 2:4
    figure

    summaryMSPath = fullfile(datadir, subjects{ii}, 'ProcessedData', 'Summary', 'microsaccades');
    summaryFigPath = fullfile(datadir, subjects{ii}, 'ProcessedData', 'Summary', 'figures');

    % for saving figures / variables
    fileName = sprintf('%s_%s_%s_%s_%i_%i', subjects{ii}, analysis_type, fieldNames{1}, fieldNames{2}, trialwise, percentileRange);
    figName = fullfile(summaryFigPath, fileName);
    summaryName = fullfile(summaryMSPath, fileName);

    if isfile(strcat(summaryName, '.mat'))
       load(strcat(summaryName, '.mat'))
       dataReady = 1;
    else
    
        dataReady = 0;
        mkdir(summaryMSPath)
        mkdir(summaryFigPath)
        
        sub_d_car = [];
        sub_d_obl = [];
        RT_marker_car = [];
        RT_marker_obl = [];
        RTcard_all = []; % all RTs including max
        RTobl_all = []; % all RTs including max
        
        rate = struct();
        rateSummary = struct();
        rtSummary = struct();
        rtTrialwise = struct();
    
        for jj = 1 : 2
            if jj== 1 & ii > 6 
                continue
            end
            for kk = 1 : 8
    
                processeddata_folder = fullfile(datadir, subjects{ii}, ...
                    'ProcessedData', protocols{jj});
                MSpath = fullfile(processeddata_folder, 'eyedata','microsaccades');
                MATpath = fullfile(processeddata_folder, 'eyedata','MATs');
    
                dirname = dir(fullfile(MATpath, sprintf('%s*_tab.mat', direction{kk})));
    
                for di=1:length(dirname) % for multiple blocks % 
    
                    directionFilename = dirname(di).name; % for multiple, it automatically sorts alphanumerically
    
                    if di == 1
                        ms_path= fullfile(MSpath,sprintf('%s_microsaccadeMatrix.mat', direction{kk}));
                    elseif di == 2
                        ms_path= fullfile(MSpath,sprintf('%s2_microsaccadeMatrix.mat', direction{kk}));
                    end

                    try % in one case, eye tracker broke for first session, so only session 2 exists (S03)
                        load(ms_path);
                    catch
                        ms_path= fullfile(MSpath,sprintf('%s2_microsaccadeMatrix.mat', direction{kk}));
                        load(ms_path);
                    end
        
                    tab_path = fullfile(MATpath, directionFilename);
                    load(tab_path)
                    
                    % Make this to filter card-oblique, or by location, or by
                    % tilt (or a combination of these)
    
                    for pp=1:length(fieldNames)
                        fieldName = fieldNames{pp};
    
                        if ~isfield(rateSummary, fieldName)
                            rateSummary.(fieldName) = []; rtSummary.(fieldName) = []; rtTrialwise.(fieldName) = []; ...
                                rate.(fieldName) = []; rtUnfiltered.(fieldName) = []; trialData.(fieldName) = [];
                        end
    
                        if strcmp(fieldName, 'cardinal')
                            locations = 0:45:315; directions = 0:90:270; tilts = [0.5, 1, 2, 4, 8];
                        elseif strcmp(fieldName, 'oblique')
                            locations = 0:45:315; directions = 45:90:315; tilts = [0.5, 1, 2, 4, 8];
                        elseif strcmp(fieldName, 'largeoffset')
                            locations = 0:45:315; directions = 0:45:315; tilts = [8];
                        elseif strcmp(fieldName, 'smalloffset')
                            locations = 0:45:315; directions = 0:45:315; tilts = [0.5, 1, 2, 4];
                        elseif strcmp(fieldName, 'horizontalLoc')
                            locations = [0 180]; directions = 0:45:315; tilts = [8];
                        elseif strcmp(fieldName, 'verticalLoc')
                            locations = [90 270]; directions = 0:45:315; tilts = [0.5, 1, 2, 4];
                        elseif strcmp(fieldName, 'easycardinal')
                            locations = 0:45:315; directions = 0:90:270; tilts = [8];
                        elseif strcmp(fieldName, 'hardoblique')    
                            locations = 0:45:315; directions = 45:90:315; tilts = [0.5, 1, 2, 4];
                        end
    
                        [s_mean,RT_marker, RT_trialwise, RT_unfiltered, binaryRate, tabTrials] = computeMSRateRT(MS_TEMP,tab,samplingRateData, locations, directions, tilts, nSampleCutOff);
    
                        if ~all(isnan(s_mean(1000:end))) % checks if the condition is probed in this loop
                            rate.(fieldName) = [rate.(fieldName);binaryRate(:,1000:nSampleCutOff)]; % all trials concatenated
                            rateSummary.(fieldName) = [rateSummary.(fieldName);s_mean(1000:nSampleCutOff)]; % summary per session
                            rtSummary.(fieldName) = [rtSummary.(fieldName), RT_marker]; 
                            rtTrialwise.(fieldName) = [rtTrialwise.(fieldName); RT_trialwise]; % *1000 to convert ms to s 
                            rtUnfiltered.(fieldName) = [rtUnfiltered.(fieldName); RT_unfiltered]; % all RTs (no filter)
                            trialData.(fieldName) = [trialData.(fieldName); tabTrials];
                        end
                    end
                      
                end
            end
        end
    end

    %% subjectwise temporal plot

    for pp=1:length(fieldNames)
        fieldName = fieldNames{pp};

        % Plotting the confidence intervals (optional)
        if trialwise
            bootFile = strcat(summaryName, sprintf('%sbootstraps.mat', fieldName));
            alt_bootFile = strcat(strrep(summaryName, '95', '68'), sprintf('%sbootstraps.mat', fieldName));
            if isfile(bootFile)
                load(bootFile)
            elseif isfile(alt_bootFile)
                load(alt_bootFile)
            else
                disp('Bootstrapping Data .. ')
                bootstrapStatistics = bootstrapData(rate, fieldName, 1); % 1 to convert to rate
                save(strcat(summaryName, sprintf('%sbootstraps.mat', fieldName)), 'bootstrapStatistics');
            end

            gausWindowSize = 100; %50;
            [lowerBound, upperBound] = findCI(bootstrapStatistics, percentileRange);
            meanBoot = mean(bootstrapStatistics,1);
            
            currRate = raster2rate(rate.(fieldName), gausWindowSize);
            upperError = upperBound - meanBoot; %upperBound - currRate;
            lowerError = meanBoot-lowerBound; %currRate-lowerBound;
            
            % replacing with meanBoot instead of currRate
            a = shadedErrorBar(x(1000:end), currRate, [upperError;lowerError], 'lineprops', {'-','Color',color{pp}/255},'transparent',1,'patchSaturation',0.2);
            hold on
            b = plot(x(1000:end), currRate,'-','LineWidth',2, 'Color', color{pp}/255);
            hold on
          
        
        % Plotting SEM across runs
        else
            semVal = std(rateSummary.(fieldName),0,1)/sqrt(size(rateSummary.(fieldName),1));
            a = shadedErrorBar(x(1000:end), mean(rateSummary.(fieldName),1),semVal, 'lineprops', {'-','Color',color{pp}/255},'transparent',1,'patchSaturation',0.2);
            hold on
            b = plot(x(1000:end), mean(rateSummary.(fieldName),1),'-','LineWidth',2, 'Color', color{pp}/255);
            hold on
        end

        xlim([-300,2000])
        ylabel('microsaccade rate (hz)', 'FontSize', (391/dpi)*4)
        xlabel('time (ms)')

        ymax = 2;
        yLimits = ylim; % get current ylim

        if ii==1
            xlim([-300 650])
            ymax = 4.5;
        elseif ii == 2
            xlim([-300 650]); %xlim([-100 850])
        elseif ii == 3
            xlim([-300 650]); %xlim([-100 880])
        elseif ii == 4
            xlim([-300 650]); %xlim([-100 670])
            ymax = 3;
        elseif ii == 5
            xlim([-300 650]); %xlim([-100 785])
        elseif ii == 6
            xlim([-300 650]); %xlim([-100 685])
            ymax = 3;
        elseif ii == 7
            xlim([-300 650]); %xlim([-100 850])
        elseif ii == 8
            xlim([-300 650]); %xlim([-100 850])
            ymax = 3;
        end

    end

     ylim([0,ymax])

     x_patch = [0, 500, 500, 0]; % x-coordinates
     y_patch = [0, 0, ymax, ymax]; % y-coordinates
     fill(x_patch, y_patch, [173, 216, 230]/255, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
     hold on

    save(strcat(summaryName, '.mat'), 'rate', 'rateSummary', 'rtSummary', 'rtTrialwise', 'rtUnfiltered', 'trialData');

    f1 = gcf;
    f1.Position = [5 916 1469 391];
    width_inch = 1469 / dpi;
    height_inch = 391 / dpi;
    
    set(gca, 'LineWidth', 1.5);
    dpi = get(0, 'ScreenPixelsPerInch');
    set(gca, 'FontName', 'Arial', 'FontSize', (391/dpi)*4);
    
    set(gcf, 'PaperOrientation', 'landscape');
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
    set(gcf, 'PaperPosition', [0 0 11.69 (11.69/width_inch)*height_inch]);

    % subtracting because for manuscript showing S02-S04 and relabelled to
    % be sequential starting from 01
    print(gcf, '-dpdf', fullfile(savedir,sprintf('supp_%s_renamedSUB%s.pdf',analysis_type, int2str(ii-1))));  % for PDF

    close all

end


%% Computing confidence / significance

function [increasing_index] = estimateRoughRebound(meanSignal, window_size)
    smoothed_data = smoothdata(nanmean(meanSignal,1), 'movmean', window_size);
    deriv1 = diff(smoothed_data(300:1000));
    X = 200; % stability
   % Find the indices where the derivative is positive
    positive_indices = find(deriv1 > 0);

    % Find the first occurrence where the derivative stays positive for at least X consecutive time points
    increasing_index = -1; % Initialize to -1 as a flag for not found
    for i = 1:length(positive_indices)
        if all(deriv1(positive_indices(i):positive_indices(i) + X - 1) > 0)
            increasing_index = positive_indices(i) + 299; % Adjust for starting index
            break; % Exit loop once found
        end
    end
   
end

