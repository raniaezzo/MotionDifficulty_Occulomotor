clc; clear all;

warning('off', 'MATLAB:interp1:NaNstrip');
warning('off', 'MATLAB:MKDIR:DirectoryExists');
% Convert EDF and Preprocess eyed
% Read JSON file from parent directory for paths and params
jsonparamfile = 'setup.json';
jsonParams = jsondecode(fileread(jsonparamfile));
datadir = jsonParams.datadir.Path; 
scriptsdir = jsonParams.scriptsdir.Path;
edf2asc = jsonParams.edf2ascdir.Path; 
savedir = jsonParams.figsavedir.Path; 

addpath(genpath(scriptsdir))
condition = jsonParams.protocols.Protocolids;
subject = jsonParams.subjects.Subjectids; 
direction = jsonParams.directions.Directionids; 
filter = [];

for ii = 1 : 8
    for jj = 1 : 2 %2 NOT 1:2 (for S07 / S08)
        for kk = 1 : 8 
            
            
             rawdata_allblocks = fullfile(datadir, subject{ii}, ...
                'RawData', condition{jj}); 
            
            
            % check for possibility of multiple blocks
            subdirectories = checkBlockn(rawdata_allblocks);
            numBlocks = length(subdirectories);

            for bi=1:numBlocks

                % for cases with repeated sessions, save both.
                if bi==1
                    blocklabel = '';
                else
                    blocklabel = num2str(bi);
                end

                rawdata_folder = fullfile(rawdata_allblocks, subdirectories{bi});
    
                csv_filepath = fullfile(rawdata_folder,sprintf('expRes%s_1dMotionAsym_Psychophysics_%s.csv', ...
                    subject{ii}, direction{kk}));
                
                if isfile(csv_filepath)
                    M_csv=csvread(csv_filepath);
                    M_csv_fixed = M_csv(M_csv(:,2)>0,:);

                    subjectarray = repmat(ii, [size(M_csv_fixed,1),1,1]);
                    
                    % fix the RT for the early subjects
                    RT_unfiltered = M_csv_fixed(:,13); % all RTs, no filter
                    if any(RT_unfiltered<.5)
                        M_csv_fixed(:,13) = M_csv_fixed(:,13)+.5;
                    end
                    
                    filter = [filter; [subjectarray, M_csv_fixed]];
                elseif ~isfile(csv_filepath) && bi~=2
                    sprintf('Sub%i Protocol%i Dir%i does not exist. Skipping..', ii, jj, kk)
                    sprintf(csv_filepath)
                end
                
            end
        end
    end
end
%%
dpi = get(0, 'ScreenPixelsPerInch');
n = length(subject);
ones_vector = zeros(n, 1);
jittered_vector = ones_vector + (rand(n, 1) - 0.5) * 0.25;

perf_offsetMag = nan(length(subject), 2, 2); % subject, dPrime, RT
perf_dir = nan(length(subject), 2, 2);

cardinalDir_labels = 5:8;
obliqueDir_labels = 1:4;

for ii = 1 : 8
    
    temp = filter(filter(:,1) == ii,:);
    temp_smallOffset = temp(temp(:,7) < 8, :); % small offset
    temp_largeOffset = temp(temp(:,7) >= 8, :); % large offset
    
    [hitRate, falseAlarmRate] = computeHitFA(temp_smallOffset);
    perf_offsetMag(ii,1,1) = calculateDprime(hitRate, falseAlarmRate);
    
    [hitRate, falseAlarmRate] = computeHitFA(temp_largeOffset);
    perf_offsetMag(ii,2,1) = calculateDprime(hitRate, falseAlarmRate);
    
    perf_offsetMag(ii,1,2) = mean(temp_smallOffset(temp_smallOffset(:,14)<3,14)); % RT for small offset
    perf_offsetMag(ii,2,2) = mean(temp_largeOffset(temp_largeOffset(:,14)<3,14)); % RT for large offset
    
    %

    temp_oblDir = temp(ismember(temp(:, 5), obliqueDir_labels), :); % oblique [1,2,3,4]
    temp_cardDir = temp(ismember(temp(:, 5), cardinalDir_labels), :); % cardinal [5,6,7,8]
    
    [hitRate, falseAlarmRate] = computeHitFA(temp_oblDir);
    perf_dir(ii,1,1) = calculateDprime(hitRate, falseAlarmRate);
    
    [hitRate, falseAlarmRate] = computeHitFA(temp_cardDir);
    perf_dir(ii,2,1) = calculateDprime(hitRate, falseAlarmRate);
    
    perf_dir(ii,1,2) = median(temp_oblDir(temp_oblDir(:,14)<5,14)); % RT for oblique
    perf_dir(ii,2,2) = median(temp_cardDir(temp_cardDir(:,14)<5,14)); % RT for cardinal
    
end

iter=1;
for mm=1:2 % dprime or RT
    
    if mm==1 
        metric = "d'";
    elseif mm==2
        metric = 'response time';
    end
    
    for aa=1:2 % tilt or direction
        
        if aa==1
            color = {[17, 119, 51],[51, 34, 136]}; 
            fieldNames = {'cardinal','oblique'};
            performance = perf_dir;
        else
            fieldNames = {'large','small'};
            color = {[0, 0, 0],[175, 175, 175]}; 
            performance = perf_offsetMag;
            
        end
        
        performance = performance(:,:,mm); 
        
        subplot(2,2,iter)
        jittered_vector = jittered_vector(randperm(length(jittered_vector)));
        x1 = ones(8,1)+jittered_vector;
        y1 = performance(:,2);
        hold on
        jittered_vector = jittered_vector(randperm(length(jittered_vector)));
        x2 = 2*ones(8,1)+jittered_vector;
        y2 = performance(:,1);

        for i = 1:length(jittered_vector)
            plot([x1(i), x2(i)], [y1(i), y2(i)], '-', 'Color', [.85 .85 .85], 'LineWidth', 1.5); % Plot line with markers
        end

        plot([1 2], [mean(performance(:,2)), mean(performance(:,1))], 'k', 'LineWidth', 2)
        hold on

        scatter(x1, y1, 200, 'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'MarkerFaceAlpha', 0.5) 
        scatter(x2, y2, 200, 'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2,  'MarkerFaceAlpha', 0.5)

        scatter(ones(1,1), mean(performance(:,2)), 350, 'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', 'w', 'LineWidth', 2); 
        scatter(2*ones(1,1), mean(performance(:,1)), 350, 'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', 'w', 'LineWidth', 2);

        if mm == 1
            sed = std(performance(:,2)-performance(:,1)) / sqrt(8);
            if aa == 1
                % just make it a floating errorbar
                yval = 4; 
            elseif aa == 2
                yval = 4; 
            end
        elseif mm == 2
            sed = std(performance(:,1)-performance(:,2)) / sqrt(8);
            yval = 0.85;
        end
        errorbar(1.5, yval, sed, 'k-', 'LineWidth', 2);
        ax1 = gca;
        ax1.XTick = 1:2;

        if mm==1
            ylabel(sprintf('%s', metric))
        elseif mm==2
            ylabel(sprintf('%s (ms)', metric))
            ylim([0.5 .9])
        end
        
        xlim([0,3])
        
        if mm==1 
            ylim([0, 5])
            ax1.XTickLabel = {};
        elseif mm==2
            ax1.XTickLabel = fieldNames;
            ylim([0.5, 0.9])
            yticks([0.5 0.6 0.7 0.8 0.9])
            yticklabels([0.5 0.6 0.7 0.8 0.9]*1000)
        end
        
        set(gca, 'LineWidth', 1.5);
        set(gca, 'FontName', 'Arial', 'FontSize', (617/dpi)*2.6);
        box off
        
        [h, p, ci, stats] = ttest(performance(:,2), performance(:,1))
        d = computeCohen_d(performance(:,2), performance(:,1), 'paired')
        
        iter = iter+1;
    end
end

f1 = gcf;
f1.Position = [67 720 903 617];

set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
set(gcf, 'PaperPositionMode', 'auto'); 

print(gcf, '-dpdf', fullfile(savedir,'fig4.pdf'));  % for PDF

%%

function [hitRate, falseAlarmRate] = computeHitFA(matrix)

    % Hits
    numClock = sum(matrix(:,12)==1); % n of clockwise trials
    numClock_correct = sum(matrix(matrix(:,12)==1, 15)==1); % num of correct clockwise trials
    hitRate = numClock_correct/numClock;
    
    % False alarms
    numCClock = sum(matrix(:,12)==0); % n of counterclockwise trials
    numClock_incorrect = sum(matrix(matrix(:,12)==1, 15)==0); % num of incorrect clockwise trials
    falseAlarmRate = numClock_incorrect/numCClock;
end


function dPrime = calculateDprime(hit_rate, false_alarm_rate)

    z_hit = norminv(hit_rate);
    z_false_alarm = norminv(false_alarm_rate);
    
    % Calculate d'
    dPrime = z_hit - z_false_alarm;
end