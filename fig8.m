clc; clear all; close all;
rng(0);

load('uniqueCombinations.mat')
jsonparamfile = 'setup.json';
jsonParams = jsondecode(fileread(jsonparamfile));
subjects = jsonParams.subjects.Subjectids; 
protocols = jsonParams.protocols.Protocolids; 
datadir = jsonParams.datadir.Path; 
savedir = jsonParams.figsavedir.Path; 

totaln = nan(40,length(subjects));
acc_count_ms = nan(40,length(subjects));
acc_count_noms = nan(40,length(subjects));

totaln_c_ms = nan(40,length(subjects));
totaln_cc_ms = nan(40,length(subjects));
totaln_c_noms = nan(40,length(subjects));
totaln_cc_noms = nan(40,length(subjects));
acc_c_count_ms =  nan(40,length(subjects));
acc_c_count_noms =  nan(40,length(subjects));
acc_cc_count_ms =  nan(40,length(subjects));
acc_cc_count_noms =  nan(40,length(subjects));

% check whether a MS occurs within a critical window of 1300-1800
% Define the range of columns to check
startCol = 1300;
endCol = 1800;

% directions in tab: 1=45, 2=225, 3=135, 4=315, 5=90, 6=270, 7=0, 8=180
mappingRefs = [45, 225, 135, 315, 90, 270, 0, 180];

%%

for sub=1:length(subjects)
    
% Define the subject ID as a string
    subjectID = subjects{sub};

    % Define the main directories
    baseDir = [datadir, subjectID, '/ProcessedData/'];
    folders = {'Full_distance_radialtangential', 'Full_distance_non_radialtangential'};

    % Initialize cell arrays to store the concatenated matrices
    concatenatedMicrosaccadeMatrices = {};
    concatenatedEventMatrices = {};

    % Loop through the main folders
    for i = 1:length(folders)
        % Define the eyedata/microsaccades subfolder
        subFolder = [baseDir, folders{i}, '/eyedata/microsaccades/'];
        subFolder2 = [baseDir, folders{i}, '/eyedata/MATs/'];

        % Get the list of all files in the subfolder
        allFiles = dir(subFolder2);

        % Extract unique identifiers (e.g., 'HL') by finding 'XX######_tab.mat' files
        uniqueIdentifiers = {};
        for j = 1:length(allFiles)
            [~, fileName, ~] = fileparts(allFiles(j).name);
            if contains(fileName, '_tab')
                uniqueIdentifiers{end+1} = fileName(1:8); % Extract the first 2 characters as the identifier
            end
        end

        if ~isempty(uniqueIdentifiers)
            % Remove duplicates
            uniqueIdentifiers = unique(uniqueIdentifiers);

            % Loop through each unique identifier
            for k = 1:length(uniqueIdentifiers)
                identifier = uniqueIdentifiers{k};

                % Load the corresponding events file
                eventsFilePattern = fullfile(subFolder, ['*', identifier(1:2), '_events.mat']);
                eventsFile = dir(eventsFilePattern);
                if ~isempty(eventsFile)
                    eventsData = load(fullfile(subFolder, eventsFile(1).name));
                    eventsMatrix = eventsData.('EVENTS'); % Load the matrix


                end

                % Load the microsaccadeMatrix file
                microsaccadeFile = fullfile(subFolder2, [identifier, '_tab.mat']);
                if exist(microsaccadeFile, 'file')
                    microsaccadeData = load(microsaccadeFile);
                    microsaccadeMatrix = microsaccadeData.('tab'); % Load the matrix

                end

                eventsSize = size(eventsMatrix,1);
                tabSize = size(microsaccadeMatrix,1);

                cutSize = min(eventsSize, tabSize);

                %if ~isempty(k) && ~isempty(j)
                    % Store the events matrix
                    concatenatedEventMatrices{end+1} = eventsMatrix(1:cutSize,:);

                    % Store the microsaccade matrix
                    concatenatedMicrosaccadeMatrices{end+1} = microsaccadeMatrix(1:cutSize,:);
                %end

            end
        end
    end

    % Optional: Concatenate all loaded matrices (for example, vertically)
    if ~isempty(concatenatedMicrosaccadeMatrices)
        finalMicrosaccadeMatrix = vertcat(concatenatedMicrosaccadeMatrices{:});
    end

    if ~isempty(concatenatedEventMatrices)
        finalEventMatrix = vertcat(concatenatedEventMatrices{:});
    end

    % figure out which trials are clockwise
    standarddirection = finalMicrosaccadeMatrix(:,10);
    targetdirection = finalMicrosaccadeMatrix(:,12);
    newstandarddirections = mappingRefs(standarddirection)';
    
    % just for 0 (edge case due to experimental code): NOT NEEDED
    % mask = targetdirection<= -352;
    % targetdirection(mask) = targetdirection(mask) + 360;
    mask = (targetdirection >= -8 & targetdirection <= -0.5);
    targetdirection(mask) = targetdirection(mask)*-1;
    
    % actual direction - standard
    ccTrials = (newstandarddirections-targetdirection) > 0;

    %%
    % Check each row for at least one '1' in the specified columns
    hasOne = any(finalEventMatrix(:, startCol:endCol), 2);
    tokenStimulus = finalMicrosaccadeMatrix(:,10:11);

    % Find unique combinations of values across the three columns
    %[uniqueCombinations, ~, idx] = unique(tokenStimulus, 'rows');

    % Initialize matrixSelect to store logical columns for each unique combination
    matrixSelect = false(size(tokenStimulus, 1), size(uniqueCombinations, 1));

    % Initialize an index vector to store the matching indices
    idx = zeros(size(tokenStimulus, 1), 1);

    % Loop through each row of tokenStimulus and find the matching index in uniqueCombinations
    for i = 1:size(tokenStimulus, 1)
        % Find the row in uniqueCombinations that matches the current row of tokenStimulus
        match = ismember(uniqueCombinations, tokenStimulus(i, :), 'rows');

        % Store the index of the matching row (if found)
        if any(match)
            idx(i) = find(match);
        else
            idx(i) = NaN; % If no match is found, you can assign NaN or any other value
        end
    end
    
    % Create logical columns for each unique combination
    for i = 1:size(uniqueCombinations, 1)
        matrixSelect(:, i) = idx == i;
    end

    for si=1:40
       msCurrent = matrixSelect(:,si) & hasOne;
       nomsCurrent = matrixSelect(:,si) & ~hasOne; 

       % Count the number of 1s in each vector
        numOnesMs = sum(msCurrent);
        numOnesNoms = sum(nomsCurrent);

        % Determine the target number of 1s (the smaller of the two counts)
        targetOnes = min(numOnesMs, numOnesNoms);

        % Adjust msCurrent if it has more 1s
        if numOnesMs > targetOnes
            % Find the indices of the 1s in msCurrent
            msIndices = find(msCurrent);

            % Randomly select indices to change from 1 to 0
            indicesToChange = randsample(msIndices, numOnesMs - targetOnes);

            % Change the selected 1s to 0s
            msCurrent(indicesToChange) = 0;
        end

        % Adjust nomsCurrent if it has more 1s
        if numOnesNoms > targetOnes
            % Find the indices of the 1s in nomsCurrent
            nomsIndices = find(nomsCurrent);

            % Randomly select indices to change from 1 to 0
            indicesToChange = randsample(nomsIndices, numOnesNoms - targetOnes);

            % Change the selected 1s to 0s
            nomsCurrent(indicesToChange) = 0;
        end

        totaln(si,sub) = targetOnes;
        acc_count_ms(si,sub) = sum(finalMicrosaccadeMatrix(msCurrent, 14));
        acc_count_noms(si,sub) = sum(finalMicrosaccadeMatrix(nomsCurrent, 14));

        % check for clockwise
        cMSidx = find(~ccTrials & msCurrent);
        cnoMSidx = find(~ccTrials & nomsCurrent);
        acc_c_count_ms(si,sub) = sum(finalMicrosaccadeMatrix(cMSidx, 14));
        acc_c_count_noms(si,sub) = sum(finalMicrosaccadeMatrix(cnoMSidx, 14));
        totaln_c_ms(si,sub) = size((finalMicrosaccadeMatrix(cMSidx, 14)),1);
        totaln_c_noms(si,sub) = size((finalMicrosaccadeMatrix(cnoMSidx, 14)),1);

        % check for counterclockwise
        ccMSidx = find(ccTrials & msCurrent);
        ccnoMSidx = find(ccTrials & nomsCurrent);
        acc_cc_count_ms(si,sub) = sum(finalMicrosaccadeMatrix(ccMSidx, 14));
        acc_cc_count_noms(si,sub) = sum(finalMicrosaccadeMatrix(ccnoMSidx, 14));
        totaln_cc_ms(si,sub) = size((finalMicrosaccadeMatrix(ccMSidx, 14)),1);
        totaln_cc_noms(si,sub) = size((finalMicrosaccadeMatrix(ccnoMSidx, 14)),1);

    end
end


%%

% Define the sets of values for the conditions
dirCardinal = [5, 6, 7, 8];
dirOblique = [1, 2, 3, 4];
tiltLarge = [8];
tiltSmall = [0.5, 1, 2, 4];

% Create a boolean column vector based on the conditions
booleanCardLarge = ismember(uniqueCombinations(:, 1), dirCardinal) & ...
                ismember(uniqueCombinations(:, 2), tiltLarge);
            
booleanCardSmall = ismember(uniqueCombinations(:, 1), dirCardinal) & ...
                ismember(uniqueCombinations(:, 2), tiltSmall);
            
booleanOblLarge = ismember(uniqueCombinations(:, 1), dirOblique) & ...
                ismember(uniqueCombinations(:, 2), tiltLarge);

booleanOblSmall = ismember(uniqueCombinations(:, 1), dirOblique) & ...
                ismember(uniqueCombinations(:, 2), tiltSmall);

            
            %%
            
n = length(subjects);
ones_vector = zeros(n, 1);
jittered_vector = ones_vector + (rand(n, 1) - 0.5) * 0.25;

wMScard_small = compute_dPrime(acc_c_count_ms, acc_cc_count_ms, totaln_c_ms, totaln_cc_ms, booleanCardSmall);
wMScard_large = compute_dPrime(acc_c_count_ms, acc_cc_count_ms, totaln_c_ms, totaln_cc_ms, booleanCardLarge);
wMSobli_small = compute_dPrime(acc_c_count_ms, acc_cc_count_ms, totaln_c_ms, totaln_cc_ms, booleanOblSmall);
wMSobli_large = compute_dPrime(acc_c_count_ms, acc_cc_count_ms, totaln_c_ms, totaln_cc_ms, booleanOblLarge);

woMScard_small = compute_dPrime(acc_c_count_noms, acc_cc_count_noms, totaln_c_noms, totaln_cc_noms, booleanCardSmall);
woMScard_large = compute_dPrime(acc_c_count_noms, acc_cc_count_noms, totaln_c_noms, totaln_cc_noms, booleanCardLarge);
woMSobli_small = compute_dPrime(acc_c_count_noms, acc_cc_count_noms, totaln_c_noms, totaln_cc_noms, booleanOblSmall);
woMSobli_large = compute_dPrime(acc_c_count_noms, acc_cc_count_noms, totaln_c_noms, totaln_cc_noms, booleanOblLarge);


%%
% combine all cardinal, and all oblique
color = {[17, 119, 51],[51, 34, 136]}; 

indvSize = 450;
meanSize = 550;

figure
hold on
x1 = ones(length(subjects),1)+jittered_vector;
y1 = wMScard_large;
x2 = 2*ones(length(subjects),1)+jittered_vector;
y2 = woMScard_large;
for i = 1:length(jittered_vector)
    plot([x1(i), x2(i)], [y1(i), y2(i)], '-', 'Color', [.85 .85 .85], 'LineWidth', 1.5); % Plot line with markers
end
hold on
plot([1 2], [mean(wMScard_large), mean(woMScard_large)], 'k', 'LineWidth', 2)
[nu, val, ci, stats] = ttest(wMScard_large-woMScard_large)
d = computeCohen_d(wMScard_large, woMScard_large, 'paired')
scatter(x1, y1, indvSize, 'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2.3, 'MarkerFaceAlpha', 0.5)
hold on
scatter(1, mean(wMScard_large), meanSize,'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', 'w', 'LineWidth', 2.3)
hold on
scatter(x2, y2, indvSize, 'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2.3, 'MarkerFaceAlpha', 0.5)
hold on
scatter(2, mean(woMScard_large), meanSize,'filled', 'MarkerFaceColor', color{1}/255, 'MarkerEdgeColor', 'w', 'LineWidth', 2.3)
hold on

x1 = 3*ones(length(subjects),1)+jittered_vector;
y1 = wMScard_small;
x2 = 4*ones(length(subjects),1)+jittered_vector;
y2 = woMScard_small;
for i = 1:length(jittered_vector)
    plot([x1(i), x2(i)], [y1(i), y2(i)], '-', 'Color', [.85 .85 .85], 'LineWidth', 1.5); % Plot line with markers
end
plot([3 4], [mean(wMScard_small), mean(woMScard_small)], 'k', 'LineWidth', 2)
[nu, val, ci, stats] = ttest(wMScard_small-woMScard_small)
d = computeCohen_d(wMScard_small, woMScard_small, 'paired')
hold on
scatter(x1, y1, indvSize, 'filled', 'MarkerFaceColor', [158, 198, 175]/256, 'MarkerEdgeColor', 'k', 'LineWidth', 2.3, 'MarkerFaceAlpha', 0.25)
hold on
scatter(3, mean(wMScard_small), meanSize,'filled', 'MarkerFaceColor', [158, 198, 175]/256, 'MarkerEdgeColor', 'w', 'LineWidth', 2.3)
hold on
scatter(x2, y2, indvSize, 'filled', 'MarkerFaceColor', [158, 198, 175]/256, 'MarkerEdgeColor', 'k', 'LineWidth', 2.3, 'MarkerFaceAlpha', 0.25)
hold on
scatter(4, mean(woMScard_small), meanSize,'filled', 'MarkerFaceColor', [158, 198, 175]/256, 'MarkerEdgeColor', 'w', 'LineWidth', 2.3)
hold on

x1 = 5*ones(length(subjects),1)+jittered_vector;
y1 = wMSobli_large;
x2 = 6*ones(length(subjects),1)+jittered_vector;
y2 = woMSobli_large;
for i = 1:length(jittered_vector)
    plot([x1(i), x2(i)], [y1(i), y2(i)], '-', 'Color', [.85 .85 .85], 'LineWidth', 1.5); % Plot line with markers
end
plot([5 6], [mean(wMSobli_large), mean(woMSobli_large)], 'k', 'LineWidth', 2)
[nu, val, ci, stats] = ttest(wMSobli_large-woMSobli_large)
d = computeCohen_d(wMSobli_large, woMSobli_large, 'paired')
hold on
scatter(x1, y1, indvSize, 'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2.3, 'MarkerFaceAlpha', 0.5)
hold on
scatter(5, mean(wMSobli_large), meanSize,'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', 'w', 'LineWidth', 2.3)
hold on
scatter(x2, y2, indvSize, 'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', 'k', 'LineWidth', 2.3, 'MarkerFaceAlpha', 0.5)
hold on
scatter(6, mean(woMSobli_large), meanSize,'filled', 'MarkerFaceColor', color{2}/255, 'MarkerEdgeColor', 'w', 'LineWidth', 2.3)
hold on

x1 = 7*ones(length(subjects),1)+jittered_vector;
y1 = wMSobli_small;
x2 = 8*ones(length(subjects),1)+jittered_vector;
y2 = woMSobli_small;
for i = 1:length(jittered_vector)
    plot([x1(i), x2(i)], [y1(i), y2(i)], '-', 'Color', [.85 .85 .85], 'LineWidth', 1.5); % Plot line with markers
end
plot([7 8], [mean(wMSobli_small), mean(woMSobli_small)], 'k', 'LineWidth', 2)
[nu, val, ci, stats] = ttest(wMSobli_small-woMSobli_small)
d = computeCohen_d(wMSobli_small, woMSobli_small, 'paired')
hold on
scatter(x1, y1, indvSize, 'filled', 'MarkerFaceColor', [143, 143, 198]/256, 'MarkerEdgeColor', 'k', 'LineWidth', 2.3, 'MarkerFaceAlpha', 0.25)
hold on
scatter(7, mean(wMSobli_small), meanSize,'filled', 'MarkerFaceColor', [143, 143, 198]/256, 'MarkerEdgeColor', 'w', 'LineWidth', 2.3)
hold on
scatter(x2, y2, indvSize, 'filled', 'MarkerFaceColor', [143, 143, 198]/256, 'MarkerEdgeColor', 'k', 'LineWidth', 2.3, 'MarkerFaceAlpha', 0.25)
hold on
scatter(8, mean(woMSobli_small), meanSize,'filled', 'MarkerFaceColor', [143, 143, 198]/256, 'MarkerEdgeColor', 'w', 'LineWidth', 2.3)
hold on

xlim([0 9])
ylim([0 6])
ylabel("sensitivity (d')")

yticks([0 1 2 3 4 5 6])
yticklabels({'0', '1', '2', '3', '4', '5', ''})
xticks(1:8)
xticklabels({'w', 'w/o', 'w', 'w/o', 'w', 'w/o', 'w', 'w/o'})

set(gca, 'LineWidth', 2);

f1 = gcf;
f1.Position = [334 778 889 513];
dpi = get(0, 'ScreenPixelsPerInch');
set(gca, 'FontName', 'Arial', 'FontSize', (513/dpi)*4.3);

set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
set(gcf, 'PaperPositionMode', 'auto'); 

set(gca, 'LineWidth', 1.5);
print(gcf, '-dpdf', fullfile(savedir,'fig8.pdf'));  % for PDF

%%
function [dPrimeVal] = compute_dPrime(ncorrectC, ncorrectCC, totalClock, totalCClock, chosenVals)

    % index the values for this condition
    ncorrectCC = sum(ncorrectCC(chosenVals,:),1);
    ncorrectC = sum(ncorrectC(chosenVals,:),1);
    totalClock = sum(totalClock(chosenVals,:),1);
    totalCClock = sum(totalCClock(chosenVals,:),1);

    nSubjs = size(ncorrectC, 2);

    for si=1:nSubjs
        % Calculate hit rate and false alarm rate
        hitRate = ncorrectCC(si) / totalCClock(si);
        falseAlarmRate = (totalClock(si) - ncorrectC(si)) / totalClock(si);

        % to control for infinite values
        if hitRate == 1, hitRate = 1 - (1 / (2 * totalCClock(si))); end
        if hitRate == 0, hitRate = 1 / (2 * totalCClock(si)); end
        if falseAlarmRate == 1, falseAlarmRate = 1 - (1 / (2 * totalClock(si))); end
        if falseAlarmRate == 0, falseAlarmRate = 1 / (2 * totalClock(si)); end

        dPrimeVal(si) = norminv(hitRate) - norminv(falseAlarmRate);
    end
    disp('~~~~~~~~~~~~')
end
