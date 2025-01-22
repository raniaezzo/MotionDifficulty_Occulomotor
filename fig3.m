clc; clear all; %close all;

% replace path to wherever json file lives
jsonparamfile = 'setup.json';
jsonParams = jsondecode(fileread(jsonparamfile));
subjects = jsonParams.subjects.Subjectids; 
protocols = jsonParams.protocols.Protocolids; 
datadir = jsonParams.datadir.Path; 
savedir = jsonParams.figsavedir.Path; 
addpath('PRET')

%%

h = pupilrf(0:4000,10.1,930,0);

figure
shiftal = 500; 
amp = 1.2;
temp = conv(h,amp);
hold on
plot([930+shiftal 930+shiftal], [0 amp],':', 'LineWidth', 5, 'Color', [.5 .5 .5])
set(gca, 'FontName', 'Arial', 'FontSize', 12);
ylabel('pupil area (psc)')
xlabel('time (ms)')
xlim([0 4500])
ylim([0, 1.6])
temp = conv(h,amp);
a1 = gca;
a1.XTick = linspace(0, 4500, 4);
plot([shiftal shiftal], [0 amp],'-', 'LineWidth', 5, 'Color', [0 .7 .5])
plot([shiftal:length(temp)], [temp(1:end-shiftal+1)], 'LineWidth', 5, 'Color', [.3 .3 .3])

box off
set(gca, 'XTick', [], 'YTick', [])

f1 = gcf;
f1.Position = [82 874 397 285];

dpi = get(0, 'ScreenPixelsPerInch');
set(gca, 'LineWidth', 1.5);

set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
set(gca, 'FontName', 'Arial', 'FontSize', (285/dpi)*5.6);
set(gcf, 'PaperPositionMode', 'auto'); 
print(gcf, '-dpdf', fullfile(savedir,'fig3A.pdf'));  % for PDF

%%

figure
x_patch = [0, 500, 500, 0]; % x-coordinates
y_patch = [0, 0, max(temp), max(temp)]; % y-coordinates
fill(x_patch, y_patch, [173, 216, 230]/255, 'FaceAlpha', 0.05, 'EdgeColor', 'none');
hold on
x_patch = [-1300, 0, 0, -1300]; % x-coordinates
y_patch = [0, 0, max(temp), max(temp)]; % y-coordinates
fill(x_patch, y_patch, [255, 255, 0]/255, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold on
x_patch = [500, 4500, 4500, 500]; % x-coordinates
y_patch = [0, 0, max(temp), max(temp)]; % y-coordinates
fill(x_patch, y_patch, [255, 0, 0]/255, 'FaceAlpha', 0.05, 'EdgeColor', 'none');
hold on
plot([1500 1500], [0 1],'-', 'LineWidth', 5, 'Color', [0 .7 .5])
hold on
plot([250 250], [0 1],'-', 'LineWidth', 5, 'Color', [0 .7 .5])
hold on
plot([-300 -300], [0 1],'-', 'LineWidth', 5, 'Color', [0 .7 .5])
hold on
plot([0 0], [0 .5],'-', 'LineWidth', 5, 'Color', [0.3 .3 .3])
hold on
plot([500 500], [0 .5],'-', 'LineWidth', 5, 'Color', [0.3 .3 .3])
hold on
plot([-15 515], [.5 .5],'-', 'LineWidth', 5, 'Color', [0.3 .3 .3])

ylabel('event impulses')
xlabel('time (ms)')
xlim([-1300 2500])
xlim([-500 2500])
a1 = gca;
a1.XTick = -1000:500:2500; 
f1 = gcf;
f1.Position = [73 983 795 285];

ylim([0 1.2])
yticks([]);
set(gca, 'LineWidth', 1.5);
box off
set(gca, 'YTick', [])

set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
set(gca, 'FontName', 'Arial', 'FontSize', (285/dpi)*5.6);
set(gcf, 'PaperPositionMode', 'auto'); 
print(gcf, '-dpdf', fullfile(savedir,'fig3B.pdf'));  % for PDF

%% load
path_exampleSubject = fullfile(datadir, 'S01', 'ProcessedData', 'Summary');
load(fullfile(path_exampleSubject, 'S01trial35_examplePupilFit.mat'))
load(fullfile(path_exampleSubject, 'S01_allpupilData.mat'))

figure
sj = pret_estimate_sj(sj,model,wnum,options);
xline(800, ':', 'Color', [1 0 0 0.05], 'LineWidth', 3)
hold on
x_patch = [0, 500, 500, 0]; % x-coordinates
y_patch = [0, 0, 200, 200]; % y-coordinates
fill(x_patch, y_patch, [173, 216, 230]/255, 'FaceAlpha', 0.15, 'EdgeColor', 'none');
ylim([0 200])
hold on
legend({'data', '', '', '', '', '', 'model'});
legend box off
f1 = gcf;
f1.Position = [53 571 795-(795*.21) 285];

% correction of psc
yticks = get(gca, 'YTick');
yticks = linspace(0,20,5) / (mean(alltrialSignalSummary(:,2))/max(alltrialSignalSummary(:,4)));
set(gca, 'YTickLabel', yticks * (mean(alltrialSignalSummary(:,2))/max(alltrialSignalSummary(:,4))));

set(gca, 'LineWidth', 1.5);
box off
xlabel('time (ms)')
ylabel('pupil area (psc)')

set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'inches', 'PaperSize', [11.69, 8.27]);
set(gca, 'FontName', 'Arial', 'FontSize', (285/dpi)*5.6);
set(gcf, 'PaperPositionMode', 'auto'); 
print(gcf, '-dpdf', fullfile(savedir,'fig3C.pdf'));  % for PDF