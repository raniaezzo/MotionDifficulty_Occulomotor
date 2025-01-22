% for pupil, retrieve values for current trial and next trial from ASC.
clc; clear all; close all;

jsonparamfile = 'setup.json';
jsonParams = jsondecode(fileread(jsonparamfile));
subjects = jsonParams.subjects.Subjectids; 
protocols = jsonParams.protocols.Protocolids; 
datadir = jsonParams.datadir.Path; 

%%
% PUPIL RESPONSE SHOULD BE DEFINED BY ASC PER TRIAL (trial_start to responseTime + 4500 ms)?

% EVENTS MATRIX (TIME ACROSS TRIAL): 
% (yes) Stimulus on (1300-1800ms), 
% (yes but not in tab file) saccade onset/latency?
% (yes) next trial stimulus on (fixation color change to black)
% (yes; in allevents.mat) [ENABLE: filter trials with ms; saccades; no ms during this period?]
% (yes; in allevents.mat) [ENABLE: filter out trials with SACCADE or EYE CLOSED?
% (yes; in allevents.mat - 0) [ENABLE: filter trials with BLINKS] - just
% interpolate linearlly

% 1 value per trial (tab file)
% Condition: tilt (-8 —> 8)
% Location (PA)
% Motion Direction (1-8)
% Session (1-16)
% response time, 
% response C/IC (feedback), 

direction = {'VU','VL','HL','HR','LL','LR','UL','UR'};
% for blink correction
pad_size = 150;
samplingRate = 1000; % hertz
stimStart = 1300;
stimEnd = 1800;
trialwise=1;
plotOn=0;

% this 0-centers the data based on this window
minWindowStart = 1000; minWindowEnd = 1800;
% this is used to calculate PSC
baselinePSCstart = 1000; baselinePSCend = 1300;

%% model fitting

% for modeling
samplingRate = 1000; % hertz
cutOff = 1300; % this is to make the time meaningful (STIM ON = 0)
event_label = {'pre_stim', 'sti on', 'response'}; % events modeled (other than boxcar)
condlabels = {'trial'}; % arbitrary
baseline = [];   
wnum = 1;

for ii =1:length(subjects)
    subj = subjects{ii};
    savePath = sprintf('%s/%s/ProcessedData/Summary/pupil', datadir, subj);

    if ~isfile(fullfile(savePath,sprintf('%s_allpupilFits.mat', subj)))
        
        load(fullfile(savePath,sprintf('%s_allpupilData.mat', subj)));
    
        [nTrials, ~] = size(allfilteredPupilData);
        
        field_names = {'eventtimes', 'boxtimes', 'samplerate', 'window', 'ampvals', 'boxampvals', 'latvals', 'tmaxval', 'yintval', 'slopeval','numparams','cost','R2','BICrel'};  % Add more field names as needed
        
        % Preallocate the struct array
        for i = 1:numel(field_names)
            [output(1:nTrials).(field_names{i})] = deal([]);  % Initialize each field with empty values
        end
        
        tic
        
        parfor trialNum=1:nTrials
        
            if mod(trialNum, 100) == 0
                sprintf('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Trial# %i ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~', trialNum)
            end
        
            tempData = allfilteredPupilData(trialNum,:); 
            
            timestamps = alltrialTimeStamps(trialNum,:);
            
            TrialEnd = timestamps(4); 
            responseTime = timestamps(3); 
            minDuringStim = timestamps(6); 
            
            event = [800-cutOff, stimStart-cutOff, responseTime-cutOff];
            
            trialwindow = [-cutOff length(tempData)-(cutOff+1)]; 
        
            options = pret_preprocess();
            options.normflag = false;
            options.blinkflag = false;
            
            sj = pret_preprocess({tempData},samplingRate,trialwindow,condlabels,baseline,options);
            
            % create model for the data
            % Pretending that we are naive to the model that we used to create our
            % data, let's create a model to actually fit the data to.
            model = pret_model();
            
            % While the trial window of our task is from -500 to 3500 ms, here we are
            % not interested in what's happening before 0. So
            % let's set the model window to fit only to the region betweeen 0 and 3500
            % ms (the cost function will only be evaluated along this interval).
            model.window = [800-cutOff trialwindow(2)];
            
            % We already know the sampling frequency.
            model.samplerate = samplingRate;
            
            % We also know the event times of our task. Let's also say that we think 
            % there will be a sustained internal signal from precue onset to response 
            % time (0 to 2750 ms).
            model.eventtimes = event;
            model.eventlabels = event_label; %optional
            model.boxtimes = {[stimStart-cutOff stimEnd-cutOff]}; %, [-300 500]}; %round(responseTime-1300)]}; %round(responseTime-1300)]}; % should I make boxcar task or stimulus related?
            model.boxlabels = {'stim'}; %, 'fixation'}; %optional
            
            % Let's say we want to fit a model with the following parameters: 
            % event-related, amplitude, latency, task-related (box) amplitude, 
            % and the tmax of the pupil response function. We turn the other parameters
            % off.
            model.yintflag = false;
            model.slopeflag = false;
            model.tmaxflag = false;
            % Now let's define the bounds for the parameters we decided to fit. We do
            % not have to give values for the y-intercept and slope because we are not
            % fitting them.
            model.ampbounds = repmat([0;200],1,length(model.eventtimes)); % changed to 150 from 100 b/c optimization frequency hits the max
            model.latbounds = [-1000 0 -500; 500 responseTime-stimStart, 1500];
    
            %event_label = {'fix enforced', 'sti on', 'response cue', 'response', 'next trial start'};
            model.boxampbounds = [0;200];
            
            % We need to fill in the values for the y-intercept and slope since we will
            % not be fitting them as parameters.
            model.yintval = 0;
            model.slopeval = 0;
            model.tmaxval = 930; 
            
            % estimate model parameters via pret_estimate_sj
            % Now let's perform the parameter estimation procedure on our subject data.
            % The mean of each condition will be fit independently. For illustration, 
            % let's run only 3 optimizations using one cpu worker (for more 
            % information, see the help files of pret_estimate and pret_estimate_sj).
            options = pret_estimate_sj();
            options.pret_estimate.optimnum = 3;
            % if you want to try fiting the parameters using single trials instead of the mean,
            % use these lines (you'll want to turn off the optimization plots for this):
            %options.trialmode = 'single';
            
            options.pret_estimate.pret_optim.optimplotflag = true; %false;
            
            sj = pret_estimate_sj(sj,model,wnum,options);
        
            output(trialNum) = sj.estim.trial;
        
        end
        toc
        
        % save it out
        save(fullfile(savePath,sprintf('%s_allpupilFits.mat', subj)), 'output');
        delete(gcp); % shut down parallel pool
    end
end

%% helpers

function padded_array = pad_zeros(array, pad_size)
    padded_array = array; % Make a copy of the original array
    zero_indices = find(array == 0); % Find indices of zeros
    
    % Iterate over zero indices
    for idx = zero_indices
        start_pad = max(1, idx - pad_size); % Calculate start padding index
        end_pad = min(length(array), idx + pad_size); % Calculate end padding index
        
        % Pad with zeros if padding is within array boundaries
        padded_array(start_pad:idx-1) = 0;
        padded_array(idx+1:end_pad) = 0;
    end
end