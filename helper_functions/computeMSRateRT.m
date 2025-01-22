function [s_mean, RT_session, RT_trialwise, RT_unfiltered, rate, tabTrials] = computeMSRateRT(MS,tab,samplerate, locations, directions, tilts, nSampleCutOff)

% this is the order of the values
locationdegreesArray = [315,135,225,45,270,90,180,0];
directiondegreesArray = [45,225,135,315,90,270,0,180];

locationIndices = find(ismember(locationdegreesArray, locations));
directionIndices = find(ismember(directiondegreesArray, directions));

ll = size(tab,1);
timepoint = nan(ll,nSampleCutOff);

% initialize the length of each trial as 0s, the rest will be nans
for i = 1 : ll
    di = tab(i,8)-tab(i,2);
    if di > nSampleCutOff 
        di = nSampleCutOff;
    end
    timepoint(i,1:di) = 0;
end

or = 0;
for j = MS(:,10)'
    or = or + 1 ;
    m_start = floor(MS(or,8)*1000/samplerate);
    m_end = floor(MS(or,9)*1000/samplerate);
    timepoint(j,m_start:m_end) = 1 ;
end
    
% Find row indices that meet the criteria
filtered_row_indices = find(ismember(tab(:, 9), locationIndices) & ...
                            ismember(tab(:, 10), directionIndices) & ...
                            ismember(tab(:, 11), tilts));

timepoint = timepoint(filtered_row_indices, :);


% these are calculated after filtering tab:
tabTrials = tab(filtered_row_indices,:);

% s_mean = nanmean(timepoint,1)*60; % not valid anymore

% this computes the rate accounting for the first timepoint of each MS.
% smoothes it using a gaussian window (for this see Palmieri, H., Fernández, A., &
% Carrasco, M. (2023).)
gausWindowSize = 100;

if ~isempty(timepoint)
    s_mean = raster2rate(timepoint, gausWindowSize);
else
    s_mean = nan(1,nSampleCutOff);
end

% NOTE: for some early sessions, the experimental code collected 
% rt without adding stimulus period. To correct this, add 500 to rts from sessions that have
% any value under 500 b/c all values should be > 500
RT_unfiltered = tabTrials(:,13)*1000; % all RTs, no filter
if any(RT_unfiltered<500)
    RT_unfiltered = 500+RT_unfiltered;
end

RT_trialwise = RT_unfiltered;

RT_trialwise(RT_trialwise>=2000) = nan; % this 2 s filter is effectively RTs < 2.5 s (once adding the stim period)

RT_session = nanmean(RT_trialwise);

rate = timepoint;

end