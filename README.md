# MotionDifficulty_Occulomotor

INSTRUCTIONS:

- Change the setup.json file to reflect the correct paths of the data. Make sure this in in your current directory in MATLAB before running.

REQUIRED SOFTWARE:

Matlab (Checked with R2021b, R2022b)

SR Research Eyelink Developers Kit (must be installed on local computer)

PRET Toolbox (included in this Repo)


~~~~~~~~~~~~~~~~~~~~~~~~~ SCRIPTS

fig2b.m: plots main sequence for subject in manuscript

Fig3.m: plots simulated data for pupil events, and an example fit for pupil model.

Fig4.m: computes d' performance and plots d' and RT for easy and difficulty conditions.

Fig5.m: plots group microsaccade rate time series for easy vs. difficulty conditions.

Fig6_7_9_10_11_12: computes and plots MS latency and peak dilation. Plots pupil amplitude. Then runs all correlations and partial correlations -- computes and plots the the MS onset latencies with respect to RT, pupil, etc. Includes ability to change the pupil event being plotted.

Fig8.m: computes d' controlling for stimulus differences, and plots this for with and without micro saccades.


figSupplementary.m: plot the MS rate per condition per subject as specified in the script (either by tilt, direction, location, etc.). Also plots the median timeseries. Shaded error is error of the across session mean.

	computeMSRateRT: creates a rate timeseries from the tab file and the MS file. Does so based on array of direction, location, and tilt values. Also computes overall RT (filtered by 2 s, and unfiltered RT. 

	shadedErrorBar.m: plots the standard error of the mean.


modelPupildata.m: Fits the pupil data (trial wise) and estimates the parameters and saves as allpupilFits.mat.


~~~~~~~~~~~~~~~~~~~~~~~~~ FILES 

Data_DI_wEYE/<SUBJECT>/ProcessedData/<condition>/eyedata/MATs/*_tab.mat

Col 1: trial attempt # (this is 1:N, including trials that are aborted). These numbers are not trialIDs, which would filter out mistrials.
Col 2: trial start based on print TRIAL_START (in eyelink sampling units)
Col 3: Unused
Col 4: trial sync time, usually the same as trial start? (Unused)
Col 5: stimulus start based on STIMULUS_ON (missing for some early subjects - but since stimulus was always 1300ms after TRIAL_START, can be retrieved)
Col 6: when EVENT_CLEAR_SCREEN prints
Col 7: Broke fixation based on 1.5 degree criteria in experiment
Col 8: TRIAL_END
Col 9: polar angle location
Col 10: motion direction (cardinal directions 5:8, oblique directions 1:4)
Col 11: tilt magnitude
Col 12: direction with tilt offset
Col 13: RT - 500
Col 14: Correct/Incorrect
Col 15: # of samples outside of the 1.5 distance from fixation
Col 16: # of samples during blink

Data_DI_wEYE/<SUBJECT>/ProcessedData/<condition>/eyedata/MATs/*_Dat_all.mat & *_Dat_stim.mat

Dat files (dat_stim for stimulus period only, dat_all for trial period only):
Col 1: Time sample
Col 2: X
Col 3: Y
Col 4: Pupil (arbitrary units)

Data_DI_wEYE/<SUBJECT>/ProcessedData/<condition>/eyedata/microsaccades/*_events.mat & *_allevents.mat
Row = Trial
Column = Binary value indicating whether Microsaccade or Not

Data_DI_wEYE/<SUBJECT>/ProcessedData/<condition>/eyedata/microsaccades/*_rate.mat & *_allevents.mat
1 x 2500 file with rate

Data_DI_wEYE/<SUBJECT>/ProcessedData/<condition>/eyedata/microsaccades/*_microsaccadeMatrix.mat
(MS_TEMP) includes the Microsaccade data filtered to be within the correct amplitude range
*Also excludes microsaccades that intersect with the end of the nSampleCutOff which is set to 2500 ms after trial onset
Col 1: onset of saccade
Col 2: end of saccade
Col 3: peak velocity of saccade (vpeak)
Col 4: horizontal component     (dx)  differences in horizontal and vertical positions between the onset and offset points of each saccade
Col 5: vertical component       (dy)
Col 6: horizontal amplitude     (dX) differences in extreme horizontal and vertical positions between the onset and offset points of each saccade. Insights into the extent or range of movement
Col 7: vertical amplitude       (dY)
Col 8: microsaccade start (temporal onset) - redundant with Col 1
Col 9: microsaccade end (temporal onset) - redundant with Col 2
Col 10: trial ID number
Col 11: Amplitude (pythagorean distance based on Col 6 and 7)
Col 12: Direction from amplitude (col 6-7)
Col 13: Direction from components (col 4-5)

Data_DI_wEYE/<SUBJECT>/ProcessedData/<condition>/eyedata/microsaccades/*_summary.mat
(contains summary of microsaccade characteristics for a given session)
Rows = [# of MS in condition, mean AMP, max AMP, min AMP, mean duration, max duration, min duration, mean velocity peak, max velocity peak, min velocity peak]
Col 1: loc 315 (obl loc)
Col 2: loc 135 (obl loc)
Col 3: loc 225 (obl loc)
Col 4: loc 45 (obl loc)
Col 5: loc 270 (card loc)
Col 6: loc 90 (card loc)
Col 7: loc 180 (card loc)
Col 8: loc 0 (card loc)
Col 9: all microsaccades in session (sum)
*4 columns within a session should be 0s, since only 4 locations are tested per session.

Data_DI_wEYE/<SUBJECT>/ProcessedData/Summary/pupil/*_allpupilData.mat
allpupilData.mat contains the filtered pupil data, timestamps, tab, and summary for each trial across the experiment. 

Data_DI_wEYE/<SUBJECT>/ProcessedData/Summary/microsaccades/*_allmsData.mat
Same matrix format as allmsData.mat. 
