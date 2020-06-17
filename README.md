# README

3 scripts of importance:
GEM_analysis_with_dialog.m - brings up dialog boxes for file inputs (don't need to set in the script)
GEM_analysis_no_dialog.m - without dialog boxes, runs faster, set inputs directly in the script
combine_CDF.m - plot and combine CDF plots of different conditions in the same plot *could be made into a function*

# FIJI TrackMate

Before using these scripts, you need tracks generated from TrackMate:
 
Helpful links:
Trackmate manual: https://imagej.net/_images/8/85/TrackMate-manual.pdf
 
Open image on Fiji 
Stack order XYTCZ, and Split Channels
Crop ROI
	
Open TrackMate 
Plugins -> Tracking -> TrackMate

Use these parameters:
Time interval - set based on time imaged/no. of frames
Select ‘DoG Detector’
Estimated blob diameter: 0.4 micron
		Threshold: ~ 0.4-0.6 (variable)
		Check ‘Use median filter’ and ‘Do sub-pixel localisation’
		-> Preview to check thresholding
Initial Thresholding: Consider all spots, note down number of spots registered
Select a view: Hyperstack Displayer
Set filters on spots: I select Uniform color
Select a tracker: Simple LAP tracker (does not allow track splitting and merging)
Simple LAP tracker (parameters based on Delarue paper)
Linking max distance: 1 micron
Gap-closing max distance: 1.6 micron
Gap-closing max frame gap: 2
Note down number of tracks found
Set filters on tracks: I select Track index
Display options
Button ‘Analysis’ - gives track statistics (displacement, speed, location)
Save ‘Analysis’ file 
Select an action: Select ‘Export tracks to XML file’ to save tracks file


# GEMquantification (using scripts)

MSD Analysis on Matlab

Download msdanalyzer class from here: https://tinevez.github.io/msdanalyzer/

Also helpful: 
https://uk.mathworks.com/matlabcentral/fileexchange/40692-mean-square-displacement-analysis-of-particles-trajectories
https://tinevez.github.io/msdanalyzer/tutorial/MSDTuto_brownian.html

Once downloaded, unzip and put whole folder in matlab bin folder 
For mac it might be: /Applications/MATLAB_R2019a.app/bin/
For windows it might be: C:\Program Files\MATLAB\R2019a\bin\

%% Initialize paths and folders
Sets paths, calls file names and opens track file
Important to know where the msdanalyser folder containing all the source code, and Fiji scripts folder is!
Makes output folder to store results

%% Parameters
Dialog box to ask for parameters (will show “standard values” but change as needed)
dT value will be value extracted from the tracking file -> based on actual acquisition time information stored in image which i think is more accurate
framesthresh = 20; (number of frames to cut off)
timethresh = 4;
fit_perc = 0.25; %for Deff calculation line fitting
dT = extractfield(md, 'frameInterval'); % Time step between acquisition (put value if known)
dimensioniality = 2;
SPACE_UNITS = 'µm';
TIME_UNITS = 's'

** Notes about thresholding: I’ve realised that if I threshold the time too short - i.e. tracks are trimmed too short - fitting at the end to get Deff does not work (ends up with NaN values), I think a safe time to trim the tracks down is when there’s at least 15-20 frames/points for each track. Fitting percentage can be changed as well, but standard curve fitting is set at 25% anyways. 

%% Initial Thresholding
First round of thresholding based on firstly the number of tracks -> discards any track with less than 20 frames (framesthresh)
Then only considers the first 4 seconds of each of these tracks and trims it down to there (timethresh)

%% Instantiate analyser
Add thresholded tracks to msd analyzer class

%% Plot trajectories
Plot tracks and saves in output folder
Also asks if you want to plot individual trajectories and save them

%% Compute displacement
Computes and stores displacement in x, displacement in y and displacement in xy in a table that will be written in output excel file at the end

%% XY Displacement velocity
Computed velocity based on xy displacement values and dT value 
Proceeds to ask for an upper and lower velocity threshold for final thresholding step -> this is to ensure we only take into account “purely” diffusive tracks as opposed to mixed or purely directed movements
This step removes the tracks that do not meet velocity parameter conditions, and re-instantiates the analyser to add the thresholded tracks into msd analyser class

%% Re-Plot trajectories
Replots all trajectories after final thresholding
Do we need individual trajectory plots here too?

%% Velocity autocorrelation
Computes velocity autocorrelation and plots mean velocity autocorrelation plot and saves in output folder
Velocity autocorrelation - time-averaged value defined over a delay domain, what I understand from this is that it is a value to tell if a particle “remembers” past movements
Defined as velocity at t2 - velocity at t1 (“v is the instantaneous velocity vector and the product is the dot product. The average is taken on all t1 and t2 such that dt=t2-t1.”)
For diffusive behaviours, dT should be around 0

%% Velocity analysis
Calculation of instantaneous velocities
Tbh I don’t think this is absolutely necessary, but resulting histogram plot of instantaneous velocities might be helpful to see distribution of velocities - diffusive behaviour should exhibit gaussian like distribution with centre around 0 

%% Determine % of track it changes direction
To count number of times the particle in each track changes direction based on angle change in its trajectory
Asks for input of critical angle to define a change in direction
Critical angle in pi 
Outputs in excel spreadsheet too

%% Plot individual 3D directories
Optional, code will ask if you want to do it
Plots in x, y, and time -> to see at which time the particle changes direction

%% Compute MSD
Computes MSD and plots all MSD curves
Then calculates and plots mean MSD curve
Saves both plots in output folder

%% Diffusion Coefficient
Fit MSD curve based on fit_perc parameter
Shows and saves Deff value in output excel spreadsheet
Plots CDF plot and saves in output folder too

%% Save data to new excel spreadsheet
Saves:
Parameters
Individual displacements
XY displacements
XY velocities
XY velocity statistics
Velocity autocorrelation
Thresholded XY velocities
Number of turns
MSD values
Weighted MSD mean values
Diffusion coefficient



