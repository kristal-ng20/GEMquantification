%% Initialize paths and folders

%add path to msdanalyser class
message = sprintf('Select path to msdanalyser class.');
uiwait(warndlg(message));
msdclasspath = uigetdir('/Applications/MATLAB_R2019a.app/bin/','Select folder holding msdanalyzer class');
addpath(msdclasspath);

%add path to Fiji scripts
message = sprintf('Select path to Fiji scripts.');
uiwait(warndlg(message));
fijipath = uigetdir('/Applications/Fiji.app/','Select folder holding Fiji scripts');
addpath(fijipath);

%get xml file
message = sprintf('Select .xml file.');
uiwait(warndlg(message));
[filename,folder] = uigetfile('*.xml','Select .xml file containing tracks');
if isequal(filename,0)
   disp('User selected Cancel');
else
   disp(['File path: ', fullfile(folder,filename)]);
end
cd(folder);

% Create holding folder for output
[filepath,name,ext] = fileparts(filename);
[parentFolder, deepestFolder] = fileparts(folder);

pname = folder;
dname = name;
mkdir(fullfile(pname,dname))
outfolder = fullfile(pname,dname);

%% Read Trackmate file

clipZ = true; % Remove Z coordinates
scaleT = true; % Use time instead of frame number
[tracks, md] = importTrackMateTracks(filename, clipZ, scaleT);

tleng1 = length(tracks);
    message = sprintf('Found %d tracks.', tleng1);
    uiwait(warndlg(message));


%% Parameters
% Prompt user input
prompt = {'Enter frames threshold:','Enter time threshold:','Enter fitting percentage:','Enter frame interval (s):'};
dlgtitle = 'Initial Parameters Input';
dims = 1;
definput = {'20','20','0.25','0.1'};
initialparametersinput = inputdlg(prompt,dlgtitle,dims,definput);

% Define parameters
framesthresh = str2double(initialparametersinput{1});
timethresh = str2double(initialparametersinput{2});
fit_perc = str2double(initialparametersinput{3}); %for Deff calculation line fitting
dT = str2double(initialparametersinput{4});
%dT = extractfield(md, 'frameInterval');% Time step between acquisition

% Provide simulation parameters (should remain unchanged)
dimensioniality = 2;
SPACE_UNITS = 'µm';
TIME_UNITS = 's';

%% Initial thresholding (frames and time)

% Tracklengths
tracklengths = cellfun('size',tracks,1);
mintrack = min(tracklengths); % minimum tracklength
maxtrack = max(tracklengths); % maximum tracklength
num_mintrack = nnz(tracklengths==mintrack); % number of occurrences of min track length
num_maxtrack = nnz(tracklengths==maxtrack); % number of occurrences of max track length


% To threshold track lengths (based on track time)
tracks_thresh = tracks(cellfun('size',tracks,1)>=framesthresh);
subset_tracks = cell(length(tracks_thresh),1);
for i = 1:length(tracks_thresh)
    idx = find((tracks_thresh{i}(:,1))-(tracks_thresh{i}(1,1))<timethresh);
   subset_tracks(i,1)= {tracks_thresh{i}(1:length(idx),:)};
end

tleng2 = length(subset_tracks);
    message = sprintf('Initial thresholding done. Found %d tracks. Instantiating analyzer.', tleng2);
    uiwait(warndlg(message));

%% Instantiate analyser
% Instantiate analyzer with subset data
ma = msdanalyzer(dimensioniality, SPACE_UNITS , TIME_UNITS);
ma = ma.addAll(subset_tracks);
disp(ma)

%% Plot trajectories
trajfig = figure;
ma.plotTracks;
ma.labelPlotTracks;

answer1 = questdlg('Would you like to save plot in output folder?', ...
	'Save tracks plot?', ...
	'Yes','Yes, select another folder','No','No');
switch answer1
    case 'Yes'
        namefig1 = 'initialtrajectoriesplot.fig';
        saveas(trajfig,fullfile(outfolder,namefig1),'fig');
        close;
    case 'Yes, select another folder'
        namefig1 = 'initialtrajectoriesplot.fig';
        selfol = uigetdir('/Users/');
        saveas(trajfig,fullfile(selfol,namefig1),'fig');
        close;
      case 'No'
        disp('Not saving.') 
        close;
end


%% Plot and save individual trajectories 

answer2 = questdlg('Would you like to plot and save individual trajectories?', ...
	'Individual Trajectories', ...
	'Yes','No','No');

% Handle response
switch answer2
    case 'Yes'
        tname = 'indiv_traectories';
        mkdir(fullfile(outfolder,tname));
        indivfolder = fullfile(outfolder,tname);
        for a=1:length(ma.tracks)
            indivfig = figure;
            ma.plotTracks(gca, a);
            saveas(indivfig, fullfile(indivfolder,sprintf('Track%d.jpeg',a)));
            delete(indivfig)
        end
    case 'No'
        disp('Not plotting.')
end

%% Compute displacement
message = sprintf('Computing displacement.');
uiwait(warndlg(message));

displacementim = cell(size(ma.tracks));
for j = 1:length(ma.tracks)
displacementim(j,1) = {diff(ma.tracks{j})};
end

for k = 1:length(displacementim)
x = displacementim{k}(:,2);
y = displacementim{k}(:,3);
displacementim(k,2) = {sqrt(x.^2 + y.^2)};
end

%% XY Displacement velocity

message = sprintf('Calculating xy velocities for thresholding.');
uiwait(warndlg(message));

velocities = cell(length(displacementim),3);
for l = 1:length(displacementim)
    dt = displacementim{l}(:,1);
    ddisp = displacementim{l,2}(:,1);
    velocities(l,1) = {ddisp./dt};
    velocities(l,2) = {mean(ddisp./dt)};
    velocities(l,3) = {std(ddisp./dt)};
end

% Prompt user input
prompt = {'Enter lower velocity threshold:','Enter upper velocity threshold:'};
dlgtitle = 'Velocity Thresholding Input';
dims = 1;
definput = {'0.01','0.6'};
velthreshinput = inputdlg(prompt,dlgtitle,dims,definput);

% Define parameters
lowvelthresh = str2double(velthreshinput{1});
upvelthresh = str2double(velthreshinput{2});
indexvel = zeros(size(ma.tracks));
for vvv = 1:length(velocities)
    idv1 = find((lowvelthresh>=velocities{vvv}(:,1))>=upvelthresh);
    %idv2 = (find((velocities{vvv}(:,1))>=velthresh))+1;
    
    if isempty(idv1)
        indexvel(vvv) = 1;
    else
        indexvel(vvv) = 0;
    end
end

% Potentially have a section here to split datasets into diffusive vs
% directed.


% Subset the data again

matrackscopy = ma.tracks;
matracksremoved = ma.tracks;
indexrem = ~indexvel;
matracksremoved(~indexrem) = {''};

subset_tracks_velthresh = subset_tracks;
subset_tracks_velthresh(~indexvel) = {''};
subset_tracks_velthresh = subset_tracks_velthresh(~cellfun('isempty',subset_tracks_velthresh));

tleng3 = length(subset_tracks_velthresh);
    message = sprintf('Subset-ing data based on threshold. Found %d tracks. Re-instantiating analyzer.', tleng3);
    uiwait(warndlg(message));

% Re-instantiate analyzer
ma = msdanalyzer(dimensioniality, SPACE_UNITS , TIME_UNITS);
ma = ma.addAll(subset_tracks_velthresh);
disp(ma)

%% Re-Plot trajectories
trajfig_thresh = figure;
ma.plotTracks;
ma.labelPlotTracks;

answer3 = questdlg('Re-plotting trajectories. Would you like to save plot in output folder?', ...
	'Save tracks plot?', ...
	'Yes','Yes, select another folder','No','No');
switch answer3
    case 'Yes'
        namefig2 = 'threshtrajectoriesplot.fig';
        saveas(trajfig_thresh,fullfile(outfolder,namefig2),'fig');
        close;
    case 'Yes, select another folder'
        namefig2 = 'threshtrajectoriesplot.fig';
        selfol = uigetdir('/Users/');
        saveas(trajfig_thresh,fullfile(selfol,namefig2),'fig');
        close;
      case 'No'
        disp('Not saving.') 
        close;
end

%% Velocity autocorrelation

ma = ma.computeVCorr;
ma.vcorr;
vautoc = figure;
ma.plotMeanVCorr

answer4 = questdlg('Performed velocity autocorrelation. Would you like to save plot in output folder?', ...
	'Save velocity autocorrelation plot?', ...
	'Yes','Yes, select another folder','No','No');
switch answer4
    case 'Yes'
        namefig3 = 'velautocorr.fig';
        saveas(vautoc,fullfile(outfolder,namefig3),'fig');
        close;
    case 'Yes, select another folder'
        namefig3 = 'velautocorr.fig';
        selfol = uigetdir('/Users/');
        saveas(vautoc,fullfile(selfol,namefig3),'fig');
        close;
      case 'No'
        disp('Not saving.') 
        close;
end

%% Velocity analysis
v = ma.getVelocities;
V = vertcat( v{:} );
vhist = figure;
histogram(V(:, 2:end), 50);
box off;
xlabel([ 'Velocity (' SPACE_UNITS '/' TIME_UNITS ')' ]);
ylabel('#');
vmean = mean(V(:,2:end));
vstd = std(V(:,2:end));

answer5 = questdlg('Performed velocity analysis and plotted histogram. Would you like to save plot in output folder?', ...
	'Save histogram?', ...
	'Yes','Yes, select another folder','No','No');
switch answer5
    case 'Yes'
        namefig4 = 'velocityhistogram.fig';
        saveas(vhist,fullfile(outfolder,namefig4),'fig');
        close;
    case 'Yes, select another folder'
        namefig4 = 'velocityhistogram.fig';
        selfol = uigetdir('/Users/');
        saveas(vhist,fullfile(selfol,namefig4),'fig');
        close;
      case 'No'
        disp('Not saving.') 
        close;
end

%% Determine % of track it changes direction

prompt1 = {'Calculating number of turns during track. Please input critical angle (in pi):'};
dlgtitle1 = 'Critical angle input';
dims = 1;
definput1 = {'0'};
angleinput = inputdlg(prompt1,dlgtitle1,dims,definput1);

% Define parameters
critangle = str2double(angleinput{1});
nturns = zeros(length(ma.tracks),2);
slopes = cell(length(ma.tracks),1);
for di = 1:length(ma.tracks)
    samplet = ma.tracks{di}(:,1);
    samplex = ma.tracks{di}(:,2);
    sampley = ma.tracks{di}(:,3);

    % Calculate number of turns
    directions = atan2(diff(sampley), diff(samplex)); % direction of each leg
    ddiff = diff(directions);             % angle at each corner 
    angle = mod(ddiff+pi, 2*pi) - pi;     % ... in range -pi to pi
    E = abs(angle) > critangle;   % turns exceeding 0 degrees
    npoints = nnz(samplex);
    nturns(di,1) = sum(E);
    nturns(di,2) = (sum(E)/npoints)*100; % Based on number of turns, calculate percentage it is "diffusive" - what would the cutoff be???
    slopes(di,1) = {diff(sampley) ./ diff(samplex)}; 
end


%% Plot individual 3D directories

answer5 = questdlg('Would you like to plot and save individual 3D trajectories?', ...
	'3D Trajectories', ...
	'Yes','No','No');

% Handle response
switch answer5
    case 'Yes'
        tname3d = '3d_traectories';
        mkdir(fullfile(outfolder,tname3d));
        indivfolder3d = fullfile(outfolder,tname3d);
        
        for b=1:length(ma.tracks);
            samplet = ma.tracks{b}(:,1);
            samplex = ma.tracks{b}(:,2);
            sampley = ma.tracks{b}(:,3);
            
            dx = gradient(samplex, 0.1);
            dy = gradient(sampley, 0.1);
            dt = gradient(samplet);

            hdc = pi/2;                             % Critical Heading Change
            hdg = atan2(dy,dx);           
            sec = 2;   % ?Look Ahead? Time (sec)
            dhdg = filter([1 zeros(1,sec*10-2) 1], 2, hdg);
            hdgidx = find(abs(dhdg) >= hdc);        % Find ?dhdg? >= pi/2
            indivfig3d = figure;
            plot3(samplex, sampley, samplet)
            hold on
            plot3(samplex(hdgidx), sampley(hdgidx), samplet(hdgidx), '.r')
            hold off
            grid on
            
            saveas(indivfig3d, fullfile(indivfolder3d,sprintf('Track%d.jpeg',b)));
            delete(indivfig3d)
        end
        
    case 'No'
        disp('Not plotting.')
end

%% Compute MSD

% Compute MSD and plot MSD Curve of all trajectories

message = sprintf('Computing MSD.');
uiwait(warndlg(message));

ma = ma.computeMSD;
ma.msd;
msdfig = figure;
ma.plotMSD

answer5 = questdlg('Would you like to save plot in output folder?', ...
	'Save MSD plot?', ...
	'Yes','Yes, select another folder','No','No');
switch answer5
    case 'Yes'
        namefig4 = 'msdplot.fig';
        saveas(msdfig,fullfile(outfolder,namefig4),'fig');
        close;
    case 'Yes, select another folder'
        namefig4 = 'msdplot.fig';
        selfol = uigetdir('/Users/');
        saveas(msdfig,fullfile(selfol,namefig4),'fig');
        close;
      case 'No'
        disp('Not saving.') 
        close;
end

% Extract MSD values for each track and calculate average MSD for each
% track (takes only first 10 values)
% m = 10
% trackMSD = {}
% for j = 1:length(ma.msd);
%     trackMSD(j,1) = {ma.msd{j}(1:m,2)}
% end
% aveMSDpertrack = cellfun(@mean,trackMSD)

% Extract MSD values for each track and calculate average MSD for each
% track (takes all values)
% trackMSD = {}
% for j = 1:length(ma.msd);
%     trackMSD(j,1) = {ma.msd{j}(:,2)}
% end
% aveMSDpertrack = cellfun(@mean,trackMSD)

% Plot ensemble mean (weighted mean of all MSD curves)

message = sprintf('Computing ensemble mean MSD.');
uiwait(warndlg(message));

ensemmsd = figure;
ma.plotMeanMSD(gca, true)
mmsd = ma.getMeanMSD;
t = mmsd(:,1);
x = mmsd(:,2);
dx = mmsd(:,3) ./ sqrt(mmsd(:,4));
errorbar(t, x, dx, 'k')

answer6 = questdlg('Would you like to save plot in output folder?', ...
	'Save ensemble mean MSD plot?', ...
	'Yes','Yes, select another folder','No','No');
switch answer6
    case 'Yes'
        namefig5 = 'ensemblemeanmsdplot.fig';
        saveas(ensemmsd,fullfile(outfolder,namefig5),'fig');
        close;
    case 'Yes, select another folder'
        namefig5 = 'ensemblemeanmsdplot.fig';
        selfol = uigetdir('/Users/');
        saveas(ensemmsd,fullfile(selfol,namefig5),'fig');
        close;
      case 'No'
        disp('Not saving.') 
        close;
end


%% Diffusion Coefficient
% Fit individual MSD curve, take mean of parameters

message = sprintf('Fitting MSD curve to calculate diffusion coefficient.');
uiwait(warndlg(message));

ma = ma.fitMSD(fit_perc); %first 30% of curve
good_enough_fit = ma.lfit.r2fit > 0.8; %make changeable??
Dmean = mean(ma.lfit.a(good_enough_fit))/2/ma.n_dim;
Dstd  =  std(ma.lfit.a(good_enough_fit))/2/ma.n_dim;

message = sprintf('Estimated diffusion coefficient D = %.3g ± %.3g (N = %d)', Dmean, Dstd, sum(good_enough_fit));
uiwait(warndlg(message));

% fprintf('Estimation of the diffusion coefficient from linear fit of the MSD curves:\n')
% fprintf('D = %.3g ± %.3g (mean ± std, N = %d)\n', ...
%     Dmean, Dstd, sum(good_enough_fit));

% Plot CDF and frequency distributions of Deff values 
message = sprintf('Plotting CDF and frequency distributions');
uiwait(warndlg(message));

plotcdf = figure; 
cdfplot(ma.lfit.a(good_enough_fit)/2/ma.n_dim)
xlabel('Deff (um2.s-1)')
ylabel('CDF(Deff)')
title('CDF Plot')

answer7 = questdlg('Would you like to save plot in output folder?', ...
	'Save CDF plot?', ...
	'Yes','Yes, select another folder','No','No');
switch answer7
    case 'Yes'
        namefig6 = 'CDFplot.fig';
        saveas(plotcdf,fullfile(outfolder,namefig6),'fig');
        close;
    case 'Yes, select another folder'
        namefig6 = 'CDFplot.fig';
        selfol = uigetdir('/Users/');
        saveas(plotcdf,fullfile(selfol,namefig6),'fig');
        close;
      case 'No'
        disp('Not saving.') 
        close;
end

freqdist = figure; 
histogram(ma.lfit.a(good_enough_fit)/2/ma.n_dim,'Normalization','pdf')
histfit(ma.lfit.a(good_enough_fit)/2/ma.n_dim,[],'kernel')

answer8 = questdlg('Would you like to save plot in output folder?', ...
	'Save Frequency distribution plot?', ...
	'Yes','Yes, select another folder','No','No');
switch answer8
    case 'Yes'
        namefig7 = 'Deff_freqdist.fig';
        saveas(freqdist,fullfile(outfolder,namefig7),'fig');
        close;
    case 'Yes, select another folder'
        namefig7 = 'Deff_freqdist.fig';
        selfol = uigetdir('/Users/');
        saveas(freqdist,fullfile(selfol,namefig7),'fig');
        close;
      case 'No'
        disp('Not saving.') 
        close;
end

%% Save data to new excel spreadsheet
message = sprintf('Saving and writing data in new Excel spreadsheet.');
uiwait(warndlg(message));

baseFileName = sprintf('%s.xlsx',name);
fullFileName = fullfile(outfolder, baseFileName);

if exist(sprintf('%s',fullFileName), 'file')==2
  delete(sprintf('%s',fullFileName));
end

% Individual displacements
individualdisplacements = vertcat( displacementim{:,1} );
writematrix(individualdisplacements,fullFileName,'Sheet','Individual_displacements')
% XY displacements
xydisplacements = vertcat( displacementim{:,2} );
writematrix(xydisplacements,fullFileName,'Sheet','XY_displacements')
% XY velocities
xyvelocities = vertcat( velocities{:,1} );
writematrix(xyvelocities,fullFileName,'Sheet','XY_velocities')
% XY velocity statistics
xyvelmean = vertcat( velocities{:,2} );
xyvelstd = vertcat( velocities{:,3} );
xyvelocitystats = horzcat(xyvelmean,xyvelstd);
xyvelocitystatst = array2table(xyvelocitystats,'VariableNames',{'Mean','StdDev'});
writetable(xyvelocitystatst,fullFileName,'Sheet','XY_velocity_stats')
% Velocity autocorrelation
vcor = vertcat( ma.vcorr{:,1} );
writematrix(vcor,fullFileName,'Sheet','Velocity_autocorrelation')
% Thresholded XY velocities
velocityanalysis = array2table(V,'VariableNames',{'Time','X','Y'});
writetable(velocityanalysis,fullFileName,'Sheet','Thresholded_XY_velocities')
% Number of turns
turnnum = array2table(nturns,'VariableNames',{'Numberofturns','percentagedirchange'});
writetable(turnnum,fullFileName,'Sheet','Number_of_turns');
% MSD values
msdvalues = vertcat( ma.msd{:,1} );
writematrix(msdvalues,fullFileName,'Sheet','MSD_values');
% Weighted MSD mean values
weightedmsdmeans = array2table(mmsd,'VariableNames',{'dT','weightedmean','weightedstd','Ndegsoffreedom'});
writetable(weightedmsdmeans,fullFileName,'Sheet','Weighted_mean_MSD');
% Diffusion coefficient
diffcoeffvals = [Dmean, Dstd, sum(good_enough_fit)];
diffcoefftab = array2table(diffcoeffvals,'VariableNames',{'Dmean','Dstd','N'});
writetable(diffcoefftab,fullFileName,'Sheet','Diffusion_coefficient');


