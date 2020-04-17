%% Initialize paths and folders

msdclasspath = uigetdir('/Applications/MATLAB_R2019a.app/','Select folder holding msdanalyzer class');
addpath(msdclasspath);

fijipath = uigetdir('/Applications/Fiji.app/','Select folder holding Fiji scripts');
addpath(fijipath);

[filename,folder] = uigetfile('*.xml','Select .xml file containing tracks');
if isequal(file,0)
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

%% Parameters
% Prompt user input
prompt = {'Enter frames threshold:','Enter time threshold:','Enter fitting percentage:','Enter frame interval (s):'};
dlgtitle = 'Initial Parameters Input';
dims = 1;
definput = {'20','20','0.25','0.1'};
initialparametersinput = inputdlg(prompt,dlgtitle,dims,definput)

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


%% Instantiate analyser
% Instantiate analyzer with subset data
ma = msdanalyzer(dimensioniality, SPACE_UNITS , TIME_UNITS);
ma = ma.addAll(subset_tracks);
disp(ma)

%% Plot trajectories
figure;
ma.plotTracks;
ma.labelPlotTracks;

[filename, pathname] = uiputfile('trajectories.fig', 'Save the file as');
    if isnumeric(filename)
       disp('User pushed cancel. Not saving anything')
    else
       save(fullfile(pname,dname))
    end


%% Plot and save individual trajectories 

answer = questdlg('Would you like to plot and save individual trajectories?', ...
	'Individual Trajectories', ...
	'Yes','No','No');
% Handle response
switch answer
    case 'Yes'
        
        tname = 'indiv_traectories';
        mkdir(fullfile(outfolder,tname));
        indivfolder = fullfile(outfolder,tname);

        for a=1:length(ma.tracks);
            h = figure;
            ma.plotTracks(gca, a);
            saveas(h, fullfile(indivfolder,sprintf('Track%d.jpeg',a)));
            delete(h)
        end
        
    case 'No'
        disp('Not plotting.')
end

%% Compute displacement

for j = 1:length(ma.tracks);
displacementim(j,1) = {diff(ma.tracks{j})};
end

for k = 1:length(displacementim);
x = displacementim{k}(:,2);
y = displacementim{k}(:,3);
displacementim(k,2) = {sqrt(x.^2 + y.^2)};
end

%% XY Displacement velocity

velocities = cell(length(displacementim),3);

for l = 1:length(displacementim);
    dt = displacementim{l}(:,1);
    ddisp = displacementim{l,2}(:,1);
    velocities(l,1) = {ddisp./dt};
    velocities(l,2) = {mean(ddisp./dt)};
    velocities(l,3) = {std(ddisp./dt)};
end

% Prompt user input
prompt = {'Enter velocity threshold:'};
dlgtitle = 'Velocity Thresholding Input';
dims = 1;
definput = {'0.4'};
velthreshinput = inputdlg(prompt,dlgtitle,dims,definput)

% Define parameters
velthresh = str2double(velthreshinput{1});


for vvv = 1:length(velocities)
    idv1 = find((velocities{vvv}(:,1))>=velthresh);
    idv2 = (find((velocities{vvv}(:,1))>=velthresh))+1;
    
    
    
    
end



tracks_thresh = tracks(cellfun('size',tracks,1)>=framesthresh);
subset_tracks = cell(length(tracks_thresh),1);

for i = 1:length(tracks_thresh)
    idx = find((tracks_thresh{i}(:,1))-(tracks_thresh{i}(1,1))<timethresh);
   subset_tracks(i,1)= {tracks_thresh{i}(1:length(idx),:)};
end


%% Velocity autocorrelation

ma = ma.computeVCorr;
ma.vcorr;
figure;
ma.plotMeanVCorr

%% Velocity analysis

v = ma.getVelocities;
V = vertcat( v{:} );
figure;
histogram(V(:, 2:end), 50);
box off;
xlabel([ 'Velocity (' SPACE_UNITS '/' TIME_UNITS ')' ]);
ylabel('#');
vmean = mean(V(:,2:end));
vstd = std(V(:,2:end));

