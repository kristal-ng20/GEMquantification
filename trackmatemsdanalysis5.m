%% Initialize

% Go to path folder in which .xml file is at, add path to msdanalyzer class and Fiji scripts
cd '/Users/kristal/Desktop/MSDanalysis'; %change dir to where the xml files are
folder = pwd;
addpath('/Applications/MATLAB_R2019a.app/bin/tinevez-msdanalyzer-325a510'); %add path to msdanalyser class
addpath('/Applications/Fiji.app/scripts'); %add path to Fiji scripts

% Define file and create holding folder for output
filename = 'GEMold_virus_4d002_Tracks.xml';

[filepath,name,ext] = fileparts(filename);
[parentFolder, deepestFolder] = fileparts(folder);

pname = folder;
dname = name;
mkdir(fullfile(pname,dname))
outfolder = fullfile(pname,dname);

% Opens Trackmate .xml file
clipZ = true; % Remove Z coordinates
scaleT = true; % Use time instead of frame number
[tracks, md] = importTrackMateTracks(filename, clipZ, scaleT);

%% Parameters
% Prompt user input
prompt = {'Enter frames threshold:','Enter time threshold:','Enter fitting percentage:','Enter frame interval (s):'};
dlgtitle = 'Initial Parameters Input';
dims = [1 35];
definput = {'20','20','0.25','0.1'};
initialparametersinput = inputdlg(prompt,dlgtitle,dims,definput)

% Define parameters
framesthresh = str2double(initialparametersinput{1});
timethresh = str2double(initialparametersinput{2});
fit_perc = str2double(initialparametersinput{3}); %for Deff calculation line fitting
dT = str2double(initialparametersinput{4});
%dT = extractfield(md, 'frameInterval');% Time step between acquisition

% Provide simulation parameters
dimensioniality = 2;
SPACE_UNITS = 'µm';
TIME_UNITS = 's';

%% Thresholding

% Tracklengths
tracklengths = cellfun('size',tracks,1);
% figure;
% histogram(tracklengths)
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

% Plot trajectories
figure;
ma.plotTracks;
ma.labelPlotTracks;

% % Plot individual trajectories (use a to select which track to look at)
% a = 2
% figure;
% ma.plotTracks(gca, a);

%% Plot and save individual trajectories (if necessary and if not too many)

tname = 'indiv_traectories';
mkdir(fullfile(outfolder,tname))
indivfolder = fullfile(outfolder,tname);

for a=1:length(ma.tracks);
    h = figure;
    ma.plotTracks(gca, a);
    saveas(h, fullfile(indivfolder,sprintf('Track%d.jpeg',a)));
    delete(h)
end


%% Compute displacement

for j = 1:length(ma.tracks);
    
%     trackdisp(j,1) = {ma.tracks{j}(:,1)}
%     trackdisp(j,2) = {ma.tracks{j}(:,2)}
%     trackdisp(j,3) = {ma.tracks{j}(:,3)}

displacementim(j,1) = {diff(ma.tracks{j})};
end
for k = 1:length(displacementim);
x = displacementim{k}(:,2);
y = displacementim{k}(:,3);
displacementim(k,2) = {sqrt(x.^2 + y.^2)};
end

% 
% slopes = diff(sampley) ./ diff(samplex)
% 
%     for k = 1:(length(ma.tracks{j})-1)
%         displacement(j,1) = mat2cell(ma.tracks{j}(k+1,:) - ma.tracks{j}(k,:));
%         displacement{j}(k,2) = (ma.tracks{j}(k+1,2) - ma.tracks{j}(k,2));
%         displacement{j}(k,3) = (ma.tracks{j}(k+1,3) - ma.tracks{j}(k,3));
%         displacement(k,4) = {
%             d = ((x2 - x1)2 + (y2 - y1)2  + (z2 - z1)2)½
    
   % trackMSD(j,1) = {ma.tracks{j}(1:m,2)}

% ma.tracks 
% 
% t x y 
% Compute displacement (from point before) 

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
dims = [1 35];
definput = {'20'};
velthreshinput = inputdlg(prompt,dlgtitle,dims,definput)

% Define parameters
velthresh = str2double(velthreshinput{1});






%% Determine % of track it changes direction

%use set 2 

samplet = ma.tracks{2}(:,1)
samplex = ma.tracks{2}(:,2)
sampley = ma.tracks{2}(:,3)

% Plotting
dx = gradient(samplex, 0.1);
dy = gradient(sampley, 0.1);
dt = gradient(samplet);

slopes = diff(sampley) ./ diff(samplex)

hdc = pi/2;                             % Critical Heading Change
hdg = atan2(dy,dx);

sec = 2;                                % ?Look Ahead? Time (sec)
dhdg = filter([1 zeros(1,sec*10-2) 1], 2, hdg);
hdgidx = find(abs(dhdg) >= hdc);        % Find ?dhdg? >= pi/2
figure(1)
plot3(samplex, sampley, samplet)
hold on
plot3(samplex(hdgidx), sampley(hdgidx), samplet(hdgidx), '.r')
hold off
grid on


% Calculate number of turns
directions = atan2(diff(sampley), diff(samplex)); % direction of each leg
ddiff = diff(directions);             % angle at each corner 
angle = mod(ddiff+pi, 2*pi) - pi;     % ... in range -pi to pi
E = abs(angle) > 0;                % turns exceeding 0 degrees
nturns = sum(E)

% Based on number of turns, calculate percentage it is "diffusive" - what
% would the cutoff be???
npoints = nnz(samplex)
turnperc = (nturns/npoints)*100


%% Compute MSD

% Compute MSD and plot MSD Curve of all trajectories
ma = ma.computeMSD;
ma.msd
figure
ma.plotMSD;


% Extract MSD values for each track and calculate average MSD for each
% track (takes only first 10 values)
m = 10
trackMSD = {}
for j = 1:length(ma.msd);
    trackMSD(j,1) = {ma.msd{j}(1:m,2)}
end
aveMSDpertrack = cellfun(@mean,trackMSD)


% Extract MSD values for each track and calculate average MSD for each
% track (takes all values)
trackMSD = {}
for j = 1:length(ma.msd);
    trackMSD(j,1) = {ma.msd{j}(:,2)}
end
aveMSDpertrack = cellfun(@mean,trackMSD)


% Plot ensemble mean (weighted mean of all MSD curves)
figure
ma.plotMeanMSD(gca, true)
mmsd = ma.getMeanMSD;
t = mmsd(:,1);
x = mmsd(:,2);
dx = mmsd(:,3) ./ sqrt(mmsd(:,4));
errorbar(t, x, dx, 'k')


%% Diffusion Coefficient
% Fit individual MSD curve, take mean of parameters
ma = ma.fitMSD(fit_perc); %first 30% of curve
good_enough_fit = ma.lfit.r2fit > 0.8; %make changeable??
Dmean = mean(ma.lfit.a(good_enough_fit))/2/ma.n_dim;
Dstd  =  std(ma.lfit.a(good_enough_fit))/2/ma.n_dim;

fprintf('Estimation of the diffusion coefficient from linear fit of the MSD curves:\n')
fprintf('D = %.3g ± %.3g (mean ± std, N = %d)\n', ...
    Dmean, Dstd, sum(good_enough_fit));

% Plot CDF and frequency distributions of Deff values 

figure; 
cdfplot(ma.lfit.a(good_enough_fit)/2/ma.n_dim)
xlabel('Deff (um2.s-1)')
ylabel('CDF(Deff)')
title('CDF Plot')

figure; histogram(ma.lfit.a(good_enough_fit)/2/ma.n_dim,'Normalization','pdf')
histfit(ma.lfit.a(good_enough_fit)/2/ma.n_dim,[],'kernel')


%% Load saved figures
% c=hgload('CDF_etoh_control.fig');
% hold on
% d = cdfplot(ma.lfit.a(good_enough_fit)/ 2 / ma.n_dim)
% d.Color = 'red';
% xlabel('Deff (um2.s-1)')
% ylabel('CDF(Deff)')
% title('CDF Plot')
% legend('EtOH Control','1uM Antimycin','Location','best')


% % Retrieve instantaneous velocities, per track
%  trackV = ma.getVelocities;
% 
%  % Pool track data together
%  TV = vertcat( trackV{:} );
% 
%  % Velocities are returned in a N x (nDim+1) array: [ T Vx Vy ...]. So the
%  % velocity vector in 2D is:
%  V = TV(:, 2:3);
% 
%  % Compute diffusion coefficient
% varV = var(V);
% mVarV = mean(varV); % Take the mean of the two estimates
% Dest = mVarV / 2 * dT;
% 
% fprintf('Estimation from velocities histogram:\n')
% fprintf('D = %.3g %s', ...
%     Dest, [SPACE_UNITS '²/' TIME_UNITS]);
% %end

%%
% % Retrieve instantaneous velocities, per track
%  trackV = ma.getVelocities;
% 
%  % Pool track data together
%  TV = vertcat( trackV{:} );
% 
%  % Velocities are returned in a N x (nDim+1) array: [ T Vx Vy ...]. So the
%  % velocity vector in 2D is:
%  V = TV(:, 2:3);
% 
%  % Compute diffusion coefficient
% varV = var(V);
% mVarV = mean(varV); % Take the mean of the two estimates
% Dest = mVarV / 2 * dT;
% 
% fprintf('Estimation from velocities histogram:\n')
% fprintf('D = %.3g %s, real value was %.3g %s\n', ...
%     Dest, [SPACE_UNITS '²/' TIME_UNITS], D, [SPACE_UNITS '²/' TIME_UNITS]);
