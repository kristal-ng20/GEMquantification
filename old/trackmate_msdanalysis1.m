%% Initialize

% Go to path folder in which .xml file is at, add path to msdanalyzer class and Fiji scripts
cd '/Users/kristal/Desktop'; %change dir to where the xml files are
addpath('/Applications/MATLAB_R2019a.app/bin/tinevez-msdanalyzer-325a510'); %add path to msdanalyser class
addpath('/Applications/Fiji.app/scripts'); %add path to Fiji scripts

% Define file and number of frames to cutoff
file = 'GEMold_virus_4d002_Tracks.xml';
frames = 20;

% Opens Trackmate .xml file
clipZ = true; % Remove Z coordinates
scaleT = true; % Use time instead of frame number
[tracks, md] = importTrackMateTracks(file, clipZ, scaleT);


%% Thresholding

% Tracklengths
tracklengths = cellfun('size',tracks,1);
% figure;
% histogram(tracklengths)
mintrack = min(tracklengths); % minimum tracklength
maxtrack = max(tracklengths); % maximum tracklength
num_mintrack = nnz(tracklengths==mintrack); % number of occurences of min track length
num_maxtrack = nnz(tracklengths==maxtrack); % number of occurences of max track length

% To threshold track lengths (based on track time)
tracks_thresh = tracks(cellfun('size',tracks,1)>=frames);
subset_tracks = cell(length(tracks_thresh),1);
for i = 1:length(tracks_thresh)
    idx = find((tracks_thresh{i}(:,1))-(tracks_thresh{i}(1,1))<1);
   subset_tracks(i,1)= {tracks_thresh{i}(1:length(idx),:)};
end
%% Compute MSD

% Provide simulation parameters
dimensioniality = 2;
SPACE_UNITS = 'µm';
TIME_UNITS = 's';
dT = extractfield(md, 'frameInterval');% Time step between acquisition
fit_perc = 0.25; %for Deff calculation line fitting

% Instantiate analyzer with subset data
ma = msdanalyzer(dimensioniality, SPACE_UNITS , TIME_UNITS);
ma = ma.addAll(subset_tracks);
disp(ma)

% Plot trajectories
figure;
ma.plotTracks;
ma.labelPlotTracks;

% Compute MSD and plot MSD Curve of all trajectories
ma = ma.computeMSD;
ma.msd
figure
ma.plotMSD;

% % Extract MSD values for each track and calculate average MSD for each
% % track (takes only first 10 values)
% m = 10
% trackMSD = {}
% for j = 1:length(ma.msd);
%     trackMSD(j,1) = {ma.msd{j}(1:m,2)}
% end
% aveMSDpertrack = cellfun(@mean,trackMSD)


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
good_enough_fit = ma.lfit.r2fit > 0.8;
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
