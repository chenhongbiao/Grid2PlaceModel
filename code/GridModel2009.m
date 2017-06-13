% (2009) Burak and Fiete - Accurate Path Integration in CAN Models of Grid Cells
% (2012) Eric A. Zilli - Models of grid cell spatial firing published 2005¨C2011
% faustmeow based on eric zilli
%% Initialization
clear; close all; clc; % use 'pack' in the server
%% To save the variables for FigureAttractorWeights:
% % run after a simulation:
% Bu09_full_netsyn = reshape((W*s')',sqrt(ncells),[]);
% Bu09_full_act = reshape(s',sqrt(ncells),[]);
% st = reshape(s,sqrt(ncells),[]);
% sbump = zeros(size(st));
% [v i] = max(s);
% [r c] = ind2sub([sqrt(ncells) sqrt(ncells)],i);
% sbump(c+(-8:8),r+(-8:8)) = st(r+(-8:8),c+(-8:8)); % grab a 17x17 window around the peak as one "bump"
% Bu09_bump_netsyn = reshape((W*reshape(sbump',[],1))',sqrt(ncells),[]);
% Bu09_bump_act = reshape(sbump',sqrt(ncells),[]);
% s1 = zeros(size(s));
% s1(i) = 1;
% Bu09_n1_netsyn = reshape((W*s1')',sqrt(ncells),[]);
% Bu09_n1_act = reshape(s1',sqrt(ncells),[]);
% snorth = (W*s')' + A.*(1+alpha*dirVects'*[0 -2]')';
% Bu09_north_netsyn = reshape((W*(snorth.*(snorth>0))')',sqrt(ncells),[]);
% Bu09_north_act = reshape((snorth.*(snorth>0))',sqrt(ncells),[]);
% figure; imagesc(Bu09_bump_netsyn);  title('Bu09_bump_netsyn');
% figure; imagesc(Bu09_bump_act);     title('Bu09_bump_act');
% figure; imagesc(Bu09_n1_netsyn);    title('Bu09_n1_netsyn');
% figure; imagesc(Bu09_n1_act);       title('Bu09_n1_act');
% figure; imagesc(Bu09_full_netsyn);  title('Bu09_full_netsyn');
% figure; imagesc(Bu09_full_act);     title('Bu09_full_act');
% figure; imagesc(Bu09_north_netsyn); title('Bu09_north_netsyn');
% figure; imagesc(Bu09_north_act);    title('Bu09_north_act');
% save Bu09_WeightFigure_vars.mat Bu09_full_netsyn Bu09_full_act Bu09_bump_netsyn Bu09_bump_act Bu09_n1_netsyn Bu09_n1_act Bu09_north_netsyn Bu09_north_act
%% Simulation parameters
livePlot = 40; % each livePlot/2 ms refresh figures of the simulation
usePeriodicNetwork = 0; % if =0, aperiodic network. if =1, periodic network.

useRealTrajectory = 1;  % if =0, just constant velocity. if =1, load trajectory
constantVelocity = 0*[0.1; 0.1]; % m/s

useCurrentW = 0; % if W exists in namespace, use it instead of loading/generating
loadWIfPossible = 1; % note that generating W is very slooow, so good for saving it
saveW = 1; % save W to disk after generating (can be >1 Gb e.g. if wSparseThresh=0)
%% Cell parameters
tau = 10; % grid cell synapse time constant, ms
if ~useRealTrajectory
  alpha = 0.10315; % input gain
else
  alpha = 50; % input gain for alpha*|v|<<1
  %alpha = 1; % for Sargolini 2006 >> max(vels)
end
%% Network/Weight matrix parameters
ncells = 128*128; % total number of cells in network
a = 1; % if >1, uses local excitatory connections
lambda = 13; % approx the periodicity of the pattern
beta = 3/(lambda^2); % width of the excitatory region
gamma = 1.1*beta; % how much inhibitory region wider than excitatory 

spikeThresh = 0.1; % threshold for plotting a cell as having spiked
%% Simulation parameters
dt = 0.5; % time step, ms
simdur = 540e3; % total simulation time, ms (540s = 9min)
stabilizationTime = 100; % no-velocity time for pattern to form, ms
tind = 0; % time step number for indexing
t = 0; % simulation time variable, ms

% dt = 0.5; % time step, ms
% simdur = 600e3; % ms - for Sargolini 2006 track
% stabilizationTime = 100; % no-velocity time for pattern to form, ms
% tind = 0; % time step number for indexing
% st = 0; % simulation time variable, ms
%% Initial conditions
s = rand(1,ncells); % activation of each cell
%% Firing field plot variables
watchCell = ncells/2-sqrt(ncells)/2; % plot which cell spatial activity
nSpatialBins = 60; % circle 180cm, so we set 60 bins, each bin is 3*3cm
minx = -0.90; maxx = 0.90; % m
miny = -0.90; maxy = 0.90; % m
occupancy = zeros(nSpatialBins);
spikes = zeros(nSpatialBins);
spikeCoords = [];

% watchCell = ncells/2-sqrt(ncells)/2; % plot which cell spatial activity
% nSpatialBins = 100;  % 100 cm * 100cm, so we set 100 bins, each bin is 1*1cm
% minx = -0.50; maxx = 0.50; % m
% miny = -0.50; maxy = 0.50; % m
% occupancy = zeros(nSpatialBins);
% spikes = zeros(nSpatialBins);
% spikeCoords = [];
%% Create 2-by-ncells preferred direction vector (radians)
dirs = [0 pi/2; pi 3*pi/2];
dirs = repmat(dirs,sqrt(ncells)/2,sqrt(ncells)/2);
dirs = reshape(dirs,1,[]);
dirVects = [cos(dirs); sin(dirs)];
%% Make x a 2-by-ncells vector of the 2D cell positions on the neural sheet
x = (0:(sqrt(ncells)-1))-(sqrt(ncells)-1)/2;
[X,Y] = meshgrid(x,x);
x = [reshape(X,1,[]); reshape(Y,1,[])];
cellSpacing = Y(2)-Y(1); % sets length of field shift in recurrent connections
ell = 2*cellSpacing; % offset of center of inhibitory output
% ell = 0; % if l = 0, then no offset of inhibitory ring, no moving pattern
cellDists = sqrt(x(1,:).^2 + x(2,:).^2); % distance from (0,0) for A below
%% Weight matrix
wSparseThresh = -1e-6;
if ~(useCurrentW && exist('W','var'))
  if usePeriodicNetwork
    fname = sprintf('data/W_Bu09_torus_n%d_l%1g.mat',ncells,ell);
  else
    fname = sprintf('data/W_Bu09_aperiodic_n%d_l%1g.mat',ncells,ell);
  end
  
  if loadWIfPossible && exist(fname,'file')
    fprintf('Attempting to load pre-generated W...\n')
    load(fname);
    fprintf('+ Loaded pre-generated W. Using a = %2g, lambda = %2g, beta = %2g, gamma = %2g, ell = %d\n',a,lambda,beta,gamma,ell)
  else
    fprintf('Generating new W. This may take a while. Notifications at 10%% intervals...\n')
    % Define inputs weights for each neuron i one at a time
    W = [];
    for i=1:ncells
      if mod(i,round(ncells/10))==0
        fprintf('Generating weight matrix. %d%% done.\n',round(i/ncells*100))
      end
      if usePeriodicNetwork
        clear squaredShiftLengths;
        % follow Guanella et al 2007's approach to the periodic distance revised by zilli
        % produce an untwisted torus (a periodic network) - 3d ring with triangular pattern 
        % but the bump spacing must be consistent with the size of the torus
        shifts = repmat(x(:,i),1,ncells) - x - ell*dirVects;
        squaredShiftLengths(1,:) = shifts(1,:).^2 + shifts(2,:).^2; % +(0,0)
        shifts = repmat(x(:,i),1,ncells) - x - sqrt(ncells)*[ones(1,ncells); zeros(1,ncells)] - ell*dirVects;
        squaredShiftLengths(2,:) = shifts(1,:).^2 + shifts(2,:).^2; % +(-1,0)
        shifts = repmat(x(:,i),1,ncells) - x + sqrt(ncells)*[ones(1,ncells); zeros(1,ncells)] - ell*dirVects;
        squaredShiftLengths(3,:) = shifts(1,:).^2 + shifts(2,:).^2; % +(1,0)
        shifts = repmat(x(:,i),1,ncells) - x - sqrt(ncells)*[zeros(1,ncells); ones(1,ncells)] - ell*dirVects;
        squaredShiftLengths(4,:) = shifts(1,:).^2 + shifts(2,:).^2; % +(0,-1)
        shifts = repmat(x(:,i),1,ncells) - x + sqrt(ncells)*[zeros(1,ncells); ones(1,ncells)] - ell*dirVects;
        squaredShiftLengths(5,:) = shifts(1,:).^2 + shifts(2,:).^2; % (0,1)
        shifts = repmat(x(:,i),1,ncells) - x + sqrt(ncells)*[ones(1,ncells); ones(1,ncells)] - ell*dirVects;
        squaredShiftLengths(6,:) = shifts(1,:).^2 + shifts(2,:).^2; % (1,1)
        shifts = repmat(x(:,i),1,ncells) - x + sqrt(ncells)*[-1*ones(1,ncells); ones(1,ncells)] - ell*dirVects;
        squaredShiftLengths(7,:) = shifts(1,:).^2 + shifts(2,:).^2; % (-1,1)
        shifts = repmat(x(:,i),1,ncells) - x + sqrt(ncells)*[ones(1,ncells); -1*ones(1,ncells)] - ell*dirVects;
        squaredShiftLengths(8,:) = shifts(1,:).^2 + shifts(2,:).^2; % (1,-1)
        shifts = repmat(x(:,i),1,ncells) - x + sqrt(ncells)*[-1*ones(1,ncells); -1*ones(1,ncells)] - ell*dirVects;
        squaredShiftLengths(9,:) = shifts(1,:).^2 + shifts(2,:).^2; % (-1,-1)
        
        % Select respective least distances:
        squaredShiftLengths = min(squaredShiftLengths);
      else
        shifts = repmat(x(:,i),1,ncells) - x - ell*dirVects;
        squaredShiftLengths = shifts(1,:).^2 + shifts(2,:).^2;
      end
      temp = a*exp(-gamma*squaredShiftLengths) - exp(-beta*squaredShiftLengths);
      temp(temp>wSparseThresh) = 0;
      W = [W; sparse(temp)];
    end

    if saveW
      save(fname,'W','a','lambda','beta','gamma','ell','-v7.3');
    end
  end
end
%% Define envelope function
if usePeriodicNetwork
  % Periodic
  A = ones(size(cellDists));
else
  % Aperiodic
  R = sqrt(ncells)/2;   % radius of main network, in cell-position units
  a0 = sqrt(ncells)/32; % envelope fall-off rate
  dr = sqrt(ncells)/2;  % radius of tapered region, in cell-position units
  A = exp(-a0*(((cellDists)-R+dr)/dr).^2);
  nonTaperedInds = find(cellDists < (R-dr));
  A(nonTaperedInds) = 1;
  % figure; imagesc(reshape(A,sqrt(ncells),sqrt(ncells)));
end
%% Make optional figure of sheet of activity
if livePlot
  h = figure('color','w','name','Activity of sheet of cells on brain''s surface');
  drawnow
end
%% Possibly load trajectory from disk
if useRealTrajectory
  load data/HaftingTraj_centimeters_seconds.mat;
  %load data/HaftingTraj2c2_cm_s.mat
  pos(3,:) = pos(3,:)*1e3; % time units from ms to s
  % interpolate down to simulation time step
  pos = [interp1(pos(3,:),pos(1,:),0:dt:pos(3,end));
         interp1(pos(3,:),pos(2,:),0:dt:pos(3,end));
         interp1(pos(3,:),pos(3,:),0:dt:pos(3,end))];
  pos(1:2,:) = pos(1:2,:)/100; % cm to m
  vels = [diff(pos(1,:)); diff(pos(2,:))]/dt; % m/s
end

% if useRealTrajectory
%   load data/11207-21060501+02_t6c1.mat;
%   % still a circle box ... 100cm diameter
%   % make s to ms for our time unit
%   t = t*1e3;
%   % interpolate down to simulation time step
%   pos = [interp1(t(~isnan(x1)),x1(~isnan(x1)),0:dt:t(end));
%          interp1(t(~isnan(y1)),y1(~isnan(y1)),0:dt:t(end));
%          interp1(t(~isnan(y1)),t(~isnan(y1)), 0:dt:t(end))];
%   pos(1:2,:) = pos(1:2,:)/100; % cm to m
%   vels = [diff(pos(1,:)); diff(pos(2,:))]/dt; % m/s
% end
%% Simulation
fprintf('Simulation starting. Press ctrl+c to end...\n')
while t<simdur
  tind = tind+1;
  t = dt*tind;
  
  % Velocity input
  if t<stabilizationTime
    v = [0; 0]; % m/s
  else
    if ~useRealTrajectory
      v = constantVelocity; % m/s
    else
      v = vels(:,tind); % m/s
    end
  end
  % curDir(tind) = atan2(v(2),v(1)); % rad
  % speed(tind) = sqrt(v(1)^2+v(2)^2);%/dt; % m/s

  % Feedforward input
  B = A.*(1+alpha*dirVects'*v)'; 
  
  % Total synaptic driving currents
  sInputs = (W*s')' + B;

  % Synaptic drive only increases if input cells are over threshold 0
  sInputs = sInputs.*(sInputs>0);

  % Synaptic drive decreases with time constant tau
  s = s + dt*(sInputs - s)/tau;

  % Save firing field information (average value of s in each spatial bin)
  if useRealTrajectory
    if s(watchCell)>spikeThresh
      spikeCoords = [spikeCoords; pos(1,tind) pos(2,tind)];
    end
    xindex = round((pos(1,tind)-minx)/(maxx-minx)*nSpatialBins)+1;
    yindex = round((pos(2,tind)-miny)/(maxy-miny)*nSpatialBins)+1;
    occupancy(yindex,xindex) = occupancy(yindex,xindex) + dt;
    spikes(yindex,xindex) = spikes(yindex,xindex) + s(watchCell);
  end

  if livePlot>0 && (livePlot==1 || mod(tind,livePlot)==1)
    if ~useRealTrajectory
      figure(h);
      set(h,'name','Activity of sheet of cells on brain''s surface');
      imagesc(reshape(s,sqrt(ncells),sqrt(ncells)));
      axis square
      set(gca,'ydir','normal')
      title(sprintf('t = %.1f ms',t))
      drawnow
    else
      figure(h);
      subplot(131);
      imagesc(reshape(s,sqrt(ncells),sqrt(ncells))); % s is firing rate
      axis square
      title('Population activity')
      set(gca,'ydir','normal')
      subplot(132);
      imagesc(spikes./occupancy); % s divide time = average firing rate
      axis square
      set(gca,'ydir','normal')
      title({sprintf('t = %.1f ms',t),'Rate map'})
      subplot(133);
      plot(pos(1,1:tind),pos(2,1:tind));
      hold on;
      if ~isempty(spikeCoords)
        plot(spikeCoords(:,1),spikeCoords(:,2),'r.')
      end
      axis square
      title({'Trajectory (blue)','and spikes (red)'})
      drawnow
    end
  end
end
%% End