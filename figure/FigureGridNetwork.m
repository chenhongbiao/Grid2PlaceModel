%% Clean
clear;close all;clc;
%% The size of network
ncells = 128*128; 
%% Create 2-by-ncells preferred direction vector (radians)
dirs = [0 pi/2; pi 3*pi/2];
dirs = repmat(dirs,sqrt(ncells)/2,sqrt(ncells)/2);
dirs = reshape(dirs,1,[]);
dirVects = [cos(dirs); sin(dirs)];
dirVects = round(dirVects);
%% Make x a 2-by-ncells vector of the 2D cell positions on the neural sheet
x = (0:(sqrt(ncells)-1))-(sqrt(ncells)-1)/2; 
[X,Y] = meshgrid(x,x);
x = [reshape(X,1,[]); reshape(Y,1,[])];
%% Make direction marker
east  = find( dirVects(1,:) == 1  & dirVects(2,:) == 0 );
west  = find( dirVects(1,:) == -1 & dirVects(2,:) == 0 );
north = find( dirVects(1,:) == 0  & dirVects(2,:) == 1 );
south = find( dirVects(1,:) == 0  & dirVects(2,:) == -1);
%% Plot it
figure;
hold on;
plot(x(1,east), x(2,east), 'k>');
plot(x(1,west), x(2,west), 'b<');
plot(x(1,north),x(2,north),'r^');
plot(x(1,south),x(2,south),'gv');
hold off;
%% End
