%% Basic reconstruction of a sphere
% this code computes the reconstruction of the raw data of the sphere acquired in the program
% "esfera.m". It creates a kwave grid and stablishes the same grid as the acquisition. To avoid
% possible errors run the code just after the acquisition or load the corresponding data.
% Later, it computes a MIP in each plane and saves the information

close all
%http://www.k-wave.org/documentation/example_pr_3D_fft_planar_sensor.php 
%% LOAD DATA Y MODIFICARLA PARA REDUCIR MATRICES
S_temp=S;
%load('') %load data

%% GRID INFORMATIONl=size(S_temp,3); %longitud de tiempo y señal
Nx = size(S_temp,2);   % number of grid points in the x (row) direction
Ny = size(S_temp,1);   % number of grid points in the y (column) direction

dt=t(2)-t(1);
dx=10e-6; %m
dy=10e-6; %m
vs=1500; %m/s
p0=1; %au

kgrid = kWaveGrid(Nx, dx, Ny, dy);
kgrid.dt=dt;
kgrid.Nt=l;

medium.sound_speed = 1500;  % [m/s]
medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;   % factor to adjust the power law for limiting frequency ranges.
medium.density=1000;        % [kg/m^3]

%% PLANAR RECONSTRUCTION
p_xyz = kspacePlaneRecon(S_temp, kgrid.dy, kgrid.dx, kgrid.dt, ...
medium.sound_speed, 'DataOrder', 'yzt', 'PosCond', true);

%% DATA VISUALIZATION
sensitivity_field1=max(p_xyz,[],1); %YX
figure (2); 
imagesc(kgrid.y_vec*1e6,kgrid.x_vec(100:301)*1e6,squeeze(sensitivity_field1(:,:,100:301))); xlabel('um'); ylabel('um');colormap('gray');
axis('square');title('XY')

sensitivity_field2=max(p_xyz,[],2); %XZ
figure (3); 
imagesc(kgrid.x_vec(100:301)*1e6, kgrid.t_array*1e6*vs,squeeze(sensitivity_field2(:,:,100:301))); xlabel('um'); ylabel('um');colormap('gray');
axis('square');title('XZ')

sensitivity_field3=max(p_xyz,[],3); 
figure (4); 
imagesc(kgrid.y_vec*1e6,kgrid.t_array*1e6*vs,squeeze(sensitivity_field3)); xlabel('um'); ylabel('um');colormap('gray');
axis('square');title('YZ')

%% DATA SAVING
% saveFolderData = '';
% fileName = datestr(now, 'yyyymmddHHMMSS');
% % 
% fileName2   = [ fileName '_NAME.mat']; 
% save([saveFolderData fileName2], 'p_xyz');
% 
% fileName2   = [ fileName 'esfera_XY.png'];
% saveas(figure(2), [saveFolderData fileName2]);
% 
% fileName2   = [ fileName 'esfera_YZ.png'];
% saveas(figure(3), [saveFolderData fileName2]);
% 
% fileName2   = [ fileName 'esfera_XZ.png'];
% saveas(figure(4), [saveFolderData fileName2]);
