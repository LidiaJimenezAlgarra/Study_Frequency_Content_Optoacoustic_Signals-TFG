%%  Basic reconstruction of a cillinder 
% this code computes the reconstruction of the raw data of a cillinder acquired 
% It creates a kwave grid and stablishes the same grid as the acquisition. To avoid
% possible errors run the code just after the acquisition or load the corresponding data.
% Later, it computes a MIP in each plane and saves the information
close all

%http://www.k-wave.org/documentation/example_pr_3D_fft_planar_sensor.php 
%% LOAD DATA 
%load('') %load data
S_temp=S_f;

%% GRID INFORMATION
l=size(S_temp,3); 
Nx = size(S_temp,2);   % number of grid points in the x (row) direction
Ny = size(S_temp,1);   % number of grid points in the y (column) direction

dx=10e-6; %m
dy=10e-6; %m
dt=abs(t(1)-t(2));
kgrid = kWaveGrid(Nx, dx, Ny, dy); %kwave grid
kgrid.dt=dt;
kgrid.Nt=l;

medium.sound_speed = 1500;  % [m/s]
medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;   % factor to adjust the power law for limiting frequency ranges.
medium.density=1000;           % [kg/m^3]

%% PLANAR RECONSTRUCTION
p_xyz = kspacePlaneRecon(S_temp, kgrid.dx, kgrid.dy, kgrid.dt, ...
medium.sound_speed, 'DataOrder', 'yzt', 'PosCond', true);

%% DATA VISUALIZATION
x_vec=linspace(-2e-3,2e-3,401);
y_vec=linspace(-1e-3,1e-3,201);
%t=linspace(0,1.0667e-06,500);
vs=1500;

sensitivity_field1=max(p_xyz,[],1); %YX
figure (2); 
imagesc(y_vec*1e3,x_vec(100:301)*1e3,squeeze(sensitivity_field1(:,:,100:301))); xlabel('mm'); ylabel('mm');colormap('gray');
axis('square');title('Plano XY')

sensitivity_field2=max(p_xyz,[],3); %XZ
figure (3); 
imagesc(y_vec*1e3,t*1e3*vs,squeeze(sensitivity_field2(:,:,:))); xlabel('mm'); ylabel('mm');colormap('gray');
axis('square');title('Plano YZ')

sensitivity_field3=max(p_xyz,[],2); 
figure (4); 
imagesc(x_vec(100:301),t*1e3*vs,squeeze(sensitivity_field3(:,:,100:301))); xlabel('mm'); ylabel('mm');colormap('gray');
axis('square');title('Plano XZ')


%% DATA SAVING
% saveFolderData = '';
% fileName = datestr(now, 'yyyymmddHHMMSS');
% % 
% fileName2   = [ fileName '_NAME.mat']; 
% save([saveFolderData fileName2], 'p_xyz');
% 
% fileName2   = [ fileName 'cillinder_XY.png'];
% saveas(figure(2), [saveFolderData fileName2]);
% 
% fileName2   = [ fileName 'cillinder_YZ.png'];
% saveas(figure(3), [saveFolderData fileName2]);
% 
% fileName2   = [ fileName 'cillinder_XZ.png'];
% saveas(figure(4), [saveFolderData fileName2]);
