close all
% run the code after the cross aqsuisition to avoid possible parameter errors.
% creates a kwave grid with the parameters of the grid of the cross
% acquisition and computes the reconstruction with the Kwave function k kspacePlaneRecon
% Later, a simple MIP of the 3 planes is displayed

%% LOAD DATA
close all
%load("") %load the data
S_temp=S_f; 

%% Parameters
l=size(S_temp,3); %l, t length 
Nx = size(S_temp,2);   % number of grid points in the x (row) direction
Ny = size(S_temp,1);   % number of grid points in the y (column) direction
dt=t(2)-t(1);
dx=10e-6; %m
dy=10e-6; %m
p0=1; %au

kgrid = kWaveGrid(Nx, dx, Ny, dy);%kwave grid
kgrid.dt=dt;
kgrid.Nt=l;

medium.sound_speed = 1500;  % [m/s]
medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;   % factor to adjust the power law for limiting frequency ranges.
medium.density=1000;       % [kg/m^3]

%% PLANAR RECONSTRUCTION
tic
p_xyz = kspacePlaneRecon(S_temp, kgrid.dy, kgrid.dx, kgrid.dt, ...
medium.sound_speed, 'DataOrder', 'yzt', 'PosCond', true);

%% RESULT VISUALIZATION MIP
x_vec=linspace(-2e-3,2e-3,401);
y_vec=linspace(-1e-3,1e-3,201);
t=linspace(0,1.0667e-06,500);
vs=1500;
sensitivity_field1=max(p_xyz,[],1); %YX
figure (2); 
imagesc(y_vec*1e3,x_vec(100:301)*1e3,squeeze(sensitivity_field1(:,:,100:301))); xlabel('mm'); ylabel('mm');colormap('gray');
axis('square');title('XY PLANE')

sensitivity_field2=max(p_xyz,[],3); %XZ
figure (3); 
imagesc(y_vec*1e3,t*1e3*vs,squeeze(sensitivity_field2(:,:,:))); xlabel('mm'); ylabel('mm');colormap('gray');
axis('square');title('YZ PLANE')

sensitivity_field3=max(p_xyz,[],2); %el máximo de la tercera dimensión 
figure (4); 
imagesc(x_vec(100:301),t*1e3*vs,squeeze(sensitivity_field3(:,:,100:301))); xlabel('mm'); ylabel('mm');colormap('gray');
axis('square');title('Plano XZ')

%% SAVING OF RAW DATA
% saveFolderData = ''; %SAVE FOLDER 
% fileName = datestr(now, 'yyyymmddHHMMSS');

% fileName2   = [ fileName '_SAVE_NAME.mat'];
% save([saveFolderData fileName2], 'VARIABLE_NAME');

% fileName2   = [ fileName '_SAVE_NAME.png'];
% saveas(figure(3), [saveFolderData fileName2]);
% 
% fileName2   = [ fileName '_SAVE_NAME.png'];
% saveas(figure(4), [saveFolderData fileName2]);