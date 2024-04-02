clear all
close all
%%http://www.k-wave.org/documentation/example_pr_3D_fft_planar_sensor.php#heading3
%% GRID 
Nx = 55;           % number of grid points in the x (row) direction
Ny = 55;           % number of grid points in the y (column) direction
Nz= 51; 
dz=30e-6; 
dx=10e-6;        % grid point spacing in the x direction [m]
dy = 10e-6;      % grid point spacing in the y direction [m]
kgrid= kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
cx=26;
cy=26;
cz=34; 
radius=1;
vs=1500;

%% MEDIUM PROPERTIES
medium.sound_speed = 1500;  % [m/s]
medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;   % factor to adjust the power law for limiting frequency ranges.
medium.density=1000;           % [kg/m^3]

%% SENSOR 
sensor.mask = zeros(kgrid.Nx, kgrid.Ny, kgrid.Nz);
sensor.mask(1:55,1:55,1) = 1;


%% BALL
ball = makeBall(Nx, Ny, Nz, cx, cy, cz, radius);

completo=ball+sensor.mask;
source.p0 = smooth(ball, true);
voxelPlot(completo, 'AxisTight', false, 'Color', [1 0 0], 'Transparency', 0.5); title('Grid kwave')

%% RAW DATA ACQUISITION
    input_args = {'PlotLayout', false, 'PlotPML', false, ...
    'DataCast', 'single', 'CartInterp', 'nearest'};

sd = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});
sd_rs =reshape(sd, 50, 50, kgrid.Nt); %reshape the volume

%% SIGNAL FILTERING
prueba_filt=zeros(201,401,564);           
fs=1/(kgrid.dt); %sampling freq
fsmin=220e6;
df=fs/564;
f = (0:df:(fs/2));    
fc=[10E6 120E6]; %CUT FREQUENCY
wn=fc/(fs/2);  
[coefb1,coefa1] = butter(2,wn,'bandpass'); %band pass filter

for i=1:401
    prueba_filt(:,i,:)=(filter(coefb1,coefa1,squeeze(prueba(:,i,:))'))';
end
for i=1:201
    prueba_filt(i,:,:)=(filter(coefb1,coefa1,squeeze(prueba(i,:,:))'))';
end
for i=1:size(prueba_filt,2) imagesc(squeeze(prueba_filt(:,i,:))');title(i); colorbar;colormap(gray);ylabel('mm');pause(0.05) ;end
plot(kgrid.t_array*vs*1e3,squeeze(prueba_filt(101,201,:)));xlabel('mm');title('Aline central')

%% VOLUME RECONSTRUCTION
 p_xyz = kspacePlaneRecon(prueba_filt, 10e-6, 10e-6, kgrid.dt, ...
   medium.sound_speed, 'DataOrder', 'yzt', 'PosCond', true, 'Plot', true);
imagesc(x_vec(176:225)*1e3,t_array*vs*1e3,squeeze(p_xyz(:,101,:)));colormap(gray);title('Bplane central reconstruido');xlabel('mm')

%% DATA VISUALIZATION
%close all
x_vec=linspace(-2e-3,2e-3,401);
y_vec=linspace(-1e-3,1e-3,201);
t_array=linspace(0,1.126e-06,564);
sensitivity_field1=max(p_xyz,[],1); %YX
figure (2); 
imagesc(y_vec*1e3,x_vec(100:301)*1e3,squeeze(sensitivity_field1(:,:,100:301))); xlabel('mm'); ylabel('mm');colormap('gray');
axis('square');title('Plano YX')

sensitivity_field2=max(p_xyz,[],2); %XZ
figure (3); 
imagesc(x_vec(100:301)*1e3,t_array*1e3*vs,squeeze(sensitivity_field2(:,:,100:301))); xlabel('mm'); ylabel('mm');colormap('gray');
axis('square');title('Plano XZ')

sensitivity_field3=max(p_xyz,[],3); %el máximo de la tercera dimensión 
figure (4); 
imagesc(y_vec*1e3,t_array*vs*1e3,squeeze(sensitivity_field3)); xlabel('mm'); ylabel('mm');colormap('gray');
axis('square');title('Plano YZ')

%% DATA SAVING
% saveFolderData = '';
% fileName = datestr(now, 'yyyymmddHHMMSS');
 
% % fileName2  = [ fileName 'sensor_kwave_sd_rs.mat']; %_kwave_pxyz_modificado.mat
% % save([saveFolderData fileName2], 'sd_rs');
% 
% fileName2  = [ fileName ' esfera_pxyz_kwave.mat']; %_kwave_pxyz_modificado.mat
% % save([saveFolderData fileName2], 'p_xyz');
% 
% fileName2   = [ fileName 'esf_YX.png'];
% saveas(figure(2), [saveFolderData fileName2]);
% 
% fileName2   = [ fileName 'esf_kwave_YZ.png'];
% saveas(figure(3), [saveFolderData fileName2]);
% 
% fileName2   = [ fileName 'esf_kwave_XZ.png'];
% saveas(figure(4), [saveFolderData fileName2]);
