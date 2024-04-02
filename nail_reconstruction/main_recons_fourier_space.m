clear all
close all
%load('') %LOAD NAIL DATA

%% GRID PARAMETERS
c=1520; %vs (m/s)
t_foc=2.637e-6; 
Fs=1e9; %sampling frequency
dt=1/Fs;
Ntime = size(S,2);

kgrid.dy=dy*1e-3; %um
kgrid.dz=dx*1e-3; %um 
kgrid.dt=1/Fs; 
medium.sound_speed=c;

xLim=xLim*1e-3; 
yLim=yLim.*1e-3; 

kgrid.Nz=round((abs(xLim(1)-xLim(2)))/kgrid.dz)+1;
kgrid.Ny=round((abs(yLim(1)-yLim(2)))/kgrid.dy)+1;

zLim(1)=(1/Fs)*size(S,2);
zLim(2)=0;
%temporal adjustment
t_0 = 0:dt:(Ntime-1)*dt; 
t_0 = t_0 + trigDelay*dt; %useg
t_sp = (t_0 - t_foc) * c * 1e6;  %m
ind=find(t_sp>=0);
S=S(:,min(ind):end);
%% BPLANES
kgrid.Nt=size(S,2);
df=Fs/kgrid.Nt;
f = (0:df:(Fs/2));  
sd_rs= reshape(S,kgrid.Nz,kgrid.Ny, kgrid.Nt);

%% LOW FREQUENCY RECONSTRUCTION
Slow=timeFiltS(double(S), [10 30]*1e6, 1/Fs); 
kgrid.Nt=size(Slow,2);
figure(1); imagesc(Slow'); colormap('gray');title('Low')
figure(2); imagesc(S'); colormap('gray');title('Normal')

sensor_data_rs = reshape(Slow, kgrid.Nz, kgrid.Ny, kgrid.Nt);

Rlow = kspacePlaneRecon(sensor_data_rs, kgrid.dz,kgrid.dy, kgrid.dt, medium.sound_speed,...
    'DataOrder', 'yzt');
                                  
%% HIGH FREQUENCY RECONSTRUCTION
Shigh=timeFiltS(S, [30 120]*1e6, 1/Fs); 
figure(3); imagesc(Shigh'); colormap('gray');title('High')

sensor_data_rs = reshape(Shigh, kgrid.Nz, kgrid.Ny, kgrid.Nt);
Rhigh = kspacePlaneRecon(sensor_data_rs, kgrid.dz, kgrid.dy, kgrid.dt, medium.sound_speed,...
    'DataOrder', 'yzt');

%% IMAGE PROCESSING
projlow=squeeze(max(-Rlow,[],3));
projhigh=squeeze(max(-Rhigh,[],3));
[alphaval] = alphacalc(projlow,projhigh);

fusx= imfuse(projlow,alphaval*projhigh,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
mini=0.11;
maxi=0.5;
fus=imadjust(fusx,[mini mini 0; maxi maxi 1],[]); 

figure(3); imagesc(fus);

%% A LINE EXTRACTION
dt=1e-9;
t=linspace(0, 663*dt,663)*c*1e6; %m
par=squeeze(fus(175:192,70,:));
figure(2);plot(t(175:192),double(par(:,1))/max(double(par(:,1))));hold on;
plot(t(175:192),double(par(:,2))/max(double(par(:,2))));title('Parallel');hold off

per=squeeze(fus(92:100,225,:));
figure(1);plot(t(92:100),double(per(:,1))/max(double(per(:,1))));hold on;
plot(t(92:100),double(per(:,2))/max(double(per(:,2))));title('Perpendicular');hold off
