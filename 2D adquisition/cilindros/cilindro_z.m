%% OPTOACOUSTIC SIGNAL ACQUISITION OF A CILINDER IN THE X AXIS
% this code computes the B plane acquisition of a cilinder in the Z axis on
% a 2D grid. It filters the signals with a 4th degree butterworth
% band-pass filter to simulate a real RSOM transducer and then allows you
% to navigate the B planes with the slider fuction created. Finally it
% saves the data in a chosen directory. In order to get a cilinder, a high
% number of sources is needed >100000
clear all;
close all;

%% Variables
tic
l=500; % signal length
dx=10e-6; %m
dy=10e-6; %m
vs=1500; %m/s
p0=1; %au
Rs=10e-6; % radius of the spheres

pdetX=(-2e-3:dx:2e-3); % X axis (m) -2 a 2 mm 
pdetY=(-1e-3:dy:1e-3); % Y axis (m) -1 a 1 mm 
Nx=length(pdetX);
Ny=length(pdetY);
N=2e5; % number of sources

t=linspace(100/l*1e-6,1*1e-6+300/l*1e-6,l); %seg 
dt=abs(t(1)-t(2));
a=-0.75e-3;
b=-1.5e-3;
z= a+(b-a).*rand(1,N); %random positions in the z axis
pfuente=zeros(3,N);% m x, y, z
pfuente(3,:)=z;

atemp=zeros(1,l); 
ptot=zeros(1,l); 
opti1=vs*t;
opti2=Rs-vs*t;
S=zeros(Ny,Nx,l);
t1opti=vs*t;
t2opti=Rs-vs*t;

%% RAW DATA ACQUISITION
for i=1:201
    for j=1:101
        pmax=-sqrt((pdetX(i))^2+(pdetY(j))^2)*sind(60)/sind(30); %NA limit
        for s=1:N 
            if  pfuente(3,s)<=pmax
            rd=sqrt((pdetX(i))^2+(pdetY(j))^2+(pfuente(3,s)^2)); %distance between source-detector      
            pin1=p0/2*(1+t1opti./rd).*heaviside(rd+t1opti).*heaviside(-rd+t2opti); %function of the signal generated
            pinr=p0/2*(1-t1opti./rd).*heaviside(-rd+t1opti).*heaviside(rd+t2opti);
            pout=p0/2*(1-t1opti./rd).*heaviside(rd-t1opti).*heaviside(Rs-rd+vs*t);
            ptot=pin1+pinr+pout;
            end 
            atemp=atemp+ptot;
            ptot=0;
        end
       S(j,i,:)=atemp;
       atemp=0;
    end
    hora=datetime(now,'ConvertFrom','datenum');
    hora=datestr(hora);
    display(['Acquiring data: ' num2str(round(i/201*100,2)) '% completed at ' hora])
end
toc
%% FREQUENCY FILTERING
fs=1/(t(2)-t(1));  %sampling freq

S(:,202:end,:)=flip(S(:,1:200,:),2); %simmetry
S(102:end,:,:)=flip(S(1:100,:,:),1); %simmetry

S_f=zeros(Ny,Nx,l);           
fsmin=220e6;
df=fs/l;
f = (0:df:(fs/2));    
fc=[10E6 120E6]; %cut frequency
wn=fc/(fs/2); 
[coefb1,coefa1] = butter(2,wn,'bandpass');  % 4th degree butterworth band pass filter 

for i=1:Nx
    S_f(:,i,:)=(filter(coefb1,coefa1,squeeze(S(:,i,:))'))';
end
for i=1:Ny
    S_f(i,:,:)=(filter(coefb1,coefa1,squeeze(S(i,:,:))'))';
end

%% DATA VISUALIZATION
% figure(1);set(gcf, 'WindowState', 'maximized'); 
% slider(S,'doble',S_f,Nx,Ny,l,t,dx,dy,'XZ');
% 
% figure(2);set(gcf, 'WindowState', 'maximized');
% slider(S,0,S_f,Nx,Ny,l,t,dx,dy,'YZ');

figure(3)
imagesc(squeeze(S_f(:,201,:))'); colorbar; colormap('gray');
xlabel('Nx:0 mm'); title('Plano YZ');
figure(4)
imagesc(squeeze(S_f(101,:,:))'); colorbar; colormap('gray');
xlabel('Ny:0 mm'); title('Plano XZ');


%% DATA SAVING
saveFolderData = 'C:\Users\ljalg\OneDrive - UAM\UAM\4� Curso\TFG\resultados\cilindro_z\2e5\'; %SAVE FOLDER 
fileName = datestr(now, 'yyyymmddHHMMSS');

fileName2   = [ fileName '_cilZ_2e5.mat'];
save([saveFolderData fileName2], 'S_f');

fileName2   = [ fileName '_YZ_2e5.png'];
saveas(figure(3), [saveFolderData fileName2]);

fileName2   = [ fileName '_XZ_2e5.png'];
saveas(figure(4), [saveFolderData fileName2]);