%% OPTOACOUSTIC SIGNAL ACQUISITION FROM A SPHERE STRUCTURE
% this code computes an acquisition of the optoacoustic signals generated by a individual sphere
% on the center of the grid. The length of the time-resolved signals is l and the spheres are assumed
% of a diameter (Rs) of 10um.
% Later, it filters the generated signals and displays their FFT 

clear all;
close all;

%% GRID 
tic
l=500; % t points
dx=10e-6; %m
dy=10e-6; %m
vs=1500; %m/s
p0=1; %u.au
Rs=10e-6; %m sphere radius

pdetX=(-2e-3:dx:2e-3); %-2 a 2 mm 
pdetY=(-1e-3:dy:1e-3); %-1 a 1 mm 
Nx=length(pdetX);
Ny=length(pdetY);

t=linspace(0,1e-6,l); %seg 
pfuente=[0,0,-1e-3]; %m x, y, z
S=zeros(Ny,Nx,l);

%% RAW DATA ACQUISITION
for i=1:201
    for j=1:101
        if (pdetX(i)<=5e-4 && pdetX(i)>=-5e-4) && (pdetY(j)>=-5e-4 && pdetY(j)<=5e-4)
            rd=sqrt((pdetX(i)-pfuente(1))^2+(pdetY(j)-pfuente(2))^2+(pfuente(3)^2));
            pin1=p0/2*(1+vs*t./rd).*heaviside(rd+vs*t).*heaviside(Rs-rd-vs*t); 
            pinr=p0/2*(1-vs*t./rd).*heaviside(-rd+vs*t).*heaviside(Rs+rd-vs*t);
            pout=p0/2*(1-vs*t./rd).*heaviside(rd-vs*t).*heaviside(Rs-rd+vs*t);
            S(j,i,:)=pin1+pinr+pout;
       else
            S(j,i,:)=zeros(1,l);
            continue
        end
    end
end
toc

%% SIGNAL PROCESSING
S(:,202:end,:)=flip(S(:,1:200,:),2); %symmetry
S(102:end,:,:)=flip(S(1:100,:,:),1); %symmetry

S_f=zeros(Ny,Nx,l);           
fs=1/(t(2)-t(1)); %sampling freq
fsmin=220e6;
df=fs/l;
f = (0:df:(fs/2));    
fc=[10E6 120E6]; %cut frequency
wn=fc/(fs/2);  
[coefb1,coefa1] = butter(2,wn,'bandpass'); % band pass filter

for i=1:Nx
    S_f(:,i,:)=(filter(coefb1,coefa1,squeeze(S(:,i,:))'))';
end
for i=1:Ny
    S_f(i,:,:)=(filter(coefb1,coefa1,squeeze(S(i,:,:))'))';
end

%% B PLANES VISUALIZATION 
% change 'XZ' por 'YZ' to change the visualized plane
% change 0 to 'doble' to visualise the filteres B planes and the normal ones
% figure(1);set(gcf, 'WindowState', 'maximized'); 
% slider(S,'doble',S_f,Nx,Ny,l,t,dx,dy,'XZ');
% 
% figure(2);set(gcf, 'WindowState', 'maximized');
% slider(S,0,S_f,Nx,Ny,l,t,dx,dy,'YZ');

figure(3)
imagesc(imrotate(squeeze(S_f(:,201,:)),-90)); colorbar; colormap('gray');
xlabel('Nx:0 mm'); title('Plano YZ');
figure(4)
imagesc(imrotate(squeeze(S_f(101,:,:)),-90)); colorbar; colormap('gray');
xlabel('Ny:0 mm'); title('Plano XZ');

%% PARAMETERS
%%%%  AXIAL  %%%%%
fs=1/(t(2)-t(1)); %sampling freq
B=110e6;
res_axial=0.88*vs/B
fs>2*B %1 for Nyquist 

%% SAVING DATA
% saveFolderData = ''; %SAVE FOLDER 
% fileName = datestr(now, 'yyyymmddHHMMSS');

% fileName2   = [ fileName '_SAVE_NAME.mat'];
% save([saveFolderData fileName2], 'VARIABLE_NAME');

% fileName2   = [ fileName '_SAVE_NAME.png'];
% saveas(figure(3), [saveFolderData fileName2]);
% 
% fileName2   = [ fileName '_SAVE_NAME.png'];
% saveas(figure(4), [saveFolderData fileName2]);