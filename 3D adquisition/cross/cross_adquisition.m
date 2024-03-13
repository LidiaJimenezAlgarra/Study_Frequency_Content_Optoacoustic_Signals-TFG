%% Tomographic acquisition of a cross structure 
% This code computes an array of sources (N) in random positions (z and z) between
% a and b.
% The length of the time-resolved signals is l and the spheres are assumed
% of a diameter (Rs) of 10um.
% Later, the signals are filtered with a 4th degree Butterworth filter and
% a sagital and coronal MIP is computed and shown

%% GRID
clear all 
close all
tic
l=500; %t points 
dx=10e-6; %m
dy=10e-6; %m
vs=1500; %m/s
p0=1; %au
Rs=10e-6; %m sources diameter

pdetX=(-2e-3:dx:2e-3); %positions of detector X
pdetY=(-1e-3:dy:1e-3); %positions of detector Y
Nx=length(pdetX);
Ny=length(pdetY);
Ntot=Nx*Ny; 

%% SOURCES
N=1e4; 
t=linspace(0,1.0667e-06,l); %seg
pfuente=zeros(3,N);% m x, y, z

%vertical
a=-0.75e-3; %
b=-1.25e-3;
z= a+(b-a).*rand(1,N/2+1); %random source positions
pfuente(3,1:N/2+1)=z;

%horizontal
a= 0.25e-3;
b=-0.25e-3;
x= a+(b-a).*rand(1,N/2);  %random source positions
pfuente(1,N/2+1:end)=x;
pfuente(3,N/2+1:end)=-1e-3;

%% SIGNAL VARIABLES
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
        pmax=-sqrt((pdetX(i))^2+(pdetY(j))^2)*sind(60)/sind(30);%NA CONDITION DEPTH
        pr=pdetX(i)+0.5e-3 ; %NA CONDITION SIDES
        pl=pdetX(i)-0.5e-3 ; 
        for s=1:N 
            if s<=N/2
                if  pfuente(3,s)<=pmax
                     rd=sqrt((pdetX(i)-pfuente(1,s))^2+(pdetY(j)-pfuente(2,s))^2+(pfuente(3,s)^2));
                     pin1=p0/2*(1+t1opti./rd).*heaviside(rd+t1opti).*heaviside(-rd+t2opti); 
                     pinr=p0/2*(1-t1opti./rd).*heaviside(-rd+t1opti).*heaviside(rd+t2opti);
                     pout=p0/2*(1-t1opti./rd).*heaviside(rd-t1opti).*heaviside(Rs-rd+vs*t);
                     ptot=pin1+pinr+pout;
                end 
                atemp=atemp+ptot;
                ptot=0;
            else 
                if (pdetY(j)>=-5e-4 && pdetY(j)<=5e-4) && (pfuente(1,s)>=pl && pfuente(1,s)<=pr) && (pdetX(i)>=-1.5e-3 && pdetX(i)<=1.5e-3)
                     rd=sqrt((pdetX(i)-pfuente(1,s))^2+(pdetY(j)-pfuente(2,s))^2+(pfuente(3,s)^2));
                     pin1=p0/2*(1+t1opti./rd).*heaviside(rd+t1opti).*heaviside(-rd+t2opti); 
                     pinr=p0/2*(1-t1opti./rd).*heaviside(-rd+t1opti).*heaviside(rd+t2opti);
                     pout=p0/2*(1-t1opti./rd).*heaviside(rd-t1opti).*heaviside(Rs-rd+vs*t);
                     ptot=pin1+pinr+pout;
                end 
                atemp=atemp+ptot;
                ptot=0;
            end
        end
       S(j,i,:)=atemp;
       atemp=0;
    end
    display(['x:' num2str(i) "mm"])
end
toc
%% SIGNAL PROCESING AND FFT
fs=1/(t(2)-t(1)); % frequency sampling

S(:,202:end,:)=flip(S(:,1:200,:),2); %symmetry 
S(102:end,:,:)=flip(S(1:100,:,:),1); %symmetry

S_f=zeros(Ny,Nx,l);           
df=fs/l;
f = (0:df:(fs/2));    
fc=[10E6 120E6]; %cut frequencies
wn=fc/(fs/2);  
[coefb1,coefa1] = butter(2,wn,'bandpass'); %4th degree butterworth filtering

for i=1:Nx
    S_f(:,i,:)=(filter(coefb1,coefa1,squeeze(S(:,i,:))'))';
end
for i=1:Ny
    S_f(i,:,:)=(filter(coefb1,coefa1,squeeze(S(i,:,:))'))';
end

%% DATA VIEWING
figure(3)
imagesc(pdetY*1e3,t*vs*1e3,squeeze(S_f(:,201,:))'); colorbar; colormap('gray');
xlabel('Nx:0 mm'); title('YZ plane');ylabel('mm')
figure(4)
imagesc(pdetX*1e3,t*vs*1e3,squeeze(S_f(101,:,:))'); colorbar; colormap('gray');ylabel('mm')
xlabel('Ny:0 mm'); title('XZ plane');

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