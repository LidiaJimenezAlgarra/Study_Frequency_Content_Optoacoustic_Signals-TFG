%% SIGNAL ACQUISITION OF OPTOACOUSTIC SIGNALS
% this code computes the acquisition of optoacoustic signals emited for an
% array of elements (N) in a line between a and b. After this, a
% visualization of the layout with the numerical aperture, signal emitors
% and the transducer. Again, it will show the signal with its FFT after filtering them. 

%% DATA
close all
clear all

vs=1500; 
p0=1; %initial pressure
Rs=10e-6; %sphere radius
N=100; %number of sources
l=7000; %length of time signal

a=-1e-3; 
b=1e-3;
z= a+(b-a).*rand(1,N); %random positions
pos_fuentes=zeros(3,N); %x, y, z 
pos_fuentes(3,:)=z; 
t=linspace(0,6*1e-6,l); %seg, time axis
ceros=zeros(1,l); 

%% DETECTOR POSITIONS(X,0,Z)
step=5; %º step
angulo= 0:step:360-step;
rad=1.1e-3; %radius of rotation

pos_det=zeros(3,360/step); 
pos_det(1,:)=rad*cosd(angulo);
pos_det(3,:)=rad*sind(angulo)+mean(z); 

angulo=deg2rad(angulo);
m2up=tan(angulo-pi/6);
m2down=tan(angulo+pi/6);

limup=-(m2up.*pos_det(1,:)-pos_det(3,:))+pos_det(3,:); %NA limit
limdown=-(m2down.*pos_det(1,:)-pos_det(3,:))-pos_det(3,:);
limdown(14:20)=-max(limup); %20 arbitrary

%% Arithmetic corrections
punt_max=max(pos_fuentes(3,:)); 
punt_min=min(pos_fuentes(3,:)); 
for i=0:360/step-1
    if limup(i+1)>punt_max
        limup(i+1)=punt_max;
    end
    if limdown(i+1)<punt_min
        limdown(i+1)=punt_min;
    end
    if limup(i+1)<punt_min
        limup(i+1)=punt_min;
    end
    if limdown(i+1)>punt_max
        limdown(i+1)=punt_max;
    end
    if i==90/step 
        limup(i+1)=punt_max;
        limdown(i+1)=punt_min;
    end
end

%% RAW DATA ACQUISITION
pt=zeros(1,l); 
pgiros=zeros(360/step,l);
for j=1:360/step/4+1 %90 degree
    for i=1:N 
        if (pos_fuentes(3,i)<limup(j) && pos_fuentes(3,i)> limdown(j)) || (pos_fuentes(3,i)<limdown(j) && pos_fuentes(3,i)>limup(j)) % NA condition 
               modulo=sqrt((pos_fuentes(1,i)-pos_det(1,j))^2+(pos_fuentes(3,i)-pos_det(3,j))^2);
               pin1=p0/2*(1+vs*t./modulo).*heaviside(modulo+vs*t).*heaviside(Rs-modulo-vs*t); %optoacoustic signal
               pinr=p0/2*(1-vs*t./modulo).*heaviside(-modulo+vs*t).*heaviside(Rs+modulo-vs*t);
               pout=p0/2*(1-vs*t./modulo).*heaviside(modulo-vs*t).*heaviside(Rs-modulo+vs*t);
               p=pin1+pinr+pout;
               pt=p+pt;
            else
                continue
         end   
     end
        pgiros(j,:)=pt;
        pt=0;     
        display(['Data acquisition angle: ' num2str(j*5) 'º'])
    end

%% DATA VISUALIZATION
figure('WindowState', 'maximized')
pause(0.5)
for i=1:19
    subplot(1,2,1), plot(t*1e3*vs,pgiros(i,:)/max(pgiros(i,:)));xlim([0 3]); ylim([-1.25 1])
    xlabel('mm');ylabel('Normalized pressure');title(['Pressure in angle:' num2str(step*i-step) 'º']);
    subplot(1,2,2),plot(linspace(0,pos_det(1,i).*1e3),linspace(limdown(i).*1e3,pos_det(3,i).*1e3),'r'); %posición límite bajo
    hold on
    plot(linspace(0,pos_det(1,i).*1e3),linspace(limup(i).*1e3,pos_det(3,i).*1e3),'g'); %posición límite verde
    plot(linspace(0,pos_det(1,i).*1e3),linspace(0,pos_det(3,i).*1e3),'b'); %posición centro
    stem(pos_det(1,i).*1e3,pos_det(3,i).*1e3); xlim([-1.5 1.5]); ylim([-1.5 1.5])%xlim([min(pos_det(1,:)-rad/2) max(pos_det(1,:)+rad/2)]);ylim([min(pos_det(3,:)-rad) max(pos_det(3,:)+rad)]); 
    title(["Layout in angle:" num2str(i*step-step) 'º']);xlabel("mm");ylabel("mm") 
    plot(pos_fuentes(1,:).*1e3,pos_fuentes(3,:).*1e3); %xlim([min(pos_det(1,:)) max(pos_det(1,:))]);ylim([min(pos_det(3,:)) max(pos_det(3,:))]);
    plot(ceros,linspace(limdown(i).*1e3,limup(i).*1e3,l),'m'); %línea de señal detectada
    hold off
    pause      
end

%% DATA FILTERING
fs=1/(t(2)-t(1)); %sampling freq
df=fs/length(t);  
dt=1/fs;    % s
f=0:df:fs/2;  % Hz,frequency points

fc=[10E6 120E6]; %cut frequencies
wn=fc/(fs/2);  
[coefb1,coefa1] = butter(2,wn,'bandpass'); %4th degree butterworth bandass filtering

sen_filt=zeros(360/step,l);
sen_filt(1,:)=filter(coefb1,coefa1,pgiros(1,:));
sen_filt(360/step/4+1,:)=filter(coefb1,coefa1,pgiros(360/step/4+1,:)); 

figure(4)
plot(t*vs*1e3,sen_filt(360/step/4+1,:)/max(sen_filt(360/step/4+1,:)));hold on
plot(t*vs*1e3,sen_filt(1,:)/max(sen_filt(1,:)));
legend('Perpendicular',"Parallel");title("Filtered signals");xlabel("mm");ylabel("Normalized pressure");xlim([0 3])
hold off

X_filt=zeros(360/step,l);
X_filt(1,:)=abs(fft(sen_filt(1,:)));
X_filt(360/step/4+1,:)=abs(fft(sen_filt(360/step/4+1,:)));
X_filt=X_filt(:,1:l/2+1); 

figure(5)
plot(f/1E6,X_filt(360/step/4+1,:)/max(X_filt(360/step/4+1,:)));xlim([0 130]); hold on
plot(f/1E6,X_filt(1,:)/max(X_filt(1,:)));
legend('Perpendicular',"Parallel");title("FFT normalized");xlabel("Frequency (MHz)")
hold off

%% DATA COMPARISONS
figure('WindowState', 'maximized')
X_filt=zeros(360/step,l);
for i=1:19
    sen_filt(i,:)=filter(coefb1,coefa1,pgiros(i,:));
    subplot(1,2,1),plot(t*vs*1e3,sen_filt(i,:)/max(sen_filt(i,:)));
    title(["Filtered signal in angle:" num2str(i*step-step) 'º']);xlabel("mm");ylabel("Normalized pressure");xlim([0 3])
    title(["FFT in angle:" num2str(i*step-step) 'º']);xlabel("mm");ylabel("mm") 

    X_filt(i,:)=abs(fft(sen_filt(i,:)));
    subplot(1,2,2),plot(f/1E6,X_filt(i,1:l/2+1)/max(X_filt(i,:))); 
    title(["FFT in angle:" num2str(i*step-step) 'º']);xlabel("MHz");ylabel("au");xlim([0 150])  
    
    pause(1)      
end