%% 1D acquisition of optoacoustic signals in a 2D grid  %%
% this code computes an array of sources (N) in random positions (y) between
% interval(1) and interval(2).
% the length of the time-resolved signals is l and the spheres are assumed
% of a diameter (Rs) of 10um.
% Later, it filters the generated signals and displays their FFT 

%% PARAMETERS
clear all;
close all;
tic
N=1e6; %number of sources
vs=1520;  %m/s
p0=1; %u.au
Rs=10e-6; %m diameter of sphere

pdet=zeros(2,3);
l=100000; %length t axis
c=zeros(1,l);
%% Transducer positions (Z,Y,0)       
pdet(1,:)=[0,1500*1e-6,0]; %0 degrees 
pdet(2,:)=[0,0,1000*1e-6];  %90 degrees 
intval=[-1000*1e-6,1000*1e-6 ]; 
a=intval(1);
b=intval(2);
y=intval(1) + (intval(2)-intval(1)).*rand(N,1);

t=linspace(0,4*1e-6,l); %T AXIS
pp=zeros(2,length(t)); 

rup=pdet(2,3)*sind(30); %NA LIMITS
rdown=pdet(2,3)*sind(-30);

%% RAW DATA ACQUISITION
for j=1:2
    for i=1:N  
        if y(i)<=rup && y(i)>=rdown && j==2 %%parallel
        psor=[0,y(i),0]; %source position
        rd=sqrt((psor(1)-pdet(j,1))^2+(psor(2)-pdet(j,2))^2+(psor(3)-pdet(j,3))^2); %distance
        pin1=p0/2*(1+vs*t./rd).*heaviside(rd+vs*t).*heaviside(Rs-rd-vs*t);
        pinr=p0/2*(1-vs*t./rd).*heaviside(-rd+vs*t).*heaviside(Rs+rd-vs*t);
        pout=p0/2*(1-vs*t./rd).*heaviside(rd-vs*t).*heaviside(Rs-rd+vs*t);
        p=pin1+pinr+pout;
        pp(j,:)=pp(j,:)+p;    
        elseif j==1 %perpendicular
        psor=[0,y(i),0];%source position
        rd=sqrt((psor(1)-pdet(j,1))^2+(psor(2)-pdet(j,2))^2+(psor(3)-pdet(j,3))^2); %distance
        pin1=p0/2*(1+vs*t./rd).*heaviside(rd+vs*t).*heaviside(Rs-rd-vs*t);
        pinr=p0/2*(1-vs*t./rd).*heaviside(-rd+vs*t).*heaviside(Rs+rd-vs*t);
        pout=p0/2*(1-vs*t./rd).*heaviside(rd-vs*t).*heaviside(Rs-rd+vs*t);
        p=pin1+pinr+pout;
        pp(j,:)=pp(j,:)+p;   
        else
            continue
        end
    end
end
%% FIGURES
figure(1);
plot(t*vs/1e-3,pp(1,:)/max(pp(1,:))); 
hold on;
plot(t*vs/1e-3,pp(2,:)/max(pp(2,:))); 
legend('Perpendicular', 'Parallel'); xlim([0 3]);
xlabel('mm');title('Received signal');ylabel(['a.u'])
%% SIGNAL FILTERING AND FFT   
fs=1/(t(2)-t(1)); %sampling freq
df=fs/length(t);
f = (0:df:(fs/2));    

fp=zeros(size(pp)); %FFT
fp(1,:)=fft(pp(1,:)); 
fp(2,:)=fft(pp(2,:));

figure (2);
plot(f*1e-6,abs(fp(1,1:length(t)/2+1))/max(abs(fp(1,1:length(t)/2+1)))); 
hold on; 
plot(f*1e-6,abs(fp(2,1:length(t)/2+1))/max(abs(fp(2,1:length(t)/2+1)))); 
legend('Perpendicular', 'Parallel');
xlabel('MHz'); title('FFT of sigals');xlim([0 150]);

fc=[10E6 120E6]; %cut frequencies
wn=fc/(fs/2);  
[coefb1,coefa1] = butter(2,wn,'bandpass'); %band pass filter

%% FILTERED SIGNAL
sen_filt=zeros(size(pp));
sen_filt(1,:)=filter(coefb1,coefa1,pp(1,:));
sen_filt(2,:)=filter(coefb1,coefa1,pp(2,:));

figure(3)
plot(t*vs/1e-3,sen_filt(1,:)/max(sen_filt(1,:))); hold on;
plot(t*vs/1e-3,sen_filt(2,:)/max(sen_filt(2,:))); 
legend('Perpendicular', 'Parallel');
xlabel('mm');title('Signal filered'); xlim([0 3]); ylim([-1 1])
hold off
ffilt=zeros(size(sen_filt));
ffilt(1,:)=abs(fft(sen_filt(1,:)));
ffilt(1,:)=ffilt(1,:)/max(ffilt(1,:));
ffilt(2,:)=abs(fft(sen_filt(2,:)));
ffilt(2,:)=ffilt(2,:)/max(ffilt(2,:));

figure (4);
plot(f*1e-6,movmean(ffilt(1,1:length(t)/2+1),6)/max(movmean(ffilt(1,1:length(t)/2+1),6))); 
hold on;
plot(f*1e-6,ffilt(2,1:length(t)/2+1)); 
legend('Perpendicular', 'Parallel');
xlabel('MHz'); title('FFT filtered') ;xlim([0 150]);
hold off

%% Layout drawing
figure
plot(1e3*linspace(intval(1),intval(2),100),zeros(1,100),'b','Linewidth',4);hold on
plot(pdet(1),pdet(2,3)*1e3,'--o','MarkerSize',5,'MarkerFaceColor','m'); title('Layout');
plot(pdet(1,2)*1e3,pdet(1),'--o','MarkerSize',5,'MarkerFaceColor','g');
plot(1e3*linspace(rup,rdown),zeros(1,100),'r','Linewidth',4);
plot(1e3*linspace(rup,0),1e3*linspace(0,pdet(2,3)),'g');
plot(1e3*linspace(rdown,0),1e3*linspace(0,pdet(2,3)),'g');
legend('Sources','Parallel',"Perpendicular",'NA condition'); xlim([-2 2]);ylim([-1.5 1.5])

hold off
toc
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