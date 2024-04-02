%% Layout representation of 2D and 3D grids with a cross reconstruction and cillinder
% you can create figures and visualize them

close all

%% VARIABLES
l=500; %length signal
dx=100e-6; %m
dy=100e-6; %m
pdetX=(-2e-3:dx:2e-3); % X (m) -2 a 2 mm
pdetY=(-1e-3:dy:3e-3); % Y (m) -1 a 1 mm 
Nx=length(pdetX);
Ny=length(pdetY);
N=500;
pfuente=zeros(3,N);% m x, y, z

%vertical
a=-0.75e-3; 
b=-1.25e-3;
z= a+(b-a).*rand(1,N/2+1); %random numbers
pfuente(3,1:N/2+1)=z;

%horizontal
a= 0.25e-3;
b=-0.25e-3;
x= a+(b-a).*rand(1,N/2); %random numbers
pfuente(1,N/2+1:end)=x;
pfuente(3,N/2+1:end)=-1e-3;

%% 
an = 0:pi/50:2*pi;
r=1e-3;
punto=zeros(3,1);
punto(1)=-0.25e-3;
punto(2)=-0.5e-3;
punto(3)=0;
x= r * cos(an)-0.25e-3;
y= r * sin(an)-0.5e-3;
posNA=zeros(3,length(an));
posNA(1,:)=x;
posNA(2,:)=y;
posNA(3,:)=-1.25e-3;

fcilx=zeros(3,100); 
fcily=zeros(3,100); 
fcilz=zeros(3,100); 

fcilx(1,:)=linspace(-0.5e-3,0.5e-3,100); %tamaño cilindro x
fcilx(3,:)=-1e-3;

fcily(2,:)=linspace(-0.75e-3,0.75e-3,100); %tamaño cilindro y
fcily(3,:)=-1e-3;

fcilz(3,:)=linspace(-0.75e-3,-1.25e-3,100); %tamaño cilindro z


%% GRID
S=zeros(3,Nx*Ny);
N=Nx*Ny;
for i=Nx:Nx:N
    S(1,i-Nx+1:i)=pdetX;
    if pdetY(i/Nx)>1e-3
       S(2,i-Nx+1:i)=0;
    else
    S(2,i-Nx+1:i)=pdetY(i/Nx);
    end
end

%% ANOTHER GRID
dx=100e-6; %m
dy=100e-6; %m
pdetX=(-2e-3:dx:2e-3); % X (m) -2 a 2 mm
pdetY=(-1e-3:dy:3e-3); % Y (m) -1 a 1 mm 

S_n=zeros(3,231);
N=231;
Nx=21;
for i=Nx:Nx:N
    S_n(1,i-Nx+1:i)=pdetX(1:21);
    if pdetY(i/Nx)>1e-3
       S_n(2,i-Nx+1:i)=0;
    else
    S_n(2,i-Nx+1:i)=pdetY(i/Nx);
    end
end

% S_n(1,:)=linspace(-2e-3,0,121);
% S_n(2,:)=linspace(-1e-3,0,121);
%% 3D layout
figure(3)
scatter3(S(1,:)*1e3,S(2,:)*1e3,S(3,:)*1e3,'blue');title('3D layout');xlabel('X axis (mm)');ylabel('Y axis (mm)');zlabel('Z axis (mm)'); hold on
%scatter3(S_n(1,:)*1e3,S_n(2,:)*1e3,S_n(3,:)*1e3,'red','filled');%grid
%scatter3(pfuente(1,:)*1e3,pfuente(2,:)*1e3,pfuente(3,:)*1e3);%cross
scatter3(fcilx(1,:)*1e3,fcilx(2,:)*1e3,fcilx(3,:)*1e3,'red','filled');%cillinder
scatter3(fcily(1,:)*1e3,fcily(2,:)*1e3,fcily(3,:)*1e3,'green','filled');%cillinder
scatter3(fcilz(1,:)*1e3,fcilz(2,:)*1e3,fcilz(3,:)*1e3,'yellow','filled');%cillinder

%scatter3(posNA(1,:)*1e3,posNA(2,:)*1e3,posNA(3,:)*1e3,'filled','magenta');
%scatter3(linspace(-0.25,0,101)*1e3,linspace(-0.5,0,101)*1e3,posNA(3,:)*1e3,'filled','magenta');

%scatter3(punto(1)*1e3,punto(2)*1e3,punto(3)*1e3,'filled','black');hold off%detector position

% %
% line(-0.25*ones(100),linspace(0.5,-0.5,100),linspace(0,-2,100),'Color','green','LineWidth',4);hold off
legend('Grid','Cilindro X','Cilindro Y','Cilindro Z'); ylim([-1.5 1.5]); xlim([-2.5 2.5]);zlim([ -1.5 0])
%legend('Grid','Real scanning points','Sources','NA','Detector'); ylim([-1.5 1.5]); xlim([-2.5 2.5])

%% 2D layout
tic
N=100000; %number of sources
vs=1500;  %m/s
p0=1; %u.au
Rs=10e-6; %m
pdet=zeros(2,3);
l=100000;
c=zeros(1,l);


pdet(1,:)=[0,1500*1e-6,0]; %0 degrees parallel
pdet(2,:)=[0,0,1000*1e-6];  %90 degrees perpendicular
intval=[-1000*1e-6,1000*1e-6 ]; 
a=intval(1);
b=intval(2);
y=intval(1) + (intval(2)-intval(1)).*rand(N,1);
t=linspace(0,4*1e-6,l); %time
pp=zeros(2,length(t)); %pressure

rup=pdet(2,3)*sind(30);%NA limits
rdown=pdet(2,3)*sind(-30);
%% PARALLEL LAYOUT
figure(1)
plot(y*1e3,c,'LineWidth',3);xlim([-2 2]); ylim([-2 2]); 
xlabel('mm','FontName','Helvetica');ylabel('mm','FontName','Helvetica');hold on
scatter(pdet(1),pdet(2,3)*1e3,'MarkerFaceColor','k','MarkerEdgeColor','k');title('2D Parallel simulation')
plot(1e3*linspace(rup,rdown),zeros(1,100),'r','Linewidth',4); %NA
plot(1e3*linspace(rup,0),1e3*linspace(0,pdet(2,3)),'g');
plot(1e3*linspace(rdown,0),1e3*linspace(0,pdet(2,3)),'g');
legend({'Sources','Detector','Sources detected with NA','NA condition'}, 'FontSize',9,'FontName','Helvetica');
set(gca, 'YTick', [-2 -1 0 1 2]); % EJE Y 
set(gca, 'XTick', [-2 -1 0 1 2]); % EJE X
hold off

%% PERPENDICULAR LAYOUT
figure(2)
plot(y*1e3,c,'LineWidth',3);hold on
plot(y*1e3,c,'LineWidth',3,'Color','r');xlabel('mm','FontName','Helvetica');ylabel('mm','FontName','Helvetica');
xlim([-2 2]); ylim([-2 2]);
scatter(pdet(1,2)*1e3,pdet(1),'MarkerFaceColor','k','MarkerEdgeColor','k');title('2D perpendicular simulation');
plot(linspace(0.5,1.5),linspace(0.5,0),'g');plot(linspace(0.5,1.5),linspace(-0.5,0),'g');
legend({'Sources','Sources detected with NA','Detector','NA condition'}, 'FontSize',9,'FontName','Helvetica');
set(gca, 'YTick', [-2 -1 0 1 2]); % EJE Y 
set(gca, 'XTick', [-2 -1 0 1 2]); % EJE X
hold off
