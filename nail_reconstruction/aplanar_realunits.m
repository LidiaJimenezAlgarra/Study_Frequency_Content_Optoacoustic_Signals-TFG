
function [R] = aplanar_realunits(R,ds,dz,zr,yr,surffile) 
 load(surffile);
 
vs=1525;
dt = 1/Fs;
dzsurf= vs*dt*1e6;

% load('R_!20150226125036_patient11_fullscan_wl1_corr_Freq10-30MHz');

%y=[yr(1):reconParams.GRID_DS:yr(end)];
%z=[zr(1):reconParams.GRID_DZ:zr(end)];
y=[yr(1):ds:yr(end)];
z=[zr(1):dz:zr(end)];

[Xq,Yq] = meshgrid(linspace(1,size(surfSmooth,2),size(R,2)),linspace(1,size(surfSmooth,1),size(R,3)));

Surface=interp2(surfSmooth,Xq,Yq);

ref=Surface(round(size(Surface,1)/2),round(size(Surface,2)/2));
%ref=Surface(1,1);
for(i=1:size(R,2))
    for(j=1:size(R,3))
    dist=(Surface(j,i)-ref)*dzsurf;
    ind=round(-dist/dz);
    %ind=Surface(i,j)-ref;
   % ind=double(ind);
    %[dum,ind]=min(abs(z-surfSmooth(j-3,i)));
    a=R(:,i,j);
   % indref-ind
   if(ind ~=0)
        a = circshift(a,round(ind));
   end
        R(:,i,j)=a;
       % R(ind,i,j)=1000;
    end
end

 %Rlow=R(70:450,:,:);
 %save('R_!20150226125036_patient11_fullscan_Freq10-30MHz_plane','R','reconParams','xr','shiftInd','xr','yr','zr');
 %Rlow=R;
 %load('R_!20150226125036_patient11_fullscan_Freq30-120MHz');

end 
  
  
  
  
 %  C=squeeze(max(R(a:end,:,:),[],1));
% 
%  a=b;
%  b=250;
%  figure(4);
%  imagesc(yr,xr,squeeze(max(R(a:b,:,:),[],1))), daspect([1 1 1]), colorbar, colormap hot;
%  
%   
%  %%%%%%%%%%%%%%%%%%%%%%%%%
% % R=RR;
%   R=R/(max(max(max(abs(R)))));
% 
%  
%  %R=Rhigh/(max(max(max(abs(Rhigh)))))+0.5*Rfull/(max(max(max(abs(Rfull)))));
%   figure(5);
%  R=R/(max(max(max(abs(R)))));
%  imagesc(xr,zr,squeeze(max(R(:,:,:),[],3))), daspect([1 1 1]), colorbar, colormap hot;
% 
%  a=90;
%  b=120;
%   figure(6);
%  imagesc(yr,xr,squeeze(max(abs(R(a:b,:,:)),[],1))), daspect([1 1 1]), colorbar, colormap gray;
%   C=squeeze(max(R(a:b,:,:),[],1));
% % imwrite(C,'psoriasis_dermalpap','jpeg'); 
%  a=b;
% 
%  figure(7);
%  imagesc(yr,xr,squeeze(max(abs(R(a:end,:,:)),[],1))), daspect([1 1 1]), colorbar, colormap gray;
%  %  C=squeeze(max(R(a:end,:,:),[],1));
%  imwrite(C,'psoriasis_below','jpeg'); 
%    
%  a=b;
%  b=250;
%  figure(8);
%  imagesc(yr,xr,squeeze(max(R(a:b,:,:),[],1))), daspect([1 1 1]), colorbar, colormap gray;
%  
%  
%  
%  
%  
%  
%  
%  
%  
%  
%  
%  for(i=1:size(R,1))
%   a= 'foramirahealth\';
%     a=[a num2str(i)];
%     imwrite(squeeze(R(i,:,:)),a,'jpeg'); 
%  end
%  
%  
% % FrangiFilter2D(I, options) 
%  