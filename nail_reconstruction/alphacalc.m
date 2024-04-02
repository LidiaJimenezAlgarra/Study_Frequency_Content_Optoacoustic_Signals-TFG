function [ alphac] = alphacalc(projlow,projhigh)
 alphaval=[0:0.01:100];
 differ=zeros(size(alphaval));
 for(i=1:length(alphaval))
 differ(i)=sum(sum((projlow - alphaval(i)*projhigh).^2));
 end
  [dum,ind]=min(differ);
  alphac=alphaval(ind);
%figure; plot(alphaval,differ);
end

