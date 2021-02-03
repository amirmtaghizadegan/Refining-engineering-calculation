clear
close all
clc
T=334.15;p=1;z=[0.250 0.400 0.200 0.150];n=4;
[pbubl,y2]=BUBP(z,T);
[pdew,x2]=Dewp(z,T);
if p>pdew && p<pbubl
gamma=gamaunifac(x2,T);
cphi1=cphi(y2,T,p);
nu2=(pbubl-p)/(pbubl-pdew);
x1=zeros(1,n);y1=zeros(1,n);
psat=antoine(T);
nu1=0;
while abs(nu2-nu1)>eps || abs(sum((y2-y1)./y2))>eps || abs(sum((x2-x1)./x2))>eps
x1=x2;y1=y2;nu1=nu2;
k=gamma.*psat./cphi1/p;
  Fx=sum(z./(1+nu1.*(k-1)))-1;
  Fy=sum((z.*k)./(1+nu1.*(k-1)))-1;
  F = Fy-Fx;
  DfDv=-sum(((z.*k-1).^2)./((1+nu1*(k-1)).^2));
  nu2=nu1-(F/DfDv);
x2=z./(1+nu2.*(k-1));
y2=k.*x2;
gamma=gamaunifac(x2,T);
cphi1=cphi(y2,T,p);
end
end
rownames={'n-hexane(1)','ethanol(2)','MCP(3)','benzene(4)'};
z=z';x=x2';y=y2';k=k';
fprintf('nu= %.4f \n',nu2)
Tab=table(z,x,y,k,'RowNames',rownames);
disp(Tab)
