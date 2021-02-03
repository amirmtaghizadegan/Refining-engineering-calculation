function [p1,y]=BUBP(x,T)
err=1;cphi1=ones(1,4);
psat=antoine(T);
gama=gamaunifac(x,T);
p1=sum(x.*gama.*psat./cphi1);
while err>=eps
y=x.*gama.*psat./cphi1./p1;
cphi1=cphi(y,T,p1);
p2=sum(x.*gama.*psat./cphi1);
err=(p2-p1)/p2;
p1=p2;
end
end