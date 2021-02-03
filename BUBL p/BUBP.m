clear
close all
clc
load data.mat
n=4;
err=1;T=334.15;x=[0.162 0.068 0.656 0.114];
psat=antoine(T);
gama=gamaunifac(x,T);
cphi1=ones(1,4);
p1=sum(x.*gama.*psat./cphi1);
while err>=eps
y=x.*gama.*psat./cphi1./p1;
cphi1=cphi(y,T,p1);
p2=sum(x.*gama.*psat./cphi1);
err=abs((p2-p1)/p2);
p1=p2;
end
fprintf('p equals to: %d\n',p1)
disp('y of each article is: ')
disp(y)