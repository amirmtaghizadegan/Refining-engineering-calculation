clear
close all
clc

p=input('input p: ');
x=input('input Xs as a line matrix: ');
cphi1=ones(1,length(x));errT=1;
load a.dat
load b.dat
load c.dat
Tsat=b./(a-log(p))-c;
T1=sum(x.*Tsat);
psat=antoine(T1,a,b,c);
gama=gamaunifac(x,T1);
psat(1)=p/sum((x.*gama./cphi1).*(psat./psat(1)));
T1=b(1)./(a(1)-log(psat(1)))-c(1);
while errT>eps
psat=antoine(T1,a,b,c);
y=(x.*gama.*psat)./(cphi1.*p);
cphi1=cphi(y,T1,p);
gama=gamaunifac(x,T1);
psat(1)=p/sum((x.*gama./cphi1).*(psat/psat(1)));
T2=b(1)./(a(1)-log(psat(1)))-c(1);
errT=(T2-T1)/T2;
T1=T2;
end
fprintf('the temprature is : %.4f \n',T1)
disp('y is :')
disp(y)