clc
clear
close all
R=82.06;
y=[0.162 0.068 0.656 0.114];
p=1;
cphi1=ones(1,4);gama=ones(1,4);
load a.dat
load b.dat
load c.dat
Tsat=b./(a-log(p))-c;
T=sum(y.*Tsat);
psat = exp(a-(b./(T+c)));
psat(1)=p*sum(((y.*cphi1)./gama).*(psat(1)./psat));
T = b(1)/(a(1)-log(psat(1)))-c(1);
psat=exp(a-b./(T+c));
cphi1=cphi(y,T,p);
x=(cphi1.*p.*y)./(gama.*psat);
x=x./sum(x);
gama=gamaunifac(x,T);
psat(1)=p*sum(((y.*cphi1)./gama).*(psat(1)./psat));
T =b(1)/(a(1)-log(psat(1)))-c(1);
psat =exp(a-b./(T+c));
cphi1=cphi(y,T,p);
x=(cphi1.*p.*y)./(gama.*psat);
x=x./sum(x);
gama2=gamaunifac(x,T);
T2=b(1)/(a(1)-log(psat(1)))-c(1);
while (abs(T2-T) ~= 0)
    psat = exp(a-b./(T2+c));
    cphi1=cphi(y,T,p);
    while (abs(gama2-gama) ~= 0)
        x = (cphi1.*P.*y)./(gama.*psat);
        x = x./sum(x);
        gama2=gamaunifac(x,T);
        gama=gama2;
    end
    psat(1) = P*sum(((y.*cphi1)./gamma).*(psat(1)./psat));
    T2 = b(1)/(a(1)-log(psat(1)))-c(1);
    
end
fprintf('the Dew point temprature is %.5f \n',T2)
fprintf('the liquid phase mole fractions are %.4f \n',y)