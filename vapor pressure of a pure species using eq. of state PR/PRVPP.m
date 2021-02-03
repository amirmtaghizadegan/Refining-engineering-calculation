clear
close all
clc
%input T: 423
%input T critical: 540.2
%input P critical: 27.40
%input omega of the specie: 0.350
sigma=1+sqrt(2);epsilon=1-sqrt(2);ohm=0.07780;psi=0.45724;zc=0.30740;
T=input('input T: ');
Tc=input('input T critical: ');
Pc=input('input P critical: ');
w=input('input omega of the specie: ');
p=input('input the initial p');
Tr=T/Tc;
alpha=(1+(0.37464+1.54226*w-0.26992*w^2)*(1-Tr^0.5))^2;
R=83.14;
a=psi*alpha*(R^2*Tc^2)/Pc;
b=ohm*R*Tc/Pc;
beta=b*p/R/T;
q=a/b/R/T;
Zl=Zliq(beta,epsilon,sigma,q);
Zv=Zvap(beta,q,epsilon,sigma);
Il=1/(sigma-epsilon)*log((Zl+sigma*beta)/(Zl+epsilon*beta));
Iv=1/(sigma-epsilon)*log((Zv+sigma*beta)/(Zv+epsilon*beta));
lnphil=Zl-1-log(Zl-beta)-q*Il;
lnphiv=Zv-1-log(Zv-beta)-q*Iv;
err=1;
while err>eps
    p=p*(exp(lnphil)/exp(lnphiv));
    beta=b*p/R/T;
    Zl=Zliq(beta,epsilon,sigma,q);
    Zv=Zvap(beta,q,epsilon,sigma);
    Il=1/(sigma-epsilon)*log((Zl+sigma*beta)/(Zl+epsilon*beta));
    Iv=1/(sigma-epsilon)*log((Zv+sigma*beta)/(Zv+epsilon*beta));
    lnphil=Zl-1-log(Zl-beta)-q*Il;
    lnphiv=Zv-1-log(Zv-beta)-q*Iv;
    err=lnphil-lnphiv;
end
pexp1=pexp(T);
fprintf('your psat calculated by PR eq. of state is: %.4f \n',p)
fprintf('your psat calculated by perrys handbook coeffs is: %.4f\n',pexp1)
fprintf('the deviation between psats is: %.4f \n',abs((pexp1-p)/pexp1))
%your psat calculated by PR eq. of state is: 3.7157 
%your psat calculated by perrys handbook coeffs is: 3.6832
%the deviation between psats is: 0.0088 input