clear
close all
clc
T=310.8;x1=transpose(0:0.001:1);
x2=[x1 1-x1];pp=zeros(1,size(x2,1));yy=zeros(size(x2,1),2);
for ii=1:size(x2,1)
x=x2(ii,:);
Tc=[190.6 425.1];
Pc=[45.99 37.96];
w=[0.012 0.200];
Tr=T./Tc;
sigma=1;epsilon=0;ohm=0.08664;psi=0.42748;
alpha=(1+(0.480+1.574.*w-0.176.*w.^2).*(1-Tr.^0.5)).^2;
R=83.14;
a=psi*alpha.*(R^2.*Tc.^2)./Pc;
b=ohm*R.*Tc./Pc;
psat=PRVPP(T,Tc,Pc,w);
%psat of methane is unknown and stimated as 46 bar
psat(1)=46;
p=sum(psat.*x);
y=psat.*x/p;check1=2;iii=0;
while abs(check1-1)>0.0001 && iii~=1000
    iii=iii+1;
    bl=sum(b.*x);bv=sum(b.*y);
    betal=bl.*p/R/T;
    betav=bv.*p/R/T;
    %for having a clean calc we initiate with a=0
    al=0;av=0;aa=zeros(length(a));
    for i=[1 2]
        for j=1:2
            aa(i,j)=sqrt(a(i)*a(j));
        end
        al=al+sum(x(i).*x.*aa(i,:));
        av=av+sum(y(i).*y.*aa(i,:));
    end
    ql=al/bl/R/T;
    qv=av/bv/R/T;
    Zl=Zliq(betal,epsilon,sigma,ql);
    Zv=Zvap(betav,qv,epsilon,sigma);
    abarv=zeros(1,length(a));abarl=zeros(1,length(a));
    qbarv=zeros(1,length(a));qbarl=zeros(1,length(a));
    for i=1:length(a)
        abarv(i)=2*sum(y.*sqrt(a(i)*a))-av;
        abarl(i)=2*sum(x.*sqrt(a(i)*a))-al;
        qbarv(i)=ql*(1+abarv(i)/av-b(i)/bv);
        qbarl(i)=ql*(1+abarl(i)/al-b(i)/bl);
    end
    Il=1/(sigma-epsilon)*log((Zl+sigma.*betal)./(Zl+epsilon.*betal));
    Iv=1/(sigma-epsilon)*log((Zv+sigma.*betav)./(Zv+epsilon.*betav));
    phiv=exp(b./bv*(Zv-1)-log(Zv-betav)-qbarv.*Iv);
    phil=exp(b./bl*(Zl-1)-log(Zl-betal)-qbarl.*Il);
    k=phil./phiv;check=0;
    y=k.*x/sum(x.*k);
    check1=sum(k.*x);
    while check1-check>eps
        check=check1;
        bv=sum(b.*y);
        betav=bv.*p/R/T;
        for i=[1 2]
            for j=1:2
                aa(i,j)=sqrt(a(i)*a(j));
            end
            av=av+sum(y(i).*y.*aa(i,:));
        end
        qv=av/bv/R/T;
        for i=1:length(a)
            abarv(i)=2*sum(y.*sqrt(a(i)*a))-av;
            qbarv(i)=ql*(1+abarv(i)/av-b(i)/bv);
        end
        Iv=1/(sigma-epsilon)*log((Zv+sigma.*betav)./(Zv+epsilon.*betav));
        phiv=exp(b./bv*(Zv-1)-log(Zv-betav)-qbarv.*Iv);
        k=phil./phiv;
        y=k.*x/sum(x.*k);
        check1=sum(k.*x);
    end
    p=p*check1;
    
end
pp(ii)=p;
yy(ii,:)=y;
end
plot(x1,pp)
hold on
plot(yy(:,1),pp)
