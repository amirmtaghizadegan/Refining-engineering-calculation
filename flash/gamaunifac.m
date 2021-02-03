function gamma=gamaunifac(x,T)
nu=[2 4 0 0 0;1 1 0 0 1;1 4 1 0 0;0 0 0 6 0];
Rk=[0.9011 0.6744 0.4469 0.5313 1];
Qk=[0.848 0.540 0.228 0.400 1.2];
a=[0 0 0 61.13 986.5;0 0 0 61.13 986.5;0 0 0 61.13 986.5;-11.12 -11.12 -11.12 0 636.1;156.40 156.40 156.40 89.60 0];
tow=exp(-a/T);
r=zeros(1,size(nu,1));q=zeros(1,size(nu,1));
e=zeros(size(nu,1),size(nu,2));beta=zeros(size(nu,1),size(nu,2));
theta=zeros(1,size(nu,2));s=zeros(1,size(nu,2));lngamac=zeros(1,size(nu,1));
for i=1:size(nu,1)
    r(i)=sum(nu(i,:).*Rk);
    q(i)=sum(nu(i,:).*Qk);
    e(i,:)=nu(i,:).*Qk/q(i);
    for k=1:size(nu,2)
        beta(i,k)=sum(e(i,:).*tow(:,k)');
        theta(k)=sum(x.*q.*e(:,k)')/sum(x.*q);
    end
end
for kk=1:size(nu,2)
    s(kk)=sum(theta.*tow(:,kk)');
end
j=r./sum(r.*x);
l=q./sum(q.*x);
for ii=1:size(nu,1)
    lngamac(ii)=(1-j(ii)+log(j(ii))-5*q(ii)*(1-j(ii)/l(ii)+(log(j(ii)/l(ii)))));
end
lngamar=q.*(1-sum(transpose(theta.*beta./s-e.*log(beta./s))));
gamma=exp(lngamac+lngamar);
end