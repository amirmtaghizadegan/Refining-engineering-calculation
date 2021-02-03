function gama=gamaunifac(x,T)
nu=0;
load data.mat
r=zeros(1,2);q=zeros(1,2);e=zeros(nk,length(k));
tow=zeros(length(k));beta=zeros(nk,length(k));
theta=zeros(1,length(k));s=zeros(1,length(k));
for i=1:nk
    r(i)=sum(nu(i,:).*R);
    q(i)=sum(nu(i,:).*Q);
    for j=1:length(k)
        e(i,j)=nu(i,j)*Q(j)/q(i);
    end
end
for i=1:length(k)
    for j=1:length(k)
        tow(i,j)=exp(-a(i,j)/T);
    end
end
for ii=1:nk
    for kk=1:length(k)
        beta(ii,kk)=sum(e(ii,:).*tow(:,kk)');
    end
end
for kkk=1:length(k)
    theta(kkk)=sum((x'.*q'.*e(:,kkk)))./sum((x'.*q'));
end
for i=1:length(k)
    s(i)=sum(theta.*tow(:,i)');
end
j=r./sum(r.*x);
l=q./sum(q.*x);
gamac=1-j+log(j)-5*q*(1-j/l+log(j/l));
gamar=q.*(1-sum(transpose(theta.*beta./s-e.*log(beta./s))));
gama=exp(gamar+gamac);
end