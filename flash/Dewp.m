function [p,x]=Dewp(y,T)
cphi1(1:length(y))=1;gama(1:length(y))=1;
errgama=1;errp=1;
psat=antoine(T);
p=1/(sum(y.*cphi1./(gama.*psat)));
x=(cphi1.*p.*y)./(psat.*gama);
gama1=gamaunifac(x,T);
p1=1/(sum(y.*cphi1./(gama1.*psat)));
while errp>eps
cphi1=cphi(y,T,p1);
while errgama>eps
x=(cphi1.*p1.*y)./(psat.*gama1);
x=x/sum(x);
gama2=gamaunifac(x,T);
errgama=(gama2-gama1)/gama2;
gama1=gama2;
end
p2=1/(sum(y.*cphi1./(gama.*psat)));
errp=(p2-p1)/p2;
p1=p2;
end
end