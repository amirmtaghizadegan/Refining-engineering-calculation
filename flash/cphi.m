function  cphi1=cphi(y,T,p)
jam1=0;
sigma=zeros(length(y));
R=83.14;
load bv.dat
B=bv;psat=antoine(T);cphi1=zeros(1,4);
for i=1:length(y)
    for j=1:length(y)
        sigma(i,j)=2*B(i,j)-B(i,i)-B(j,j);
    end
end
for i=1:length(y)
    for j=1:length(y)
        for k=1:length(y)
            jam2=y(j)*y(k)*(2*sigma(j,i)-sigma(j,k));
            jam1=jam1+jam2;
        end
    end
    cphi1(i)=exp((B(i,i)*(p-psat(i))+0.5*p*jam1)/(R*T));
    jam1=0;
end
end
