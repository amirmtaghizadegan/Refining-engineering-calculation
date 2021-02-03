function [n,T,x,a,b,c]=inputbublp()
n=input('how many articles are in system? ');
T=input('enter T: ');
x=zeros(1,n);a=zeros(1,n);b=zeros(1,n);c=zeros(1,n);
for i=1:n
fprintf('input x for article %d: ',i)
x(i)=input('');
end
for i=1:n
fprintf('input A,B,C(antoine eqs. constants) for fraction %d:\n ',i)
a(i)=input('enter a: ');
b(i)=input('enter b: ');
c(i)=input('enter c: ');
end