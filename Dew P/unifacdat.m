function unifacdat()

nk=input('input number of articles in system: ');
k=input('input your sub group numbers as a line matrix like [1 2 33]: ');
nu=zeros(nk,length(k));a=zeros(length(k));
for i=1:nk
    fprintf('enter number of subgroups for article number %d: \n',i)
    nu(i,:)=input('');
end
R=input('input Rk as a line matrix: ');
Q=input('input Qk as a line matrix: ');
fprintf('use table h.2 and enter all a(m,k)) \n')
for i=1:length(k)
    for j=1:length(k)
        fprintf('a(%d,%d): ',k(i),k(j))
        a(i,j)=input('');
    end
save('data.mat','nk','k','nu','R','Q','a','-v7.3')
end