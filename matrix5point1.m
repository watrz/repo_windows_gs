%A1=ones(1,65536);
%A1=4*A1;
function A=matrix5point1(n)
A2=ones(1,n^2);
A2=-1*A2;
for i=1:n
i=(n-1)*i;
A2(i)=-2;
end
for i=1:(n-1)
i=n*i;
A2(i)=0;
end
A3=ones(1,(n-1)*n);
for i=(n-1)^2:(n-1)*n
A3(i)=-2;
end
%以下是检索的方法
C=zeros(n^2);
C(n^2,n^2)=4;
for i=1:(n^2-1)
C(i,i)=4;
C(i+1,i)=A2(i);
C(i,i+1)=A2(n-i);
end
for i=1:((n-1)^2-1)
C(i+n,i)=A3(i);
C(i,i+n)=A3((n-1)^2-1);
end