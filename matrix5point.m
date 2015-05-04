function A=matrix5point(n)
%=========以下是先定义一些矩阵块===========================
%n=3;
A1=-1*eye(n);
A2=-2*eye(n);
A3=zeros(n);
%==========实现两边上下的块矩阵======================
%A1(n,1)=-1;
%A2(n,1)=-1;
%A11=A1';
%A22=A2';
%==========实现当中的矩阵================================================
triabove =[-2,-1*ones(1,n-2)];
tribelow =[-1*ones(1,n-2),-2];
A31=diag(triabove,1) + diag(tribelow, -1) + 4*eye(n);
%========
%A31=diag(-1*eye(n-1), 1) + diag(-1*eye(n-1), -1) + 4*eye(n);
%==========实现整个矩阵==================================================
%利用元胞数组，两次转换即可
A_cell=cell(n);
%============解决第一行=======
% for i=1:n
%     if i==1
%         A_cell(1,i)={A31};
%     elseif i==2
%         A_cell(1,i)={A2};
%     else
%         A_cell(1,i)={zeros(n)};
%     end 
% end
% %==========解决最后一行==========================
% for i=1:n
%     if i==n-1
%         A_cell(n,i)={A2};
%     elseif i==n
%         A_cell(n,i)={A31};
%     else
%         A_cell(n,i)={zeros(n)};
%     end
% end
% %===========解决当中的（n-2）行==========================
% for i=2:n-1
%     for j=1:n
%         if j==i
%             A_cell(i,j)={A31};
%         elseif j==(i-1)||(i+1)
%             A_cell(i,j)={A1};
%         else
%             A_cell(i,j)={zeros(n)}
%         end
%     end
% end
for i=1:n
    A_cell(i,i)={A31};
    if i==1
        A_cell(i,i+1) = {A2};
    elseif i==n
        A_cell(i,i-1) = {A2};
    else
        A_cell(i,i-1) = {A1};
        A_cell(i,i+1) = {A1};
    end
end

for i=1:n
    for j=1:n
        if isempty(A_cell{i,j})
            A_cell{i,j} = zeros(n);
        end
    end
end
A=cell2mat(A_cell);
% A=cell2mat(A_cell);