function y=main_gs(he,num,N_2d)%这个是主程序
format long
%num=8;%晶格大小
num1=1000000;%迭代次数
e=1.6*10^(-19);%电子电量
m_e=9.1*10^(-31);%电子质量
u0=4*pi*10^(-7);%真空磁导率
Ms=8.0*10^5;%磁化强度大小
gamma=-e/(2*m_e);%
A=1.3*10^(-11);%交换常数
Ku=5.0*100;%各向异性常数
alpha=0.2;%吉尔伯特耗散常数
delta_t=2*10^(-12);%时间间隔
delta_x=4*10^(-9);%空间间隔
L=1*10^(-6);%区域面积
H0=400;%外加最大磁场，单位（/Oe）
epsilon=A/(u0*Ms^2*L^2);
Q=Ku/(u0*Ms^2);
%=====================================
fivematrix=matrix5point(num);%此处于润泽负责，由他的子函数调入,但最终作为了备用
h_s=zeros(3,num^2);%退磁场，由子函数调入，但是因为只要某一个时间点的所有格点的磁矩都求出来后，那么下一步计算的所有退磁场就已定，
%所以我们要在每一次整层磁矩计算完成后，计算一遍退磁场，所以退磁场应该在循环内计算
h_e=repmat(he',[num^2,1]);%外场
%=====================================
%下面定义包含有磁矩初始条件的迭代储存矩阵
m1=[ones(num^2,1),ones(num^2,num1)];
m2=[zeros(num^2,1),ones(num^2,num1)];
m3=[zeros(num^2,1),ones(num^2,num1)];
%=======================================
%%此处矩阵求逆需要检查
g0=((eye(num^2)-epsilon*delta_t*fivematrix))^(-1);%计算第一个逆因子
g00=((eye(num^2)-alpha*epsilon*delta_t*fivematrix))^(-1);%计算第二个逆因子
g1=0;
g2=0;
g3=0;
g1_star=0;
g2_star=0;
g3_star=0;
m1_star1=0;
m2_star1=0;
m3_star1=0;
m1_star2=0;
m2_star2=0;
m3_star2=0;
m1_star=[];
m2_star=[];
m3_star=[];
f=ones(num^2,1);
%=======================================
%=======================================
for n=1:num1
    h_s=Demagnet(num,m1(:,n),m2(:,n),m3(:,n),N_2d);
    
    %=============第一步=================================================
    %此处我们有可能是调入退磁场程序，每循环一遍需要计算一下每个格点的退磁场大小
    for k=1:num^2
        f=-Q*[zeros(num^2,1),m2(:,n),m3(:,n)]+h_s'+h_e;%众向排列
        g1=g0(k,:)*(m1(:,n)+delta_t*f(:,1));
        g2=g0(k,:)*(m2(:,n)+delta_t*f(:,2));
        g3=g0(k,:)*(m3(:,n)+delta_t*f(:,3));
        m1_star1=m1(k,n)+(g2*m3(k,n)-g3*m2(k,n));
        m1_star=[m1_star;m1_star1];
    end
    
    for k=1:num^2
        f=-Q*[zeros(num^2,1),m2(:,n),m3(:,n)]+h_s'+h_e;%众向排列
        g3=g0(k,:)*(m3(:,n)+delta_t*f(:,3));
        g1_star=g0(k,:)*(m1_star+delta_t*f(:,1));
        m2_star1=m2(k,n)+(g3*m1_star(k)-g1_star*m3(k,n));
        m2_star=[m2_star;m2_star1];
    end
  
    for k=1:num^2
        f=-Q*[zeros(num^2,1),m2(:,n),m3(:,n)]+h_s'+h_e;%众向排列
        g1_star=g0(k,:)*(m1_star+delta_t*f(:,1));
        g2_star=g0(k,:)*(m2_star+delta_t*f(:,2));
        m3_star1=m3(k,n)+(g1_star*m2_star(k)-g2_star*m1_star(k));
        m3_star=[m3_star;m3_star1];
    end
    
    for k=1:num^2
        %%=============第二步================================================  
        
        f_star=-Q*[zeros(num^2,1),m2_star,m3_star]+h_s'+h_e;%此处的退磁场在理论上与第一步不一样，暂时编程是用的与上面一样的数值，怎样变化还未处理
        m1_star2=g00(k,:)*(m1_star+alpha*delta_t*f_star(:,1));
        m2_star2=g00(k,:)*(m2_star+alpha*delta_t*f_star(:,2));   
        m3_star2=g00(k,:)*(m3_star+alpha*delta_t*f_star(:,3));
        
        %%=============第三步=================================================
        
        model=1/sqrt(m1_star2^2+m2_star2^2+m3_star2^2);
        m1(k,n+1)=model*m1_star2;
        m2(k,n+1)=model*m2_star2;
        m3(k,n+1)=model*m3_star2;
    end
    m1_star=[];
    m2_star=[];
    m3_star=[];
end
abs(mean(real(m1(:,10)))-mean(real(m1(:,11))))/mean(real(m1(:,10)))
y=real([m1;m2;m3])

