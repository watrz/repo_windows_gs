function y=main_gs(he,num,N_2d)%�����������
format long
%num=8;%�����С
num1=1000000;%��������
e=1.6*10^(-19);%���ӵ���
m_e=9.1*10^(-31);%��������
u0=4*pi*10^(-7);%��մŵ���
Ms=8.0*10^5;%�Ż�ǿ�ȴ�С
gamma=-e/(2*m_e);%
A=1.3*10^(-11);%��������
Ku=5.0*100;%�������Գ���
alpha=0.2;%�������غ�ɢ����
delta_t=2*10^(-12);%ʱ����
delta_x=4*10^(-9);%�ռ���
L=1*10^(-6);%�������
H0=400;%������ų�����λ��/Oe��
epsilon=A/(u0*Ms^2*L^2);
Q=Ku/(u0*Ms^2);
%=====================================
fivematrix=matrix5point(num);%�˴����������������Ӻ�������,��������Ϊ�˱���
h_s=zeros(3,num^2);%�˴ų������Ӻ������룬������ΪֻҪĳһ��ʱ�������и��Ĵžض����������ô��һ������������˴ų����Ѷ���
%��������Ҫ��ÿһ������žؼ�����ɺ󣬼���һ���˴ų��������˴ų�Ӧ����ѭ���ڼ���
h_e=repmat(he',[num^2,1]);%�ⳡ
%=====================================
%���涨������дžس�ʼ�����ĵ����������
m1=[ones(num^2,1),ones(num^2,num1)];
m2=[zeros(num^2,1),ones(num^2,num1)];
m3=[zeros(num^2,1),ones(num^2,num1)];
%=======================================
%%�˴�����������Ҫ���
g0=((eye(num^2)-epsilon*delta_t*fivematrix))^(-1);%�����һ��������
g00=((eye(num^2)-alpha*epsilon*delta_t*fivematrix))^(-1);%����ڶ���������
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
    
    %=============��һ��=================================================
    %�˴������п����ǵ����˴ų�����ÿѭ��һ����Ҫ����һ��ÿ�������˴ų���С
    for k=1:num^2
        f=-Q*[zeros(num^2,1),m2(:,n),m3(:,n)]+h_s'+h_e;%��������
        g1=g0(k,:)*(m1(:,n)+delta_t*f(:,1));
        g2=g0(k,:)*(m2(:,n)+delta_t*f(:,2));
        g3=g0(k,:)*(m3(:,n)+delta_t*f(:,3));
        m1_star1=m1(k,n)+(g2*m3(k,n)-g3*m2(k,n));
        m1_star=[m1_star;m1_star1];
    end
    
    for k=1:num^2
        f=-Q*[zeros(num^2,1),m2(:,n),m3(:,n)]+h_s'+h_e;%��������
        g3=g0(k,:)*(m3(:,n)+delta_t*f(:,3));
        g1_star=g0(k,:)*(m1_star+delta_t*f(:,1));
        m2_star1=m2(k,n)+(g3*m1_star(k)-g1_star*m3(k,n));
        m2_star=[m2_star;m2_star1];
    end
  
    for k=1:num^2
        f=-Q*[zeros(num^2,1),m2(:,n),m3(:,n)]+h_s'+h_e;%��������
        g1_star=g0(k,:)*(m1_star+delta_t*f(:,1));
        g2_star=g0(k,:)*(m2_star+delta_t*f(:,2));
        m3_star1=m3(k,n)+(g1_star*m2_star(k)-g2_star*m1_star(k));
        m3_star=[m3_star;m3_star1];
    end
    
    for k=1:num^2
        %%=============�ڶ���================================================  
        
        f_star=-Q*[zeros(num^2,1),m2_star,m3_star]+h_s'+h_e;%�˴����˴ų������������һ����һ������ʱ������õ�������һ������ֵ�������仯��δ����
        m1_star2=g00(k,:)*(m1_star+alpha*delta_t*f_star(:,1));
        m2_star2=g00(k,:)*(m2_star+alpha*delta_t*f_star(:,2));   
        m3_star2=g00(k,:)*(m3_star+alpha*delta_t*f_star(:,3));
        
        %%=============������=================================================
        
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

