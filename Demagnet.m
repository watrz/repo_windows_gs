function H_s=Demagnet(num,mm1,mm2,mm3,N_2d) 
%此处输入的mm1，mm2,mm3为 (num*num)*1的大小的矩阵，表示第n时刻每个磁矩的分量，所以这个数据对应于主程序中存储磁矩的数据矩阵的某一列，表示每个格点处的磁矩大小
%当然还需要调入函数demag
%输出的是一个3*（num*num）的矩阵，表示在某个时刻所有格点上的退磁场；
k=-4*pi;%这是式子前面的系数；
%num=256;%格点数
%N_2d=demag(num);%元胞数组，里面是每个点对另外点的退磁矩阵；
N_beforefft=zeros(9,num*num);%此处是FFT之前的数据，是将第i个格点与第i个格点本身的退磁矩阵算上了，当然令其为零,
N_afterfft=zeros(9,num*num);%这是fft之后的数据存储
mm_before=[mm1,mm2,mm3];%输入的磁矩数据整合
mm_after=zeros(num^2,3);%磁矩FFT之后的数据
mcross=zeros(3,num^2);%K空间退磁矩阵和磁矩的点乘后的数据
H_s1=zeros(num^2,3);%这是最后求得的退磁场
%下面是将元胞数组的数据提取出来排成(9*num^2)*num^2的大小，方便进行傅里叶变换
%for i=1:num^2
    for j=1:num^2
        for k=1:9
            k1=rem(k,3);%求余数
            temp=~k1;%避免整除3时导致下面的列指标为0；
            k2=(k-k1)/3;%求商
            N_beforefft(k,j)=N_2d{1,j}(k2+(k1~=0),k1+temp*3);
        end
    end
%end
%=============先对上述矩阵进行FFT，我们是分行进行的========================================
%对提取数据后的退磁矩阵进行快速傅里叶变换
for k=1:9
    N_afterfft(k,:)=fft(N_beforefft(k,:));%这里的数据是横着排列的
end
%对磁矩数据进行傅里叶变换
for k=1:3
    mm_after(:,k)=fft(mm_before(:,k));%这是竖着排列的三列
end
%变到K空间后的点乘，对于每个格点输出num*num个3*1的K空间的磁矩
for k=1:3
    mcross(k,:) = N_afterfft((k-1)*3+1,1).*(mm_after(:,1)')+N_afterfft((k-1)*3+2,2).*(mm_after(:,2)')+N_afterfft((k-1)*3+3,3).*(mm_after(:,3)');
end

H_s1=[ifft(mcross(1,:));ifft(mcross(2,:));ifft(mcross(3,:))];
H_s=H_s1;

