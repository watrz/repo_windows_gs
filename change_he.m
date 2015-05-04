tic;
num=8;
N_2d=demag(num);
he=zeros(150,3);
temp=[102:2:400];
he=[temp*cos(pi/180);zeros(1,150);temp*sin(pi/180)];
for k=1:1
    hee=he(:,k);
    yy=main_gs(hee,num,N_2d);
end
toc;