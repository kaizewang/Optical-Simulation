%配曲法演示
R1=3.550;R2=3.550;R3=60;a1=0.15;r=sqrt(a1*(2*R1-a1));
d2=0.2;
N=100;%采点个数
n1=1;n2=1.5166;n3=1.6256;n4=1;theta=zeros(1,N);
y0=linspace(-1,1,N);
deltax=zeros(1,N);
deltay=zeros(1,N);
syms x1 x2
f1=R1^2-(x2^2+(x1-R1+a1)^2);
f2=R2^2-x2^2-(x1+R2-a1)^2;
f3=R3^2-x2^2-(x1-a1-d2+R3)^2;
for i=1:N
    p1=[-4,y0(i)];
    k=n1*[1,0];%第一个k
    p=p1;
    F1=@(x)[k(1)*(x(2)-p(2))-k(2)*(x(1)-p(1));R1^2-(x(2)^2+(x(1)-R1+a1)^2)];
    p2=fsolve(F1,[-a1,y0(i)]);
    clc;
    w=double(subs(gradient(f1),[x1 x2],p2));
    w=w/sqrt(w'*w);%标准化
    G1=k*w;G2=sqrt(n2^2-n1^2+G1^2);
    k=k+(G2-G1)*w';%第二个k
    p=p2;
    F2=@(x)[k(1)*(x(2)-p(2))-k(2)*(x(1)-p(1));R2^2-(x(2)^2+(x(1)+R2-a1)^2)];
    p3=fsolve(F2,[a1 y0(i)]);
    clc;
    w=double(subs(gradient(f2),[x1 x2],p3));
    w=w/sqrt(w'*w);%标准化
    G1=-k*w;G2=sqrt(n3^2-n2^2+G1^2);
    k=k-(G2-G1)*w';%第三个的k
    p=p3;
    F3=@(x)[k(1)*(x(2)-p(2))-k(2)*(x(1)-p(1));R3^2-(x(2)^2+(x(1)+R3-a1-d2)^2)];
    p4=fsolve(F3,[a1+d2 y0(i)]);
    clc;
    w=double(subs(gradient(f3),[x1 x2],p4));
    w=w/sqrt(w'*w);%标准化
    G1=-k*w;G2=sqrt(n4^2-n3^2+G1^2);
    k=k-(G2-G1)*w';%出射的k
    v0=8;%p5位置初始化
    p5=[v0,k(2)*(v0-p4(1))/k(1)+p4(2)];
    deltax(i)=-p4(2)*k(1)/k(2)+p4(1);%收集数据
    if i==N/2
        x0=deltax(i);
    end
    if mod(i,20)==1
    figure(1);
    line([p1(1) p2(1) p3(1) p4(1) p5(1)],[p1(2) p2(2) p3(2) p4(2) p5(2)]);
    end
    hold on;
    grid on;
    axis equal;
end
deltax=deltax-x0;
figure(1)
 drawing2(R1,R2,a1,R3,d2);
figure(2)
plot(deltax,y0);
ylabel('入射高度H');
xlabel('轴向偏差deltax');
grid on;
grid minor;