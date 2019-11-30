%模拟球差
%初值设定
M=4;
nabla=zeros(1,M);
for j=1:1:M
R1=100*j;a=0.1;r=sqrt(a*(2*R1-a));
n1=1;n3=1;n2=1.55;y0=0;
N=40;%采点个数
phi=(n2-n1)/100+(n2-n3)/100;
R2=(n2-n3)/(phi-(n2-n1)/R1);
ff1=n1/phi;
u0=2*ff1;
gamma=atan((r-y0)/u0);%发散角范估计
d=2*a;
[v0,yy0]=ideal(u0-a,n1,n2,n3,R1,R2,d,y0);
theta=linspace(-0.9*gamma,0.9*gamma,N);
deltax=zeros(1,N);
deltay=zeros(1,N);
syms x1 x2
f1=R1^2-(x2^2+(x1-R1+a)^2);
f2=R2^2-x2^2-(x1+R2-a)^2;
for i=1:N
    p1=[-u0,y0];
    k=n1*[cos(theta(i)),-sin(theta(i))];%第一个k
    p=p1;
    F1=@(x)[k(1)*(x(2)-p(2))-k(2)*(x(1)-p(1));R1^2-(x(2)^2+(x(1)-R1+a)^2)];
    p2=fsolve(F1,[-a,y0]);
    clc;
    w=double(subs(gradient(f1),[x1 x2],p2));
    w=w/sqrt(w'*w);%标准化
    G1=k*w;G2=sqrt(n2^2-n1^2+G1^2);
    k=k+(G2-G1)*w';%第二个k
    p=p2;
    F2=@(x)[k(1)*(x(2)-p(2))-k(2)*(x(1)-p(1));R2^2-(x(2)^2+(x(1)+R2-a)^2)];
    p3=fsolve(F2,[a y0]);
    clc;
    w=double(subs(gradient(f2),[x1 x2],p3));
    w=w/sqrt(w'*w);%标准化
    G1=-k*w;G2=sqrt(n3^2-n2^2+G1^2);
    k=k-(G2-G1)*w';%出射的k
    x0=v0;%p4位置初始化
    p4=[x0,k(2)*(x0-p3(1))/k(1)+p3(2)];
    deltax(i)=-p3(2)*k(1)/k(2)+p3(1)-v0;%收集数据
end
nabla((j-1)/0.5+1)=R1/R2;
plot(theta,deltax);
xlabel('轴向夹角theta');
ylabel('横向偏差deltax');
hold on;
grid on;
grid minor;
end
legend('1','2','3','4');
