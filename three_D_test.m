%3D situation
%模拟球差
%初值设定
R1=100;R2=100;a=20;r=sqrt(a*(2*R1-a));
n1=1;n3=1;n2=1.55;y0=20;
N=100;%采点个数
phi=(n2-n1)/R1+(n2-n3)/R2;
ff1=n1/phi;
u0=2*ff1;
d=2*a;
c=2*sqrt(R1^2-y0^2)+2*a-2*R1;
gamma=atan(c/u0);%发散角范估计
[v0,yy0]=ideal(u0-a,n1,n2,n3,R1,R2,d,y0);
alpha=linspace(-0.8*gamma,0.8*gamma,N);
deltax=zeros(1,N);
deltay=zeros(1,N);
syms x1 x2 x3
f1=R1^2-(x2^2+(x1-R1+a)^2+x3^2);
f2=R2^2-x2^2-(x1+R2-a)^2-x3^2;
zeta=pi/2;
for i=1:N
    p1=[-u0,y0,0];
    k=n1*[cos(alpha(i))*sin(zeta),cos(zeta),sin(alpha(i))*sin(zeta)];%第一个k
    p=p1;         
    F1=@(x)[k(1)*(x(3)-p(3))-k(3)*(x(1)-p(1));R1^2-(x(2)^2+(x(1)-R1+a)^2+x(3)^2);...
        k(1)*(x(2)-p(2))-k(2)*(x(1)-p(1))];
    p2=fsolve(F1,[-a,y0,0]); 
    w=double(subs(gradient(f1),[x1 x2 x3],p2));
    w=w/sqrt(w'*w);%标准化
    G1=abs(k*w);G2=sqrt(n2^2-n1^2+G1^2);
    k=k+(G2-G1)*w';%第二个k
    p=p2;
    F2=@(x)[k(1)*(x(2)-p(2))-k(2)*(x(1)-p(1));...
        k(1)*(x(3)-p(3))-k(3)*(x(1)-p(1));...
        R2^2-(x(2)^2+(x(1)+R2-a)^2+x(3)^2)];
    p3=fsolve(F2,[a y0 0]);
    w=double(subs(gradient(f2),[x1 x2 x3],p3));
    w=w/sqrt(w'*w);%标准化
    G1=abs(k*w);G2=sqrt(n3^2-n2^2+G1^2);
    k=k+(G2-G1)*(-w');%出射的k
    x0=v0;%p4位置初始化
    p4=[x0,k(2)*(x0-p3(1))/k(1)+p3(2),...
        k(3)*(x0-p3(1))/k(1)+p3(3)];
    deltaz(i)=p4(3);
    if i==N/2
        z0=deltaz(i);
    end
    if mod(i,20)==1
    figure(1);
    line([p1(1) p2(1) p3(1) p4(1)],[p1(2) p2(2) p3(2) p4(2)],...
        [p1(3) p2(3) p3(3) p4(3)]);
    end
    hold on;
    grid on;
    axis equal;
end
deltaz=deltaz-z0;
figure(1)
 drawing3D(R1,a);
 figure(2)
 plot(alpha,deltaz);
