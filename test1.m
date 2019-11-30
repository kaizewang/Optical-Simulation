%ģ�����
%��ֵ�趨
R1=100;R2=100;a=5;r=sqrt(a*(2*R1-a));
n1=1;n3=1;n2=1.55;y0=0;
N=100;%�ɵ����
phi=(n2-n1)/R1+(n2-n3)/R2;
ff1=n1/phi;
u0=2*ff1;
gamma=atan((r-y0)/u0);%��ɢ�Ƿ�����
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
    k=n1*[cos(theta(i)),-sin(theta(i))];%��һ��k
    p=p1;
    F1=@(x)[k(1)*(x(2)-p(2))-k(2)*(x(1)-p(1));R1^2-(x(2)^2+(x(1)-R1+a)^2)];
    p2=fsolve(F1,[-a,y0]);
    clc;
    w=double(subs(gradient(f1),[x1 x2],p2));
    w=w/sqrt(w'*w);%��׼��
    G1=k*w;G2=sqrt(n2^2-n1^2+G1^2);
    k=k+(G2-G1)*w';%�ڶ���k
    p=p2;
    F2=@(x)[k(1)*(x(2)-p(2))-k(2)*(x(1)-p(1));R2^2-(x(2)^2+(x(1)+R2-a)^2)];
    p3=fsolve(F2,[a y0]);
    clc;
    w=double(subs(gradient(f2),[x1 x2],p3));
    w=w/sqrt(w'*w);%��׼��
    G1=-k*w;G2=sqrt(n3^2-n2^2+G1^2);
    k=k-(G2-G1)*w';%�����k
    x0=v0;%p4λ�ó�ʼ��
    p4=[x0,k(2)*(x0-p3(1))/k(1)+p3(2)];
    deltay(i)=p4(2)-yy0;
    deltax(i)=-p3(2)*k(1)/k(2)+p3(1)-v0;%�ռ�����
    if mod(i,10)==1
    figure(1);
    line([p1(1) p2(1) p3(1) p4(1)],[p1(2) p2(2) p3(2) p4(2)]);
    end
    hold on;
    grid on;
    axis equal;
end
figure(1)
 drawing(r,R1,R2,a);
figure(2)
subplot(2,1,1);
plot(theta,deltax);
xlabel('����н�theta');
ylabel('����ƫ��deltax');
grid on;
grid minor;
figure(2)
subplot(2,1,2);
plot(theta,deltay);
xlabel('����н�theta');
ylabel('����ƫ��deltax');
grid on;
grid minor;