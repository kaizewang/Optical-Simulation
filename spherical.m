%单球面像差模拟
for i=1:2
    r=100*(-1)^(i);
n=1.33;N=100;
phi=(n-1)/r;
s=2/phi;
s0=n/(phi-1/s);
h=linspace(0,abs(r)/10,N);
x1=sqrt(r.^2-h.^2);x2=(r-sign(r)*x1)+s;
alpha=atan(h./x2);beta=asin(h./abs(r));
gamma=asin(sin(sign(r)*alpha+beta)/n);
theta=beta-gamma;
d=sign(r)*(h./tan(theta)+abs(r)-x1-abs(s0));
y=d.*tan(beta-gamma);
figure(1)
plot(d,h);
xlabel('轴向位移d');
ylabel('横向伸展高度h');
grid on;
hold on;
figure(2)
plot(y,h);
xlabel('横向位移H');
ylabel('横向伸展高度h');
grid on;
hold on;
clc;
end
figure(1)
legend('Concave','Convex');
figure(2)
legend('Concave','Convex');
