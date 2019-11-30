function [x y]=circle(a,R,alpha,beta)
phi=linspace(alpha,beta,1000);
x=a(1)+R*cos(phi);
y=a(2)+R*sin(phi);
end