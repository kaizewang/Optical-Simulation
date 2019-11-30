function [v,y]=ideal(u0,n1,n2,n3,R1,R2,d,y0)
phi1=(n2-n1)/R1;
phi2=(n2-n3)/R2;
u11=n2/(phi1-n1/u0);
u1=d-u11;
u22=n3/(phi2-n2/u1);
v=u22+d/2;
y=y0*n1/n2*u11/u0*n2/n3*u22/u1;
end