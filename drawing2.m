function drawing2(R1,R2,a,R3,d2)
%»­Í¸¾µ
r=sqrt(a*(2*R1-a));
alpha=asin(r/R1);
[q1, q2]=circle([R1-a 0],R1,pi-alpha,pi+alpha);
plot(q1,q2,'color','r','linewidth',2);
hold on;
[s1, s2]=circle([a-R2 0],R2,-alpha,alpha);
plot(s1,s2,'r','linewidth',2);
beta=asin(r/R3);
[t1,t2]=circle([a+d2-R3 0],R3,-beta,beta);
plot(t1,t2,'r','linewidth',2);
x11=q1(1);y11=q2(1);x12=q1(length(q1));y12=q2(length(q1));
x32=t1(1);y32=t2(1);x31=t1(length(t1));y31=t2(length(t1));
x=[x12 x32];y=[ y12 y32];
plot(x,y,'r','linewidth',2);
x=[x11 x31];y=[y11 y31];
plot(x,y,'r','linewidth',2);
end

