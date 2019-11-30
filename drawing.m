function drawing(r,R1,R2,a)
%»­Í¸¾µ
alpha=asin(r/R1);
[q1, q2]=circle([R1-a 0],R1,pi-alpha,pi+alpha);

plot(q1,q2,'color','r','linewidth',2);
[q1, q2]=circle([a-R2 0],R2,-alpha,alpha);

plot(q1,q2,'r','linewidth',2);


