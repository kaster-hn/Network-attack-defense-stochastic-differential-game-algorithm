function dx= stoc3_3_3_1(t,x) %Slove the node state transition equations
global p_A_H p_A_M p_A_L p_D_H p_D_L p_D_M
global k c length mu va vd sd sa c_M
global alpha Q T r_1 r_2 r_3 r_4
global c_D1 c_D2 c_D3 c_A1 c_A2 c_A3
global e_A_H1 e_A_M1 e_A_L1 e_D_H1 e_D_M1 e_D_L1
global e_A_H e_A_M e_A_L e_D_H e_D_L e_D_M

ndd1=normrnd(0,2*e_D_H1);ndd2=normrnd(0,2*e_D_M1);ndd3=normrnd(0,2*e_D_L1);
nda1=normrnd(0,2*e_A_H1);nda2=normrnd(0,2*e_A_M1);nda3=normrnd(0,2*e_A_L1);

e_D_H=sd*ndd1+e_D_H1;e_D_M=sd*ndd2+e_D_M1;e_D_L=sd*ndd3+e_D_L1;
e_A_H=sa*nda1+e_A_H1;e_A_M=sa*nda2+e_A_M1;e_A_L=sa*nda3+e_A_L1;

%% Calculate the result of the attack-defense confrontation
a=p_A_H*e_A_H+p_A_M*e_A_M+p_A_L*e_A_L;
d=p_D_H*e_D_H+p_D_L*e_D_L+p_D_M*e_D_M;
tau=a-d;
if(tau<=0)
    tau_NI=0;tau_NR=-tau ;
else
    tau_NR=0;tau_NI=tau;
end
tau_IM=tau_NI;tau_IR=tau_NR;tau_IN=tau_IR;


%%
if(t>=length*(k-1))
    c(k,1)=p_D_H;
    c(k,2)=p_D_M;
    c(k,3)=p_D_L;
    c(k,4)=p_A_H;
    c(k,5)=p_A_M;
    c(k,6)=p_A_L;
    c(k,7)=t;
    c(k,8)=tau;
    k=k+1;
end

y(4,T+1)=zeros();
y(1,T)=0;y(2,T)=0;y(3,T)=0;y(4,T)=c_M;
%% Substitute the result of the attack-defense confrontation in equations (2) and (11)¨C(14)
dx=zeros(8,1);
dx(1)=-tau_NI*alpha*pi*x(2)*x(1)/Q*x(2)/Q*va-mu*tau_NR*x(1)*(x(1)+x(3))/Q*vd+(1-mu)*tau_IN*x(2)*(x(1)+x(3))/Q*vd;
dx(2)=(tau_NI*alpha*pi*x(2)*x(1)/Q*x(2)/Q*va-tau_IM*x(2)*x(2)/Q*va-(1-mu)*tau_IN*x(2)*(x(1)+x(3))/Q*vd-mu*tau_IR*x(2)*(x(1)+x(3))/Q*vd);
dx(3)=mu*tau_NR*x(1)*(x(1)+x(3))/Q*vd+mu*tau_IR*x(2)*(x(1)+x(3))/Q*vd;
dx(4)=(tau_IM*x(2)*x(2)/Q*va);
for t=T:-1:2
    y(1,t-1)=y(1,t)-(r_1+y(1,t)-y(2,t))*tau_NI*alpha*pi*x(2)*x(2)/Q/Q-mu*tau_NR*(y(1,t)-y(3,t)-r_2)*(2*x(1)+x(3))/Q...
        +tau_IR*x(2)*(r_2-y(2,t)+y(3,t))*mu/Q+tau_IN*x(2)*(r_4+y(1,t)-y(2,t))*(1-mu)/Q...
        -1/2*(p_D_H*e_D_H)^2*c_D1-1/2*(p_D_M*e_D_M)^2*c_D2-1/2*(p_D_L*e_D_L)^2*c_D3;
    y(2,t-1)=y(2,t)+tau_IR*(r_2-y(2,t)+y(3,t))*mu*(x(1)+x(3))/Q+tau_IN*(r_4+y(1,t)-y(2,t))*(1-mu)*(x(1)+x(3))/Q...
        +2*tau_NI*(-r_1-y(1,t)+y(2,t))*alpha*x(2)*x(1)/Q/Q+2*tau_IM*(-r_3-y(2,t)-y(4,t))*x(2)/Q...
        +1/2*(p_A_H*e_A_H)^2*c_A1+1/2*(p_A_M*e_A_M)^2*c_A2+1/2*(p_A_L*e_A_L)^2*c_A3;
    y(3,t-1)=y(3,t)+tau_NR*(r_2-y(1,t)+y(3,t))*mu*x(1)/Q+tau_IR*(r_2-y(2,t)+y(3,t))*mu*x(2)/Q...
        +tau_IN*(r_4+y(1,t)-y(2,t))*(1-mu)*x(2)/Q...
        -1/2*(p_D_H*e_D_H)^2*c_D1-1/2*(p_D_M*e_D_M)^2*c_D2-1/2*(p_D_L*e_D_L)^2*c_D3;
    y(4,t-1)=y(4,t)-0;
end

pd(1)=p_D_H;pd(2)=p_D_M;pd(3)=p_D_L;
pa(1)=p_A_H;pa(2)=p_A_M;pa(3)=p_A_L;

Aeqd=[1 1 1];beqd=[1];
lbd=[0;0;0];ubd=[1;1;1];
Aeqa=[1 1 1];beqa=[1];
lba=[0;0;0];uba=[1;1;1];

%% Calculate the optimal strategies of the defender and attacker
qd=fmincon(@(pd)ssd3_3_3_1(x,tau,pd,pa),pd,[],[],Aeqd,beqd,lbd,ubd);
qa=fmincon(@(pa)ssa3_3_3_1(x,tau,pd,pa),pa,[],[],Aeqa,beqa,lba,uba);


p_D_H=qd(1);p_D_M=qd(2);p_D_L=qd(3);
p_A_H=qa(1);p_A_M=qa(2);p_A_L=qa(3);


t;
k;
end

