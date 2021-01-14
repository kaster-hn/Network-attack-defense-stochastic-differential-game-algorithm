clc
clear all
% main function
global p_A_H p_A_M p_A_L p_D_H p_D_M p_D_L
global e_A_H1 e_A_M1 e_A_L1 e_D_H1 e_D_M1 e_D_L1 
global k c length mu va vd
global c_N c_I c_R c_M
global alpha Q r_1 r_2 r_3 r_4 
global c_D1 c_D2 c_D3 c_A1 c_A2 c_A3 
global sa sd ns

length=4;
for j=1:1:3
    if(j==1)
        ns=0;sd=ns;sa=ns;length=8;
    else if(j==2)
            ns=0.01;sd=ns;sa=ns+0.02;
        else
            ns=0.06;sd=ns;sa=ns+0.02;
        end
    end
    
t1=clock;

alpha=10;  
r_1=4;r_2=7;r_3=10;r_4=4;
e_A_H1=0.82;e_A_M1=0.45;e_A_L1=0.3;e_D_H1=0.76;e_D_M1=0.5;e_D_L1=0.267;
c_A1=4.4;c_A2=2.4;c_A3=1.5;c_D1=4;c_D2=2.7;c_D3=1.3;
c_N=2;c_I=2;c_R=1;c_M=3;c_MT=6;va=0.7;vd=0.7;mu=0.8;k=1;
p_A_H=1;p_A_M=0;p_A_L=0;
p_D_H=1;p_D_M=0;p_D_L=0;
Q=1000;I=10;
x0=[Q-I,I,0,0,0,0,0,0];

td=100*length;
ts=0:length:td;

[t,x]=ode45(@(t,x)stoc3_3_3_1(t,x),ts,x0);
figure('Name',' ')
plot(t,x(:,1),'b-',t,x(:,2),'r-',t,x(:,3),'g-',t,x(:,4),'k-')
legend('N','I','R','M')
xlabel('t')
ylabel('Number of nodes')

xx=c(:,7);
figure('Name',' ')
plot(xx,c(:,1),'b',xx,c(:,2),'g',xx,c(:,4),'r',xx,c(:,5),'m',xx,c(:,8),'k')
legend('p_D^H','p_D^M','p_A^H','p_A^M','\tau')
xlabel('t')
ylabel('Strategies')

mx=x(:,1:4);
mc(:,1:2)=c(:,1:2);
mc(:,3:4)=c(:,4:5);
mc(:,5)=c(:,8);

t2=clock;
t=etime(t2,t1)

end
