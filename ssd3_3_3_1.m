function hd=ssd3_3_3_1(x,tau,pd,pa)%Solve the optimal defense strategy
global alpha Q r_1 r_2 r_3 r_4 mu
global e_A_H e_A_M e_A_L e_D_H e_D_L e_D_M
global c_D1 c_D2 c_D3 vd va

%% Calculate the defender's payoff
if(tau<=0)
    hd=-((pd(1)*e_D_H+pd(2)*e_D_M+pd(3)*e_D_L-pa(1)*e_A_H-pa(2)*e_A_M-pa(3)*e_A_L)...
        *(r_2*mu*x(1)*vd*(x(1)+x(3))/Q+mu*r_2*x(2)*vd*(x(1)+x(3))/Q+(1-mu)*r_4*x(2)*vd*(x(1)+x(3))/Q)...
        -1/2*(pd(1)*e_D_H)^2*c_D1*(x(1)+x(3))-1/2*(pd(2)*e_D_M)^2*c_D2*(x(1)+x(3))-1/2*(pd(3)*e_D_L)^2*c_D3*(x(1)+x(3)));
else
    hd=-((pa(1)*e_A_H+pa(2)*e_A_M+pa(3)*e_A_L-pd(1)*e_D_H-pd(2)*e_D_M-pd(3)*e_D_L)...
        *(((-r_1)*alpha*pi*x(2)*x(1)*va*x(2))/(Q*Q)+(-r_3)*x(2)*va*x(2)/Q)...
        -1/2*(pd(1)*e_D_H)^2*c_D1*(x(1)+x(3))-1/2*(pd(2)*e_D_M)^2*c_D2*(x(1)+x(3))-1/2*(pd(3)*e_D_L)^2*c_D3*(x(1)+x(3)));
end

end
