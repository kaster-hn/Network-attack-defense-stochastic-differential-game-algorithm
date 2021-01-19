function ha=ssa3_3_3_1(x,tau,pd,pa)%Solve the optimal attack strategy
global alpha Q r_1 r_2 r_3 r_4 mu vd va
global e_A_H e_A_M e_A_L e_D_H e_D_L e_D_M
global c_A1 c_A2 c_A3 

if(tau<=0)
 ha=((pd(1)*e_D_H+pd(2)*e_D_M+pd(3)*e_D_L-pa(1)*e_A_H-pa(2)*e_A_M-pa(3)*e_A_L)...
     *(r_2*mu*x(1)*vd*(x(1)+x(3))/Q+mu*r_2*x(2)*vd*(x(1)+x(3))/Q+(1-mu)*r_4*x(2)*vd*(x(1)+x(3))/Q)...
     +1/2*(pa(1)*e_A_H)^2*c_A1*x(2)+1/2*(pa(2)*e_A_M)^2*c_A2*x(2)+1/2*(pa(3)*e_A_L)^2*c_A3*x(2));
else 
 ha=((pa(1)*e_A_H+pa(2)*e_A_M+pa(3)*e_A_L-pd(1)*e_D_H-pd(2)*e_D_M-pd(3)*e_D_L)...
     *(((-r_1)*alpha*pi*x(2)*va*x(1)*x(2))/(Q*Q)+x(2)*(-r_3)*va*x(2)/Q)...
     +1/2*(pa(1)*e_A_H)^2*c_A1*x(2)+1/2*(pa(2)*e_A_M)^2*c_A2*x(2)+1/2*(pa(3)*e_A_L)^2*c_A3*x(2));
end

end
