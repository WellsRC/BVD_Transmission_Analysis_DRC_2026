function dxdt = SIR_ODE_Alt(t,x,R0,delta_I)

dxdt=zeros(2,1);

dxdt(1)=-delta_I.*R0.*(1-0.1.*cos(t.*pi.*0.015)).*x(1).*x(2);
dxdt(2)=delta_I.*R0.*(1-0.1.*cos(t.*pi.*0.015)).*x(1).*x(2)-delta_I.*x(2);

end