% inertia matrix

function f=M(x2,x3,n)
%beam lenth (ft)
l=4;
%beam mass (lbfs^2/ft)
m=0.0205;
% payload mass (lbfs^2/ft)
ML=0.0031;
% joint inertiar (lbfs^2*ft)
I0=1;
% beam inertia reletive to joint (lbfs^2*ft)
J0=0.109;
% payload inertia (lbfs^2*ft)
Jp=0.00852;
% beam cross area (ft^2)
A=0.000976;
% Yong's modulus*beam inertia (lbf*ft^2)
EI=28.6;
% density (lbf*s^2/ft^4)
rho=5.25;


f(1,1,n)=J0+Jp+ML*l^2+I0+ML*(phi2(1)*x2+phi3(1)*x3)^2;
f(1,2,n)=ML*l*phi2(1)+Jp*dphi2(1)+rho*A*l^2*quadl(@kappa2,0,1);
f(1,3,n)=ML*l*phi3(1)+Jp*dphi3(1)+rho*A*l^2*quadl(@kappa3,0,1);
f(2,1,n)=f(1,2,n);
f(2,2,n)=rho*A*l*quadl(@chi22,0,1)+ML*phi2(1)^2+Jp*dphi2(1)^2;
f(2,3,n)=rho*A*l*quadl(@chi23,0,1)+ML*phi2(1)*phi3(1)+Jp*dphi2(1)*dphi3(1);
f(3,1,n)=f(1,3,n);
f(3,2,n)=f(2,3,n);
f(3,3,n)=rho*A*l*quadl(@chi33,0,1)+ML*phi3(1)^2+Jp*dphi3(1)^2;
end