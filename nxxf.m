% nxxf is a transcendental equation

function eq=nxxf(beta)
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
M=ML/(rho*A*l);
J=Jp/(rho*A*l^3);
eq=(1+cos(beta).*cosh(beta))-M.*beta.*(sin(beta).*cosh(beta)-cos(beta).*sinh(beta))-J.*beta.^3.*(sin(beta).*cosh(beta)+cos(beta).*sinh(beta))+M.*J.*beta.^4.*(1-cos(beta).*cosh(beta));
end





