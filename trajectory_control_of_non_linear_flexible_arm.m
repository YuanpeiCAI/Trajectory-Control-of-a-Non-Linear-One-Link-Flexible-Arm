clear
clc





%--------------------------parameter settings------------------------------
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
% frequency of the first mode (Hz)
f1=2.12;
% frequency of the second mode (Hz)
f2=14.3;
% spring coefficients associated with the first mode (lbf/ft)
k1=5.54;
% spring coefficients associated with the second mode (lbf/ft)
k2=198.56;
% smapling time
t=0.001;
% ending time
Te=4;
% middle time 
T=2;
% time array
i=t:t:Te;
% counting number
n=1;
j=1;


% -------------------------desired motion----------------------------------
% velocity profile
while(n*t<=Te)
    if(i(n)<=T)
        dyd(n)=(90/T)*(1-cosd(360*i(n)/T));
    else
        dyd(n)=0;
    end
    n=n+1;
end


% %---------------------------computing mode-------------------------------
% C2(2) and C2(3) represents the proportion in the first two mode
% beta(2) and beta(3) represents the parameter related to frequency
% omega(2) and omega(3) represents the frequency
M=ML/(rho*A*l);
J=Jp/(rho*A*l^3);

% beta=fsolve((1+cos(beta)*cosh(beta))-M*beta*(sin(beta)*cosh(beta)
% -cos(beta)*sinh(beta))-J*beta^3*(sin(beta)*cosh(beta)+cos(beta)*sinh(beta))
% +M*J*beta^4*(1-cos(beta)*cosh(beta)),8)
% f1=2.12;
% beta=sqrt(sqrt(rho*A*(2*pi*f1)^2*l^4/EI))
beta(2)=1.6099;
beta(3)=3.211;
j=2;
xi=t:t:1;
while(j<=3)
    n=1;
%     natural frequency
    omega(j)=sqrt(beta(j)^4*EI/(rho*A*l^4));
%     appied boundary condition and get the coefficient of the mode
    C2(j)=(cos(beta(j))+cosh(beta(j))-M*beta(j)*(sin(beta(j))-sinh(beta(j))))...
        /(sin(beta(j))-sinh(beta(j))+M*beta(j)*(cos(beta(j))-cosh(beta(j))));
%     while(n<=1000)
%         phi(j,n)=sin(beta(j)*xi(n))-sinh(beta(j)*xi(n))+C2(j)*(cos(beta(j)*xi(n))-cosh(beta(j)*xi(n)));
%         n=n+1;
%     end
    j=j+1;
end
% plot(xi,phi(2,:))
% hold on 
% plot(xi,phi(3,:))

y=rho*A*l^2*quadl(@kappa3,0,1)




  


% -------------------------essential matrix--------------------------------
% the 1st flexible mode

M(1,1)=J0+ML*l^2+I0+ML*()

















