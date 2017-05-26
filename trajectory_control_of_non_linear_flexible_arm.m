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
syms n
n=1;
j=1;
% sign of open loop control
open=1;


% -------------------------desired motion----------------------------------
% velocity profile
% dyd is the desired velocity of actuator
% yd is the desired displacement of actuator
% ddyd is the desired acceleration of actuator
while(n*t<=Te)
    if(i(n)<=T)
        yd(n)=(90/T)*(i(n)-1/pi*sind(360*i(n)/T));
        dyd(n)=(90/T)*(1-cosd(360*i(n)/T));
        ddyd(n)=(90/T)*pi*sind(360*i(n)/T);
    else
        yd(n)=yd(n-1);
        dyd(n)=0;
        ddyd(n)=0;
    end
    n=n+1;
end
plot(i,dyd)



% %---------------------------computing mode-------------------------------
% C2(2) and C2(3) represents the proportion in the first two mode
% beta(2) and beta(3) represents the parameter related to frequency
% omega(2) and omega(3) represents the frequency
M=ML/(rho*A*l);
J=Jp/(rho*A*l^3);
% beta=fsolve(@nxxf,3)
% f=sqrt(beta^4*EI/(rho*A*l^4))/(2*pi)
beta(2)=1.6099;
beta(3)=3.2114;
j=2;
% xi=t:t:1;
while(j<=3)
%     n=1;
%     natural frequency
    omega(j)=sqrt(beta(j)^4*EI/(rho*A*l^4));
%     appied boundary condition and get the coefficient of the mode
    C2(j)=(cos(beta(j))+cosh(beta(j))-M*beta(j)*(sin(beta(j))-sinh(beta(j))))...
        /(sin(beta(j))-sinh(beta(j))+M*beta(j)*(cos(beta(j))-cosh(beta(j))));
%     C2(j)=(J*beta(j)^3*(cos(beta(j))-cosh(beta(j)))+sin(beta(j))+sinh(beta(j)))...
%         /(J*beta(j)^3*(sin(beta(j))+sinh(beta(j)))-(cos(beta(j))+cosh(beta(j))))
%     while(n<=1000)
%         phi(j,n)=sin(beta(j)*xi(n))-sinh(beta(j)*xi(n))+C2(j)*(cos(beta(j)*xi(n))-cosh(beta(j)*xi(n)));
%         n=n+1;
%     end
    j=j+1;
end
% plot(xi,phi(2,:))
% hold on 
% plot(xi,phi(3,:))







% -------------------------essential matrix--------------------------------
% the 1st flexible mode
% M matrix is the inertia matrix
% K matrix is the equivalent-spring constant matrix
% F matrix is the damping matrix
% x(1) is theta
% x(2) and x(3) are delta(1) and delta(2) respectively
% x(4) is the derivitive of theta 
% x(5) and x(6) are derivitive of delta(1) and delta(2) in terms of time
syms xi;
dphi2(xi)=diff(phi2(xi),xi,1);
dphi3(xi)=diff(phi3(xi),xi,1);
ddphi2(xi)=diff(dphi2(xi),xi,1);
ddphi3(xi)=diff(dphi3(xi),xi,1);
% inertia matrix
M(1,2)=ML*l*phi2(1)+Jp*dphi2(1)+rho*A*l^2*quadl(@kappa2,0,1);
M(1,3)=ML*l*phi3(1)+Jp*dphi3(1)+rho*A*l^2*quadl(@kappa3,0,1);
M(2,1)=M(1,2);
M(2,2)=rho*A*l*quadl(@chi22,0,1)+ML*phi2(1)^2+Jp*dphi2(1)^2;
M(2,3)=rho*A*l*quadl(@chi23,0,1)+ML*phi2(1)*phi3(1)+Jp*dphi2(1)*dphi3(1);
M(3,1)=M(1,3);
M(3,2)=M(2,3);
M(3,3)=rho*A*l*quadl(@chi33,0,1)+ML*phi3(1)^2+Jp*dphi3(1)^2;
% stiffness matrix
% % when I calculated integral of zeta, first I run the program to calculated 
% % ddphi and then I paste the expression in the zeta profile
% syms a;
% K(1,1)=int(ddphi2(xi)^2,0,1);
K(1,1)=quadl(@zeta2,0,1);
% K(2,2)=int(ddphi3(xi)^2,0,1);
K(2,2)=quadl(@zeta3,0,1);
K(1,2)=0;
K(2,1)=0;


% damping matrix
F(1,1)=0;
F(2,1)=0;



% ------------------------desired flexible mode----------------------------
xd(1,1)=0;
xd(2,1)=0;
xd(3,1)=0;
xd(4,1)=0;
n=1;
while(n*t<=Te)
    % nonlinear term
    n2(1,1)=-ML*dyd(n)^2*(phi2(1)*phi2(1)*xd(1,n)+phi2(1)*phi3(1)*xd(2,n));
    n2(2,1)=-ML*dyd(n)^2*(phi2(1)*phi3(1)*xd(1,n)+phi3(1)*phi3(1)*xd(2,n)); 
    %D is acceleration
    D(:,n)=[M(2,2) M(2,3); M(3,2) M(3,3)]^-1*([0; 0]-[K(1,1)*xd(1,n); K(2,2)*xd(2,n)]...
        -[F(1,1)*xd(3,n); F(2,1)*xd(4,n)]-[n2(1,1); n2(2,1)]-[M(2,1)*ddyd(n); M(3,1)*ddyd(n)]);
    %interation
    dx1d(n)=xd(3,n);
    dx2d(n)=xd(4,n);
    dx3d(n)=D(1,n);
    dx4d(n)=D(2,n);
    xd(1,n+1)=xd(1,n)+dx1d(n)*t;
    xd(2,n+1)=xd(2,n)+dx2d(n)*t;
    xd(3,n+1)=xd(3,n)+dx3d(n)*t;
    xd(4,n+1)=xd(4,n)+dx4d(n)*t;
    n=n+1;
end
% i=t:t:Te+t;
% plot(i,xd(1,:))
% hold on;
plot(i,D(1,:));
hold on;
plot(i,xd(2,:))


%-------------------joint-based inversion control--------------------------
if(open==1)
    n=1;
    while(n*t<=Te)
        % inertia matrix
        M(1,1)=J0+Jp+ML*l^2+I0+ML*(phi2(1)*xd(1,n)+phi3(1)*xd(2,n))^2;
        % nonlinear term
        n1=2*ML*dyd(n)*(phi2(1)*xd(1,n)+phi3(1)*xd(2,n))*(phi2(1)*xd(3,n)+phi3(1)*xd(4,n));
        n2(1,1)=-ML*dyd(n)^2*(phi2(1)*phi2(1)*xd(1,n)+phi2(1)*phi3(1)*xd(2,n));
        n2(2,1)=-ML*dyd(n)^2*(phi2(1)*phi3(1)*xd(1,n)+phi3(1)*phi3(1)*xd(2,n)); 
        B=M^-1;
        u(n)=n1+(M(1,1)-[M(1,2) M(1,3)]*[M(2,2) M(2,3); M(3,2) M(3,3)]^-1....
            *[M(2,1); M(3,1)])*(ddyd(n)+[B(1,2) B(1,3)]*([n2(1,1); n2(2,1)]...
            +[K(1,1)*xd(1,n); K(2,2)*xd(2,n)]+[F(1,1)*xd(3,n); F(2,1)*xd(4,n)]));
        n=n+1;
    end
end


% discritization
% initialization
x(1,1)=0;
x(2,1)=0;
x(3,1)=0;
x(4,1)=0;
x(5,1)=0;
x(6,1)=0;
n=1;
while(n*t<=Te)
    % inertia matrix
    M(1,1)=J0+Jp+ML*l^2+I0+ML*(phi2(1)*x(2,n)+phi3(1)*x(3,n))^2;
    % nonlinear term
    n1=2*ML*x(4)*(phi2(1)*x(2,n)+phi3(1)*x(3,n))*(phi2(1)*x(5,n)+phi3(1)*x(6,n));
    n2(1,1)=-ML*x(4)^2*(phi2(1)*phi2(1)*x(2,n)+phi2(1)*phi3(1)*x(3,n));
    n2(2,1)=-ML*x(4)^2*(phi2(1)*phi3(1)*x(2,n)+phi3(1)*phi3(1)*x(3,n)); 
    %D is acceleration
    D=M^-1*([u(n); 0; 0]-[0; K(1,1)*x(2,n); K(2,2)*x(3,n)]...
        -[0; F(1,1)*x(5,n); F(2,1)*x(6,n)]-[n1; n2(1,1); n2(2,1)]);
    %interation
    dx1(n)=x(4,n);
    dx2(n)=x(5,n);
    dx3(n)=x(6,n);
    dx4(n)=D(1);
    dx5(n)=D(2);
    dx6(n)=D(3);
    x(1,n+1)=x(1,n)+dx1(n)*t;
    x(2,n+1)=x(2,n)+dx2(n)*t;
    x(3,n+1)=x(3,n)+dx3(n)*t;
    x(4,n+1)=x(4,n)+dx4(n)*t;
    x(5,n+1)=x(5,n)+dx5(n)*t;
    x(6,n+1)=x(6,n)+dx6(n)*t;
    n=n+1;

end
a=t:t:(Te+0.001);
plot(a,x(1,:))












