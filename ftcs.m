% hold on
x0 = 0;
x1 = 1;
t0 = 0;
t1 =0.4;
dx = 1 / 40;
dt = 1 / 1000;
v = 0.1;
nx = (x1-x0)/dx;
nt = (t1-t0)/dt;

U = zeros(nt,nx);

%IC
for i=1:1:nx
    x =(x0 +(i)*dx);
    U(1,i) = sin(pi*x);
end

%BC
for i=1:1:nt
%Left BC
    U(i,1) = 0;
%Right BC
    U(i,nx) = 0;
end

for i=1:1:nt-1
    for j=2:1:nx-1
        U(i+1,j)=U(i,j) + (dt)*(-U(i,j)*(U(i,j+1)-U(i,j-1))/(2*dx) + v*(U(i,j-1)-2*U(i,j)+U(i,j+1))/(dx^2) );
    end
end


x = linspace(0,1,40);
plot(x,U(end,:),'-o'), hold on
ylim([0,1]);

xlabel('x')
ylabel('u(t,x)')
% text(0.5,0.05,'t_{f} = 0.5')
% text(0.5,0.2,'t_{f} = 0.1')
% text(0.5,0.45,'t_{f} = 0.4')
% text(0.5,0.45,'v = 1')
title('v=0.1,\Deltat=0.001,\Deltax=0.025 FTCS')


%Exact solution

C.L = 1;
C.D = 0.1;
C.eta = 1;
C.K = 1;
C.Ip = 1;
m = 0;
x = linspace(0,1,40);
t = linspace(0,0.4,1000);
eqn = @(x,t,u,dudx) transistorPDE(x,t,u,dudx,C);
ic = @(x) transistorIC(x,C);
sol = pdepe(m,eqn,ic,@transistorBC,x,t);
plot(x,sol(end,:))
ylim([0,1]);

legend('Numerical','Analytical')

%
L2norm =sqrt((dx)*sum((U(end,:) -sol(end,:)).^2));
Linfty = norm((U(end,:) -sol(end,:)),inf);


%
function [c,f,s] = transistorPDE(x,t,u,dudx,C)
D = C.D;
eta = C.eta;
L = C.L;

c = 1;
f = D*dudx;
s = -(eta/L*u)*dudx;
end

function u0 = transistorIC(x,C)

u0 = sin(pi*x);
end

function [pl,ql,pr,qr] = transistorBC(xl,ul,xr,ur,t)
pl = ul;
ql = 0;
pr = ur;
qr = 0;
end
