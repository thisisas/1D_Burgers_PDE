% hold on
x0 = 0;
x1 = 1;
t0 = 0;
t1 = 0.001;
dx = 1 / 40;
dt = 1 / 10000;
v = 1;
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
        U(i+1,j)=U(i,j) + ((dt)*(U(i,j)-U(i,j-1)))/(dx)+(v*(dt)*(U(i,j+1)-2*U(i,j)+U(i,j-1)))/(dx^2);
    end
end


x = linspace(0,1,40);
plot(x,U(end,:),'-'), hold on
xlabel('x')
ylabel('u(t,x)')
% text(0.5,0.5,'t_{f} = 0.')
% text(0.5,0.05,'t_{f} = 0.5')
% text(0.4,0.25,'t_{f} = 0.1')
% text(0.4,0.5,'t_{f} = 0.05')
% text(0.4,0.95,'t_{f} = 0.001')
title('v=1, \Deltat=0.0001,\Deltax=0.025')



