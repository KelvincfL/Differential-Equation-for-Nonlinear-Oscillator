format long
initial_x = [-2,-0.1,0.1,0.5,0.5];   % choice of initial condition for x(0)
initial_y = [2,-0.1,0.1,1.5,-1.5];     % choice of initial condition for dx/dt(0)
a = 1.5;
b = 1;      % coefficients in equation
n = 5000000;      % number of total steps
endpoint = 500;   % end point of t 
t = linspace(0,endpoint,n);  % starting and ending point of t
h = t(2)-t(1);    % step size
for i = 1:5
    x = zeros(n,1);
    y = zeros(n,1);
    f = @(t,x,y) y;     % formula for dx/dt
    g = @(t,x,y) a*x-x^3+(b-x^2)*y;     %formula for dx^2/dt^2
    x(1) = initial_x(i);     % initial value of x
    y(1) = initial_y(i);     % initial value of y=dx/dt
    for j=1:n-1
        k1 = f(t(j),x(j),y(j));     % k is the RK4 for x
        l1 = g(t(j),x(j),y(j));     % l is the RK4 for y=dx/dt

        k2 = f(t(j)+0.5*h,x(j)+0.5*h*k1,y(j)+0.5*h*l1);
        l2 = g(t(j)+0.5*h,x(j)+0.5*h*k1,y(j)+0.5*h*l1);

        k3 = f(t(j)+0.5*h,x(j)+0.5*h*k2,y(j)+0.5*h*l2);
        l3 = g(t(j)+0.5*h,x(j)+0.5*h*k2,y(j)+0.5*h*l2);

        k4 = f(t(j)+h,x(j)+h*k3,y(j)+h*l3);
        l4 = g(t(j)+h,x(j)+h*k3,y(j)+h*l3);

        x(j+1) = x(j)+(1/6)*h*(k1+2*k2+2*k3+k4); %updating values of x
        y(j+1) = y(j)+(1/6)*h*(l1+2*l2+2*l3+l4);    %updating values of y
    end
    plot(x,y,'DisplayName', sprintf('(%g,%g)', initial_x(i), initial_y(i)))
    hold on
end
xlabel('x')
ylabel('dx/dt')
legend()
xlim([-2,2])
ylim([-2,2])