format long
b_r = [-0.1,-0.05,-0.01,-0.001];      % range of b to be plotted
m = numel(b_r);
a = 10;      % value of coef a
n = 50000000;      % number of total steps
endpoint = 5000;   % end point of t 
t = linspace(0,endpoint,n);  % starting and ending point of t
h = t(2)-t(1);    % step size
for i = 1:m
    b = b_r(i);
    x = zeros(n,1);
    y = zeros(n,1);
    f = @(t,x,y) y-a*(x^3/3-x);     % formula for x_dot
    g = @(t,x,y) -x+1+b;     %formula for y_dot
    x(1) = 0.5;     % initial value of x
    y(1) = 0.5;     % initial value of y=dx/dt
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
    plot(x(n-2000000:n),y(n-2000000:n),'DisplayName', sprintf('b=%g', b_r(i)))
    hold on
end
% to plot F(x)
x = linspace(-2.1,2.1,5000);
y = a.*((x.^3)./3-x);
plot(x,y,'DisplayName','F(x)')
xlabel('x')
ylabel('y')
legend()