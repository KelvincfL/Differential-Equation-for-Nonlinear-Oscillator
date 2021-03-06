format long
b_r = [-0.0005,-0.0002];      % range of b to be plotted
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
    f = @(t,x,y) y;     % formula for dx/dt
    g = @(t,x,y) 1+b-x-a*(x^2-1)*y;     %formula for dx^2/dt^2
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
    % we plot t=4990 to 5000 for a=1
    % we plot t=4900 to 5000 for a=5
    % we plot t=4800 to 5000 for a=10
    % this is larger than in q8 since when b is small the period is large
    plot(x(n-2000000:n),y(n-2000000:n),'DisplayName', sprintf('b=%g', b_r(i)))
    hold on
end
xlabel('x')
ylabel('dx/dt')
legend()