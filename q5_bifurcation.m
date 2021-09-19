format long
b = 0.3;      % coefficients in equation
n = 5000000;      % number of total steps
nplot = 50;         % number of points we plot
endpoint = 500;   % end point of t 
t = linspace(0,endpoint,n);  % starting and ending point of t
h = t(2)-t(1);    % step size
% i loop is used to find the numerical integral
for a = 0.1:0.0001:0.5      % increase a by 0.0001 every iteration
    x = zeros(n,1);
    y = zeros(n,1);
    f = @(t,x,y) y;     % formula for dx/dt
    g = @(t,x,y) -a*y+x-x^3+b*cos(t);     %formula for dx^2/dt^2
    x(1) = 0.1;     % initial value of x
    y(1) = 0.1;     % initial value of y=dx/dt
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
    
   % below is to find value of x,y at t=2k*pi
   max_multiple_pi = floor(endpoint/(2*pi));    % max number of multiples of pi that we've integrated numerically
   x_points = zeros(max_multiple_pi,1); % store the x and y points corresponding to t= 2k*pi
   y_points = zeros(max_multiple_pi,1);
    for k=1:max_multiple_pi
        x_points(k) = x(round(2*k*pi/h));    % find nearest point of x(t) where t=2k*pi
        y_points(k) = y(round(2*k*pi/h));   
    end
    % plot the final points at multiple of 2pi
    plot(a*ones(nplot,1), x_points(max_multiple_pi-nplot+1:max_multiple_pi),...
        '.', 'markersize', 2,'MarkerFaceColor','k','MarkerEdgeColor','k');
    hold on
end
