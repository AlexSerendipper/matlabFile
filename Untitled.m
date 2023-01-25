x0 = [1];

x = 0:0.1:4*pi;
y = cos(x);
% plot(y);

x_data = [];
y_data = [];
x_data = [x_data,x(20),x(30),x(40),x(50),x(60),x(70)];
y_data = [y_data,y(20),y(30),y(40),y(50),y(60),y(70)];

fun = @(x,x_data) x(1)*cos(x_data);

[x,resnorm]=lsqcurvefit(fun,x0,x_data,y_data);

plot(x_data,fun(x,x_data),x_data,y_data,'*');


