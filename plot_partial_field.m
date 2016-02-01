X_dim = X(1);
Y_dim = X(2);

%x_lb = 0.5759;
%x_ub = 0.6196;
%y_lb = 2.0420;
%y_ub = 2.094;

x_lb = 0.1778*pi;
x_ub = 0.2*pi;
y_lb = 1.327*pi;
y_ub = 1.355*pi;

dat = zeros(X_dim, Y_dim);
X_span = linspace(x_lb, x_ub, X_dim);
Y_span = linspace(y_lb, y_ub, Y_dim);
for i = 1:X_dim
    for j = 1:Y_dim
        dat(i,j) = X(2+(i-1)*X_dim + j);
    end
end
r = 1;
[theta, phi] = meshgrid(X_span, Y_span);
x1 = r.*cos(theta).*cos(phi);
y1 = r.*sin(theta).*cos(phi);
z1 = r.*sin(phi);
surf(x1,y1,z1,dat, 'Edgecolor', 'None')
colorbar