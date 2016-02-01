X_dim = X(1);
Y_dim = X(2);
dat = zeros(X_dim, Y_dim);
theta = linspace(0,pi, X_dim);
phi = linspace(0, 2*pi, Y_dim);
for i = 1:X_dim
    for j = 1:Y_dim
        dat(i,j) = X(2+(i-1)*X_dim + j);
    end
end

x = sin(theta')*cos(phi);
y = sin(theta')*sin(phi);
z = cos(theta')*ones(1,Y_dim);

title('Python code method')
surf(x,y,z, dat, 'EdgeColor', 'None')
colorbar
shading flat

