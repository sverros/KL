%h = 1E-3;
%a = 0.9;
%b = 1.0;

h = 1E-10;
a = 0.99999;
b = 1.0;
N = round((b-a)/h);
radius = 6371.0;

pts = zeros(N+1,1);
pts(1) = a;

for i=1:N
    pts(i+1) = a + i*h;
end

pts(end)

%%
integral1 = 0;
integral2 = 0;

l = 0;
op = zeros(N+1,1);
op(1) = exp(-3*radius*acos(pts(1))/8.5)*legendreP(l,pts(1));
for i=2:N
    op(i) = exp(-3*radius*acos(pts(i))/8.5)*legendreP(l,pts(i));
    integral1 = integral1 + (op(i)+op(i-1))*h*pi;
end

l=70;

for i=1:N
    op = exp(-3*radius*acos(pts(i))/8.5)*legendreP(l,pts(i));
    %op1 = exp(-3*radius*acos(pts(i+1))/8.5)*legendreP(l,pts(i+1));
    %integral2 = integral2 + (op+op1)*h*pi;
end