function C = compute_C(num, num_points, radius)
l_range = num-1;
C = zeros(num,1);
a = -1.0;
b = 1.0;
pts = zeros([num_points+1,1]);
pts(1) = a;
for j=1:num_points
        pts(j+1) = a + ((b-a)*j)/num_points;
end
for i=0:l_range
    int_approx = 0.0;
    
    for j=1:num_points
        int_approx = int_approx+(func(pts(j), i) + func(pts(j+1),i));
    end
    C(i+1) = int_approx*(2.0*pi*radius/num_points);
end

end

function output = func(x,l)
    rho = exp(-3.0*acos(x)/8.5);
    L = legendreP(l, x);
    output = rho*L;
end
 