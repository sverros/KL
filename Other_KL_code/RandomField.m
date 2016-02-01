function [x,y,z,T] = RandomField(alpha, N, L)
radius = 1.0;
T = zeros(2*N,N);
theta = linspace(0,pi,N);
phi = linspace(0, 2*pi, 2*N)';
T = (1/sqrt(2*pi))*randn(1,1)*1/2*sqrt(1/pi)*T;
%T = (1/sqrt(2*pi))*0.01*ones(1,1)*1/2*sqrt(1/pi)*T;
C = compute_C(L+1, 5000, radius);
%C = (1:L+1).^(-alpha);
for l=0:L
    %beta = ones(1,2*l+1)*0.01;
    beta = randn(1,2*l+1);
    Y = Y_lm_full(l,N);
    for m=-l:l
        T = T + sqrt(C(l+1))*beta(1, l+m+1)*Y(:,:,l+m+1);
    end
end

x = (cos(phi)*sin(theta));
y = (sin(phi)*sin(theta));
z = (ones(2*N,1)*cos(theta));

surf(x,y,z,T)
shading flat
colorbar

end