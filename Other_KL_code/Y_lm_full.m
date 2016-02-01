function Y=Y_lm_full(l,N)

theta = linspace(0, pi, N);
phi = linspace(0, 2*pi, 2*N)';
Y = zeros(2*N, N, 2*l+1);

L = sqrt(1/(2*pi))*legendre(l, cos(theta), 'norm');

for m=-l:l
    if (m>0)
        Y(:,:, l+m+1) = sqrt(2)*cos(phi*m)*L(m+1, :);
    elseif (m<0)
        Y(:,:, l+m+1) = sqrt(2)*sin(phi*m)*L(-m+1, :);
    else       
        Y(:,:, l+m+1) = ones(2*N,1)*L(1, :);
    end
end

end