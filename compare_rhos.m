x = linspace(0, pi, 1000);
plot( x, (x+1).^(-2.1), x, (x+1).^(-4.1), x, exp(-3*x/8.5), x, exp(-3*6371*x/8.5))
legend('2.1', '4.1', 'exp', 'exp_radius')

%%

l = linspace(0,100,101);
al1 = (l+1).^(-0.001);
al2 = (l+1).^(-0.01);
al3 = (l+1).^(-0.9);
al4 = (l+1).^(-1);

al5 = (l+1).^(-2);
al6 = (l+1).^(-4);


plot(l, al1, l, al2, l, al3, l, al4, l, al5, l, al6)
legend('1', '2', '3', '4', '5', '6')

%%
l1 = linspace(0, 0.01, 100);
l2 = linspace(0.01, .4, 10000);
%plot(l, (l).^(-2), l, exp(-3*l*6371/8.5), l, (l+1).^(-2))
hold on
plot(l1, ones(100,1))
plot(l2, exp(-3*(l2-0.01)*6371/8.5))
%legend('al', 'exp')