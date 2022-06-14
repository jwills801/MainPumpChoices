chi = .1;
k = 1.4;
gamma = .001:.001:.999;

slope1 = 1./gamma *(1+chi)^(1/k)/((1+chi)^(1/k)-1);
slope2 = 1./gamma *(1-chi)^(1/k)/((1-chi)^(1/k)-1);

figure, plot(gamma,slope1,gamma,-slope2), ylim([0 200])
ylabel('Slope'), xlabel('V_G over V_T at P_{rail}', 'Interpreter','Latex')
legend('Slope if \Delta V > 0','Negative of Slope if \Delta V < 0', 'Interpreter','Latex')