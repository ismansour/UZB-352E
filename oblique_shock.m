function[M_2, beta, po2_o1, p_s] = oblique_shock(M_1, theta)

syms beta
gamma =1.4;

equ = tan(theta*pi/180)== ((2*cot(beta*pi/180)*(M_1^2*sin(beta*pi/180)^2 - 1))/((gamma + cos(2*beta*pi/180))*M_1^2 + 2));

sol_beta = [vpasolve(equ, beta,20), vpasolve(equ, beta,85)];
beta = double(sol_beta(1));

Mn_1 = M_1*sind(beta);
[Mn_2, p_s, T_s, po2_o1, po1_1, po2_1] = normal_shock(Mn_1, gamma);
M_2 = Mn_2/ sind(beta-theta);

end