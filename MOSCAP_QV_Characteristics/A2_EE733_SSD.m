%%%% ASSIGNMENT 2 : MOSCAP Charge Profile %%%%%

epsilon_Si = 11.8; % Dielectric constant of Silicon
epsilon_ox = 3.9; % Dielectric constant of Silicon
epsilon_0 = 8.85 * 10^(-14); % Permittivity of free space in F/cm
n_i = 1.5 * 10^10; % Intrinsic carrier concentration per cm^3
N_a = 1 * 10^17; % Acceptor doping concentration per cm^3
q = 1.6 * 10^(-19); % Charge of an electron
V_t = 0.026; % Thermal voltage in volts
t_ox = (2+5) * 10^(-7); % Oxide thickness in cm
h = 0.01; % Voltage step size

%%%% Total Charge ( C / cm^2 ) with delta-depletion approximation %%%%

C_ox = (epsilon_ox*epsilon_0)/t_ox; % Oxide capacitance per unit area in F/cm^2
V_g = -5:h:5; % Gate voltage
psi_f = V_t * log(N_a/n_i); % Fermi potential in V
V_Th = sqrt(4*q*N_a*epsilon_Si*epsilon_0*psi_f)/C_ox + 2 * psi_f; % Threshold voltage in V
Q_a = -1*(C_ox * V_g(1,1:round((5/h)+1))); % Charge in accumulation region in C/cm^2
Q_d = abs((q*epsilon_Si*epsilon_0*N_a) * (1-sqrt(1+(2*C_ox^2*V_g(1,round((5/h)+2):round(((5+V_Th)/h)+1)))/(q*epsilon_Si*epsilon_0*N_a)))/C_ox); % Charge in depletion region in C/cm^2
Q_i = C_ox * (V_g(1,round(((5+V_Th)/h)+1):length(V_g))-V_Th) + sqrt(4*q*N_a*epsilon_Si*epsilon_0*psi_f); % Charge in inversion region in C/cm^2
Q_i(1,1) = Q_d(1,length(Q_d)); % for continuous curve

%%%% Newton-Raphson for determining Surface Potential, psi_s %%%%

psi_s = zeros(1,length(V_g));
for k = 1:length(V_g)
    V = V_g(1,k);
    if V >= 0
        x = 0.01;
    else
        x = -0.01;
    end
    for i = 1:1000
        g = 2*q*epsilon_Si*epsilon_0*V_t*(N_a*(exp(-x/V_t)-1)+(n_i^2/N_a)*(exp(x/V_t)-1)+(N_a*x)/V_t);
        e = 2*q*epsilon_Si*epsilon_0*(-N_a*exp(-x/V_t)+(n_i^2/N_a)*exp(x/V_t)+N_a);
        m = C_ox*(V-x);
        f = g-m^2;
        f_1 = e - 2*C_ox^2*(x-V);
        x = x - (f/f_1);
    end
    psi_s(1,k) = x;
end

%%%% Total Charge ( C / cm^2 ) without delta-depletion approximation %%%%

Q_s = sqrt(2*q*epsilon_Si*epsilon_0*V_t*(N_a*(exp(-psi_s/V_t)-1)+(n_i^2/N_a)*(exp(psi_s/V_t)-1)+(N_a*psi_s)/V_t));
Q_s(1,round((5/h)+1)) = 0; % No charge for zero gate voltage

%%%% Plot of Total Charge %%%%
figure;
semilogy(V_g(1,1:round((5/h)+1)),Q_a,'Displayname','Accumulation Region approximation');
hold on;
semilogy(V_g(1,round((5/h)+2):round(((5+V_Th)/h)+1)),Q_d,'Displayname','Depletion Region approximation');
semilogy(V_g(1,round(((5+V_Th)/h)+1):length(V_g)),Q_i,'Displayname','Inversion Region approximation');
semilogy(V_g,Q_s,'Displayname','Total Charge ( without approximation )');
xlabel('Gate Voltage ( V )');
ylabel('Total Charge ( C / cm^2 )');
title('Surface Charge Density vs. Gate Voltage');
legend;
grid on;
hold off;