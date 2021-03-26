% --- Frame Approximation with Bounded Coefficients --- %
% Exploding Coefficients!
% Ben Adcock - Mohsen Seifi, Department of Math @ Simon Fraser University
%   -------------------------------------------------------------------   %

c = input(['Choose one of the following functions' ... 
'\n 1: f_1(t) = 1/(1+75t^2) \n 2: f_2(t) = 1/(.57-t)' ...    
'\n 3: f_3(t) = exp(sin(20*t+.5))*sqrt(1+t)+cos(10*t) \n']);
switch c
    case 1
        f = @(t) 1./(1+75.*t.^2);
        name = 'f_1';
    case 2
        f = @(t) 1./(.57-t);
        name = 'f_2';
    case 3
        f = @(t) exp(sin(20*t+.5)).*sqrt(1+t)+cos(10*t);
        name = 'f_3';
    otherwise
        disp('Enter 1 to choose f_1, 2 to choose f_2, or 3 to choose f_3');
end
%% Setting Parameters - Declaring Variables 

c_ASVD1 = 15;  eps_TSVD = 1e-15;   lambda_tikhonov = 1e-15; 

omega = 0.5;    

NMax = 300; step = 5; NVect = (step:step:NMax)';   NL = length(NVect);

R = 100*NMax;  K = 10*NMax;
RiemannGrid = linspace(-omega, omega, R)'; ...
    ErrGrid = linspace(-omega, omega, K)';  fGrid = f(ErrGrid);

Err_TSVD = zeros(NL ,1);        coefNorm_TSVD = zeros(NL ,1); 
Err_Tikhonov = zeros(NL ,1);    coefNorm_Tikhonov = zeros(NL ,1); 
Err_Backslash = zeros(NL ,1);   coefNorm_Backslash = zeros(NL ,1);
Err_QR = zeros(NL ,1);          coefNorm_QR = zeros(NL ,1);
Err_ASVD1 = zeros(NL ,1);       coefNorm_ASVD1 = zeros(NL ,1);  
%% Setting the Error and the Super-Gram Matrix

E = SystemMatrixLeg(0:1:(NMax-1), ErrGrid)/sqrt(2);
Phi = SystemMatrixLeg(0:1:(NMax-1), RiemannGrid)/sqrt(2);

G_NMax = (2*omega/R)*(Phi'*Phi); y = (2*omega/R)*Phi'*f(RiemannGrid);
%% Main Computation

for j = 1:NL
    N = NVect(j, 1);
    E_N = E(:, 1:N);  y_N = y(1:N); ...
        G_N = G_NMax(1:N, 1:N); [U_N, S_N, V_N] = svd(G_N); ...
        [Q_N, R_N] = qr(G_N);
    
    x_TSVD = TSVD(U_N, S_N, V_N, y_N, eps_TSVD);
    coefNorm_TSVD(j, 1) = norm(x_TSVD, 2);
    Err_TSVD(j, 1) = norm(E_N*x_TSVD-fGrid)/sqrt(K);
    
    x_Tikhonov = Tikhonov(U_N, S_N, V_N, y_N, lambda_tikhonov);
    coefNorm_Tikhonov(j, 1) = norm(x_Tikhonov, 2);
    Err_Tikhonov(j, 1) = norm(E_N*x_Tikhonov-fGrid)/sqrt(K);
    
    x_Backslash = G_N\y_N;
    coefNorm_Backslash(j, 1) = norm(x_Backslash, 2);
    Err_Backslash(j, 1) = norm(E_N*x_Backslash-fGrid)/sqrt(K);
    
    x_QR = R_N\(Q_N'*y_N);
    coefNorm_QR(j, 1) = norm(x_QR, 2);
    Err_QR(j, 1) = norm(E_N*x_QR-fGrid)/sqrt(K);
    
    [x_ASVD1, ~] = ASVD1(U_N, S_N, V_N, y_N, c_ASVD1, eps_TSVD);
    coefNorm_ASVD1(j, 1) = norm(x_ASVD1, 2);
    Err_ASVD1(j, 1) = norm(E_N*x_ASVD1-fGrid)/sqrt(K);
end
%% Setting up plotting parameters
default_position = [675, 548, 570, 505];

default_color = [0    0.4470    0.7410 ;
    0.8500    0.3250    0.0980 ;
    0.9290    0.6940    0.1250 ;
    0.4940    0.1840    0.5560 ;
    0.4660    0.6740    0.1880 ;
    0.3010    0.7450    0.9330 ;
    0.6350    0.0780    0.1840 ];

NTick = 0:50:300;

Err_Ylim = [1e-10, 1e0];    Err_YTick = 10.^(-10:2:0);
coef_Ylim = [1e-1, 1e7]; coef_YTick = 10.^(0:2:6);
%% Error Plot
f1 = figure(1);

semilogy(NVect, Err_Tikhonov, '-', 'LineWidth', 2.2, 'Color', ...
    default_color(1, :), 'DisplayName', 'Tikhonov');
hold on
semilogy(NVect, Err_Backslash, '-', 'LineWidth', 2.2, 'Color', ...
    default_color(2, :), 'DisplayName', 'Backslash');
hold on
semilogy(NVect, Err_QR, '-', 'LineWidth', 2.2, 'Color', ...
    default_color(3, :), 'DisplayName', 'QR');
hold on
semilogy(NVect, Err_ASVD1, '-', 'LineWidth', 2.2, 'Color', ...
    default_color(4, :), 'DisplayName', 'ASVD1');
hold on
semilogy(NVect, Err_TSVD, 'k--', 'LineWidth', 1.5, 'DisplayName', 'TSVD');

legend('Location', 'northeast', 'Interpreter', 'LaTex', 'LineWidth', ...
    1.1, 'FontSize', 18);

ax = gca;   ax.FontSize = 15;   ax.LineWidth = 2;

ax.XTick = NTick;
ax.XGrid = 'on';
ax.XMinorTick = 'off';   ax.XMinorGrid = 'on';

ax.YTick = Err_YTick;   ax.YLim = Err_Ylim; 
ax.YGrid = 'on'; 
ax.YMinorTick = 'off';  ax.YMinorGrid = 'on';

ax.Units = 'normalized';
f1.InnerPosition = default_position;
f1.OuterPosition = default_position;

hold off

saveas(figure(1), strcat(name, '_exploding_err'), 'epsc');
saveas(figure(1), strcat(name, '_exploding_err.fig'));
%% Coef. Plot

f2 = figure(2);

semilogy(NVect, coefNorm_Tikhonov, '-', 'LineWidth', 2.2, 'Color', ...
    default_color(1, :), 'DisplayName', 'Tikhonov');
hold on
semilogy(NVect, coefNorm_Backslash, '-', 'LineWidth', 2.2, 'Color', ...
    default_color(2, :), 'DisplayName', 'Backslash');
hold on
semilogy(NVect, coefNorm_QR, '-', 'LineWidth', 2.2, 'Color', ...
    default_color(3, :), 'DisplayName', 'QR');
hold on
semilogy(NVect, coefNorm_ASVD1, '-', 'LineWidth', 2.2, 'Color', ...
    default_color(4, :), 'DisplayName', 'ASVD1');
hold on
semilogy(NVect, coefNorm_TSVD, 'k--', 'LineWidth', 1.5, 'DisplayName', ...
    'TSVD');

legend('Location', 'northeast', 'Interpreter', 'LaTex', 'LineWidth', ...
    1.1, 'FontSize', 18);

ax = gca;   ax.FontSize = 15;   ax.LineWidth = 2;

ax.XTick = NTick;
ax.XGrid = 'on';
ax.XMinorTick = 'off';   ax.XMinorGrid = 'on';

ax.YTick = coef_YTick;   ax.YLim = coef_Ylim; 
ax.YGrid = 'on'; 
ax.YMinorTick = 'off';  ax.YMinorGrid = 'on';

ax.Units = 'normalized';
f2.InnerPosition = default_position;
f2.OuterPosition = default_position;

hold off

saveas(figure(2), strcat(name, '_exploding_coef'), 'epsc');
saveas(figure(2), strcat(name, '_exploding_coef.fig'));

