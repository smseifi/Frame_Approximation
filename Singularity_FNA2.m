% --- Adaptive Frame Numerical Approximation --- %
% FNA2-type Computation - To capture singularities
% Ben Adcock - Mohsen Seifi, Simon Fraser University, Department of Math

%   -------------------------------------------------------------------   %
c = input(['Choose one of the following functions \n' ... 
'1: f_1 = exp(sin(15*t+.5))+log(t).*cos(t) \n' ... 
'2: f_2(t) = exp(sin(15*t+.5))+log(t).*cos(20*t) \n' ...    
'3: f_3(t) = exp(sin(15*t+.5))+log(t).*cos(40*t) \n']);

switch c
    case 1
        f = @(t) exp(sin(15*t+.5))+log(t).*cos(t);
        name = 'f_1'; 
    case 2
        f = @(t) exp(sin(15*t+.5))+log(t).*cos(20*t);
        name = 'f_2';
    case 3
        f = @(t) exp(sin(15*t+.5))+log(t).*cos(40*t);
        name = 'f_3';
    otherwise
        disp('Enter 1 to choose f_1, 2 to choose f_2, or 3 to choose f_3');
end

C1Vect = [1.1, 2, 5, 15, 100]';   C2Vect = [2.1, 3, 5, 15, 100]';
C1L = length(C1Vect);   C2L = length(C2Vect);

epsVect = [1e-8, 1e-14, 1e-15]';  eL = length(epsVect);
tol = epsVect(3, 1);
%%
NMax = 300; step = 5;   NVect = (step:step:NMax)';  NL= length(NVect);
gamma = 2;  MVect = gamma*NVect; k = 5;

K = 10*NMax; ErrGrid = linspace(1e-15, 1, K)';  fGrid = f(ErrGrid);

Err_TSVD = zeros(eL, NL); coefNorm_TSVD = zeros(eL, NL);
Err_ASVD1 = zeros(C1L, NL);  coefNorm_ASVD1 = zeros(C1L, NL); 
Err_ASVD2 = zeros(C2L, NL);  coefNorm_ASVD2 = zeros(C2L, NL);

Rho_Err = SystemMatrixLeg(0:1:NMax-1, 2*ErrGrid-1);
Psi_Err = zeros(K, k);
for i = 1:k
    Psi_Err(:, i) = log(ErrGrid).*Rho_Err(:, i);
end


for j = 1:NL
    
    N = NVect(j, 1); 
    E = [Rho_Err(:, 1:N), Psi_Err];
    
    M = MVect(j, 1);    mVect = (1:1:M)';
    ChebGrid = (cos(pi*(2*mVect-1)/(2*M))+1)/2;
    
    ChebGridA = [1 ; ChebGrid];
    w = ChebGridA(1:M)-ChebGrid; 
    
    Rho = SystemMatrixLeg(0:1:N-1, 2*ChebGrid-1);
    Psi = zeros(M, k);
    for i = 1:k
        Psi(:, i) = log(ChebGrid).*Rho(:, i);
    end
    Phi = diag(sqrt(w))*[Rho, Psi];     [U, S, V] = svd(Phi);
   
    l = f(ChebGrid).*sqrt(w); 
    
    for i_T = 1:eL
        x_TSVD = TSVD(U, S, V, l, epsVect(i_T, 1));
        coefNorm_TSVD(i_T, j) = norm(x_TSVD);
        Err_TSVD(i_T, j) = norm(E*x_TSVD-fGrid)/sqrt(K);
    end
 
    for i_1 = 1:C1L
        C = C1Vect(i_1, 1)*norm(l, 2);
        [x_ASVD1, ~] = ASVD1(U, S, V, l, C, tol);
        coefNorm_ASVD1(i_1, j) = norm(x_ASVD1);
        Err_ASVD1(i_1, j) = norm(E*x_ASVD1-fGrid)/sqrt(K);
    end
    
    for i_2 = 1:C2L
        C = (C2Vect(i_2, 1)*norm(l, 2))^2;
        [x_ASVD2, ~] = ASVD2(U, S, V, l, C, tol);
        coefNorm_ASVD2(i_2, j) = norm(x_ASVD2);
        Err_ASVD2(i_2, j) = norm(E*x_ASVD2-fGrid)/sqrt(K);
    end
end
%% 

j = find(epsVect == tol); 

default_position = [675, 548, 570, 505];
default_color = [0    0.4470    0.7410 ;
    0.8500    0.3250    0.0980 ;
    0.9290    0.6940    0.1250 ;
    0.4940    0.1840    0.5560 ;
    0.4660    0.6740    0.1880 ;
    0.3010    0.7450    0.9330 ;
    0.6350    0.0780    0.1840 ];

Err_Ylim = [1e-15, 1e4];  Err_YTick = 10.^(-14:2:4);
coef_Ylim = [1e0, 1e11];   coef_YTick = 10.^(0:2:10);
NTick = 0:50:300;

%% ASVD1_Err plot
f1 = figure(1);
for i = 1:C1L
semilogy(NVect, Err_ASVD1(i, :), '-', 'Color',default_color(i, :), ...
    'LineWidth', 2, 'DisplayName', ...
    strcat('$c = $', num2str(C1Vect(i, 1))));
hold on
end
semilogy(NVect, Err_TSVD(j, :), 'k--', 'LineWidth', 1, ...
    'DisplayName', 'TSVD');

legend('Location', 'northeast', 'Interpreter', 'LaTex', 'LineWidth', ...
    1.1, 'FontSize', 22);

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
saveas(figure(1), strcat(name, '_Sing_FNA2_ASVD1_err'), 'epsc');
saveas(figure(1), strcat(name, '_Sing_FNA2_ASVD1_err.fig'));
%% ASVD1_Coef plot
f3 = figure(3);
for i = 1:C1L
semilogy(NVect, coefNorm_ASVD1(i, :), '-','Color',default_color(i, :), ...
    'LineWidth', 2, 'DisplayName', ...
    strcat('$c = $', num2str(C1Vect(i, 1))));
hold on
end

semilogy(NVect, coefNorm_TSVD(j, :), 'k--', 'LineWidth', 1, ...
    'DisplayName', 'TSVD');

legend('Location', 'northeast', 'Interpreter', 'LaTex', 'LineWidth',  ...
    1.1, 'FontSize', 22);

ax = gca;   ax.FontSize = 15;   ax.LineWidth = 2;

ax.XTick = NTick;
ax.XGrid = 'on';
ax.XMinorTick = 'off';   ax.XMinorGrid = 'on';

ax.YTick = coef_YTick;   ax.YLim = coef_Ylim; 
ax.YGrid = 'on'; 
ax.YMinorTick = 'off';  ax.YMinorGrid = 'on';
   
ax.Units = 'normalized';
f3.InnerPosition = default_position;
f3.OuterPosition = default_position;

hold off

saveas(figure(3), strcat(name, '_Sing_FNA2_ASVD1_coef'), 'epsc');
saveas(figure(3), strcat(name, '_Sing_FNA2_ASVD1_coef.fig'));
%% ASVD2_Err plot
f2 = figure(2);
for i = 1:C2L
semilogy(NVect, Err_ASVD2(i, :), '-', 'Color',default_color(i, :), ...
    'LineWidth', 2, 'DisplayName', ...
    strcat('$c = $', num2str(C2Vect(i, 1))));
hold on
end
semilogy(NVect, Err_TSVD(j, :), 'k--', 'LineWidth', 1, ...
    'DisplayName', 'TSVD');

legend('Location', 'northeast', 'Interpreter', 'LaTex', 'LineWidth', ...
    1.1, 'FontSize', 22);

ax = gca;   ax.FontSize = 15;   ax.LineWidth = 2;

ax.XTick = NTick;
ax.XGrid = 'on';
ax.XMinorTick = 'off';   ax.XMinorGrid = 'on';

ax.YTick = Err_YTick;   ax.YLim = Err_Ylim; 
ax.YGrid = 'on'; 
ax.YMinorTick = 'off';  ax.YMinorGrid = 'on';
   
ax.Units = 'normalized';
f2.InnerPosition = default_position;
f2.OuterPosition = default_position;

hold off

saveas(figure(2), strcat(name, '_Sing_FNA2_ASVD2_err'), 'epsc');
saveas(figure(2), strcat(name, '_Sing_FNA2_ASVD2_err.fig'));

%% ASVD2_Coef plot

f4 = figure(4);
for i = 1:C2L
semilogy(NVect, coefNorm_ASVD2(i, :), '-','Color',default_color(i, :), ...
    'LineWidth', 2, 'DisplayName', ...
    strcat('$c = $', num2str(C2Vect(i, 1))));
hold on
end

semilogy(NVect, coefNorm_TSVD(j, :), 'k--', 'LineWidth', 1, ...
    'DisplayName', 'TSVD');

legend('Location', 'northeast', 'Interpreter', 'LaTex', 'LineWidth',  ...
    1.1, 'FontSize', 22);

ax = gca;   ax.FontSize = 15;   ax.LineWidth = 2;

ax.XTick = NTick;
ax.XGrid = 'on';
ax.XMinorTick = 'off';   ax.XMinorGrid = 'on';

ax.YTick = coef_YTick;  ax.YLim = coef_Ylim; 
ax.YGrid = 'on'; 
ax.YMinorTick = 'off';  ax.YMinorGrid = 'on';
   
ax.Units = 'normalized';
f4.InnerPosition = default_position;
f4.OuterPosition = default_position;

hold off

saveas(figure(4), strcat(name, '_Sing_FNA2_ASVD2_coef'), 'epsc');
saveas(figure(4), strcat(name, '_Sing_FNA2_ASVD2_coef.fig'));