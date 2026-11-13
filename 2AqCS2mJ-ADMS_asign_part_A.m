clearvars 
close all
clc
%% Known parameters 
L = 1.2 ;    %[m]
b = 0.04 ;   %[m]
h = 0.008 ;  %[m]
rho = 2700 ; %[kg/m^3]
E = 68*1e9 ; %[N/m^2]
%% Finding firsts naturals frequencies of H (roots of function)
m = rho*b*h ;    %[kg/m]
J = (b*h^3)/12 ; %[m^4]
EJ = E*J ;      %[N*m^2]
gamma = @(w) (m*w^2/EJ)^0.25;
H = @(w)[ 1                0                1                0             ;
    0                1                0                1             ;
    -cos(gamma(w)*L) -sin(gamma(w)*L) cosh(gamma(w)*L) sinh(gamma(w)*L) ; 
    sin(gamma(w)*L) -cos(gamma(w)*L) sinh(gamma(w)*L) cosh(gamma(w)*L) ; ] ;
detH = @(w) det(H(w)) ;
w_sample = linspace(0, 1000, 5000);
detH_val = arrayfun(detH, w_sample);
% Balayage pour localiser les changements de signe
index = find(detH_val(1:end-1).*detH_val(2:end) < 0);  
nat_puls = zeros(numel(index),1);
for k = 1:numel(index)
    nat_puls(k) = fzero(detH, [w_sample(index(k)), w_sample(index(k)+1)]);
end
nat_freq = nat_puls/(2*pi) ;
%% Mode shapes - shapes functions
Nmode = 4 ;
x_grid = linspace(0, L,500);
figure; hold on
for i = 1:Nmode
    wi  = nat_puls(i)  ;
    gi  = gamma(wi) ; 
    Hi = H(wi) ; 
    z  = [-Hi(1:3,1:3) \ Hi(1:3,4) ; 1] ;  %normalized with last const.
    modeShape = z(1)*cos(gi*x_grid) + z(2)*sin(gi*x_grid) + z(3)*cosh(gi*x_grid) + z(4)*sinh(gi*x_grid) ;
    normalized_mode =  modeShape / max(abs(modeShape));
    plot(x_grid,normalized_mode,'LineWidth',1.8,'DisplayName',sprintf('Mode %d : f = %.2f Hz',i,nat_freq(i)));
end
% Labels and plot appearance
yline(0, '--', 'HandleVisibility', 'off');   
xlabel('Beam position x [m]');
ylabel('Normalized displacement \Phi_i(x)');
title('First 4 Mode Shapes of the Cantilever Beam');
legend('Location', 'best');
hold off;
%% FRF 

xk = L ;                            %[m] - Position of force applied
xj = [0.2; 0.4; 0.5; 0.8; 1; 1.2] ; %[m] - avoiding mode shapes nodes   
ksi = 0.01*ones(1,Nmode) ;          % damping of 1%
wn = nat_puls ;                     % Vector of 4 natural frequencies
gam= (m*wn.^2/EJ).^0.25 ;           % Vector of 4 gamma values associated to w_i
BC  = -(cos(gam*L) + cosh(gam*L)) ./ (sin(gam*L) + sinh(gam*L));                       % Boundary constant values and
Phi = @(x,i) cos(gam(i)*x) - cosh(gam(i)*x) + BC(i)*(sin(gam(i)*x) - sinh(gam(i)*x));  % analytical solution
mi = zeros(1,Nmode);                % modal mass of mode 4 1st modes              
for i = 1:Nmode 
    mi(i) = trapz(x_grid, m*Phi(x_grid,i).^2);
end 
function G = frf(xj, xk, w_samp , Phi, mi, ksi, wn, Nmode)
G = zeros(numel(xj), numel(w_samp));                                    % M rows of FRF 
for i = 1:Nmode                                                         % by N columns of evaluating points for plot,
    num = Phi(xj, i) .* Phi(xk, i) / mi(i) ; 
    den = -w_samp.^2 + 1i*2*ksi(i)*wn(i)*w_samp + wn(i)^2;
    G = G + (num ./ den);
end
end
G_exp = frf(xj, xk, w_sample, Phi, mi, ksi, wn, Nmode);

figure; grid on    
subplot(2,1,1)            
semilogy(w_sample/(2*pi), abs(G_exp(1,:)), ...
    w_sample/(2*pi), abs(G_exp(3,:)), ...
    w_sample/(2*pi), abs(G_exp(end,:)), 'LineWidth', 1.2) ; 
ylabel('|G_{jk}| [m/N]');
legend('x_j = 0.2 m', 'x_j = 0.4 m', 'x_j = 1.2 m (colocated with xk)', 'Location', 'best');
title('FRF - Loaded at the end : x_k = L');

subplot(2,1,2)            
plot(w_sample/(2*pi), unwrap(angle(G_exp(1,:))), ...
    w_sample/(2*pi), unwrap(angle(G_exp(3,:))), ...
    w_sample/(2*pi), unwrap(angle(G_exp(end,:))), 'LineWidth', 1.2); 
xlabel('Frequency [Hz]'); ylabel('[rad]');


%% Experimental part : identification des 4 premiers modes (n = 1)

% Data preparation
% G_exp est nFRF x M frequences -> prendre .' pour avoir M pulsations x nFRF FRF (frf colonne par colonne)
nFRF = size(G_exp, 1);                 % Nombre de FRF

%% Structure du vecteur des inconnues (identique pour chaque mode)
ind = zeros(5, nFRF);              % indices de chaque bloc dans sol
for ii = 1:5
    ind(ii,:) = 2 + ii + (nFRF-1)*(ii-1) + (0:(nFRF-1));
end
nX = 2 + 5*nFRF;

options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt');

sol_list     = zeros(nX, Nmode);     % colonne k = parametres identifies du mode k
idx_fit_cell = cell(1, Nmode);       % indices de la bande de fitting du mode k
w_fit_cell   = cell(1, Nmode);       % pulsations de la bande du mode k
G_fit_cell   = cell(1, Nmode);       % FRF identifiee sur la bande (S_k x nFRF)
err_rel      = zeros(1, Nmode);      % erreur relative FRF dans la bande

for k = 1:Nmode
    % Bande de fitting autour du mode k (+/- 20 %, sans toucher w = 0 a cause de R_L/w^2)
    w_target = nat_puls(k);
    w_min = max(0.8*w_target, w_sample(2));
    w_max = min(1.2*w_target, w_sample(end));
    idx_fit = find(w_sample >= w_min & w_sample <= w_max);
    w_fit   = w_sample(idx_fit);            % vecteur colonne
    FRF_exp = G_exp(:, idx_fit).';          % sous-matrice PSEUDO EXP (S x nFRF)
    idx_fit_cell{k} = idx_fit;
    w_fit_cell{k}   = w_fit;

    % Numerical FRF + fonction d'erreur (chainage du cours)
    FRF_num = @(x) compute_FRF_num(x, w_fit, nFRF, ind);
    fun1   = @(x) FRF_num(x);
    fun2   = @(x) fun1(x) - FRF_exp ;
    err_re = @(x) real(fun2(x));
    err_im = @(x) imag(fun2(x));
    err    = @(x) [err_re(x); err_im(x)];

    % Initial guess : w_i et xi_i au pic, A_r au pic de chaque FRF, residus nuls
    w0  = w_target;
    xi0 = 0.01 ;
    A0  = max(abs(FRF_exp)) * (2*xi0*w0^2);   % 1 x nFRF
    X0  = [xi0, w0, A0, zeros(1, 4*nFRF)];    % 2 + 5*nFRF elements

    % minimization : les nFRF FRF sont ajustees simultanement
    sol = lsqnonlin(@(x) err(x), X0, [], [], options) ;
    sol_list(:,k) = sol;
    G_fit_cell{k} = FRF_num(sol);             % fit uniquement sur la bande w_fit
    err_rel(k)    = norm(fun2(sol),'fro') / norm(FRF_exp,'fro');
end

% Parametres identifies
xi_id  = sol_list(1,:);                        % 1 x Nmode
w_id   = sol_list(2,:);
f_id   = w_id/(2*pi);                          % Hz
A_id   = sol_list(ind(1,:), :);                % nFRF x Nmode : colonne k = deformee du mode k
% A_j ~ Phi(xj)*Phi(L)/m_i : le signe de Phi(L) alterne d'un mode a l'autre
s_tip  = sign(Phi(L, 1:Nmode));                % 1 x Nmode
Phi_id = A_id .* s_tip ./ max(abs(A_id), [], 1);   % deformees normalisees, signe aligne

%% Item 4 : qualite de l'identification - FRF exp sur tout le domaine, les 4 fits sur leur bande
ph_exp = unwrap(angle(G_exp), [], 2);          % phase deroulee sur tout le domaine (nFRF x M)
cmap   = lines(Nmode);

figure('Name', 'Identification - FRF experimentale vs identifiee'); grid on
subplot(2,1,1); hold on
for r = 1:nFRF
    if r == 1
        semilogy(w_sample/(2*pi), abs(G_exp(r,:)), 'b-', 'LineWidth', 0.8, 'DisplayName', 'G_{jk}^{EXP} (domaine complet)');
    else
        semilogy(w_sample/(2*pi), abs(G_exp(r,:)), 'b-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end
for k = 1:Nmode
    for r = 1:nFRF
        if r == 1
            semilogy(w_fit_cell{k}/(2*pi), abs(G_fit_cell{k}(:,r)), 'Color', cmap(k,:), 'LineWidth', 1.6, ...
                'DisplayName', sprintf('Fit mode %d : %.2f Hz', k, f_id(k)));
        else
            semilogy(w_fit_cell{k}/(2*pi), abs(G_fit_cell{k}(:,r)), 'Color', cmap(k,:), 'LineWidth', 1.6, 'HandleVisibility', 'off');
        end
    end
end
ylabel('|G_{jk}| [m/N]');
title('FRF (x_k = L) : experimentales sur tout le domaine, fits sur leur bande');
legend('Location', 'best'); grid on

subplot(2,1,2); hold on
for r = 1:nFRF
    if r == 1
        plot(w_sample/(2*pi), ph_exp(r,:), 'b-', 'LineWidth', 0.8, 'DisplayName', 'G_{jk}^{EXP} (domaine complet)');
    else
        plot(w_sample/(2*pi), ph_exp(r,:), 'b-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end
for k = 1:Nmode
    for r = 1:nFRF
        ph_fit = unwrap(angle(G_fit_cell{k}(:,r)));
        ph_fit = ph_fit + (ph_exp(r, idx_fit_cell{k}(1)) - ph_fit(1));   % recalage de la phase
        if r == 1
            plot(w_fit_cell{k}/(2*pi), ph_fit, 'Color', cmap(k,:), 'LineWidth', 1.6, ...
                'DisplayName', sprintf('Fit mode %d : %.2f Hz', k, f_id(k)));
        else
            plot(w_fit_cell{k}/(2*pi), ph_fit, 'Color', cmap(k,:), 'LineWidth', 1.6, 'HandleVisibility', 'off');
        end
    end
end
xlabel('Frequency [Hz]'); ylabel('Phase [rad]'); legend('Location', 'best'); grid on

%% Item 5 : comparaison des deformees - 4 figures (theorique vs identifiee)
for k = 1:Nmode
    phi_th = Phi(x_grid, k) / max(abs(Phi(x_grid, k)));
    figure('Name', sprintf('Mode shape - mode %d', k));
    plot(x_grid, phi_th, 'b-', 'LineWidth', 1.8, 'DisplayName', 'Modele (theorique)');
    hold on
    plot(xj, Phi_id(:,k), 'ro', 'MarkerSize', 8, 'LineWidth', 1.8, 'DisplayName', 'Identifie');
    xlabel('Beam length [m]'); ylabel(sprintf('Mode shape \\Phi_%d(x)', k));
    title(sprintf('Mode %d : f_{id} = %.2f Hz, \\xi_{id} = %.3f %%', k, f_id(k), 100*xi_id(k)));
    legend; grid on
end

%% Item 5 : tableau comparatif simule vs identifie
fprintf('\nTableau comparatif : parametres simules vs identifies\n');
fprintf('Mode |  f_sim [Hz]  |   f_id [Hz]  |  err f [%%] | xi_sim |   xi_id   | err xi [%%] |  err Phi\n');
fprintf('-----+--------------+--------------+------------+--------+-----------+------------+----------\n');
for k = 1:Nmode
    phi_th_xj = Phi(xj, k) / max(abs(Phi(xj, k)));
    err_phi   = norm(Phi_id(:,k) - phi_th_xj) / norm(phi_th_xj);
    fprintf('%4d | %12.4f | %12.4f | %9.4f | %6.4f | %9.5f | %9.4f | %.2e\n', ...
        k, nat_freq(k), f_id(k), 100*abs(f_id(k)-nat_freq(k))/nat_freq(k), ...
        ksi(k), xi_id(k), 100*abs(xi_id(k)-ksi(k))/ksi(k), err_phi);
end
fprintf('Erreur relative FRF dans chaque bande : '); fprintf('%.2e ', err_rel); fprintf('\n');

%  Fonction : modele FRF a un mode + residus complexes


function G_num = compute_FRF_num(x, w_s, nFRF, ind)
xi_i = x(1);
w_i  = x(2);
A    = x(ind(1,:));                          
RL   = x(ind(2,:)) + 1i * x(ind(3,:));        
RH   = x(ind(4,:)) + 1i * x(ind(5,:));       

S = numel(w_s);
G_num = zeros(S, nFRF);
den = -w_s.^2 + 1i*2*xi_i*w_i*w_s + w_i^2;     % denominateur commun (S x 1)

for r = 1:nFRF
    G_num(:,r) = A(r)./den + RL(r)./w_s.^2 + RH(r);
end
end

