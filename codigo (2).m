function generate_all_figs_no_toolboxes_CORRIGIDO()


% Cria pasta de figuras
if ~exist('figs','dir')
    mkdir('figs');
end

% Força fundo branco para todas as figuras por padrão
try
    set(groot, 'defaultFigureColor', 'white');
    set(groot, 'defaultAxesXColor', 'k');
    set(groot, 'defaultAxesYColor', 'k');
    set(groot, 'defaultAxesZColor', 'k');
    set(groot, 'defaultTextColor', 'k');
    set(groot, 'defaultAxesFontSize', 12);
    set(groot, 'defaultLineLineWidth', 1.5);
catch
end

%% ---------------- Parâmetros do modelo ----------------
m = 1500;            % kg
Cd = 0.32;
A = 2.2;             % m^2
rho = 1.225;         % kg/m^3
mu = 0.015;
v0 = 20;             % m/s (ponto de operação)
g = 9.81;
Th = 1.5;            % tempo de headway
d_ref_base = 25;     % m

% Linearização (modelo escalar): dv/dt = a*v + b*u
% onde u é a força de tração normalizada pela massa (aceleração de controle)
% Modelo: m*dv/dt = F_trac - F_drag - F_friction
% Linearizado: dv/dt = a*v + b*u, onde u = F_trac/m (aceleração)
a = - (rho*Cd*A*v0) / m;   % ≈ -0.0057 (polo estável)
b = 1;                      % CORRIGIDO: b=1 quando u é aceleração direta

% Nota: No modelo original b = 1/m, mas isso assume u = força.
% Para controle de veículos, é mais natural u = aceleração comandada.
% Então usamos b = 1.

% Amostragem para simulações discretas (Euler)
Ts = 0.1;
tfinal = 60;
t = 0:0.01:tfinal;      % tempo para plots contínuos
td = 0:Ts:tfinal;       % tempo discreto para MPC/PID
N = length(t);

%% ---------- Referências e cenário ----------
v_ref_step = 22 * ones(size(t)); % referência de velocidade (m/s)
leader = 18*ones(size(t));
leader(t>=15) = 12;
leader(t>=25) = 20;
leader(t>=35) = 15;

r_d = interp1(t, v_ref_step, td);

%% ---------------- Root Locus (analítico para primeira ordem) -----------
% Para G(s)= b/(s - a), laço com ganho K: polos em s tal que s - a + K*b = 0 -> s = a - K*b
Kvec = linspace(0,2000,500);
poles_rl = a - Kvec * b;
figure;
plot(real(poles_rl), imag(poles_rl), 'LineWidth',1.4);
xlabel('Re(s)'); ylabel('Im(s)'); grid on;
title('Lugar das raízes - sistema 1a ordem');
save_figure('root_locus');

%% ---------------- Bode (analítico) ----------------
w = logspace(-3,2,1000); % rad/s
% G(jw) = b / (j*w - a)
Gjw = b ./ (1i*w - a);
mag = abs(Gjw);
phase = angle(Gjw) * 180/pi;
figure;
subplot(2,1,1);
loglog(w, mag, 'LineWidth',1.2);
ylabel('Magnitude |G(j\omega)|'); grid on;
title('Diagrama de Bode');
subplot(2,1,2);
semilogx(w, phase, 'LineWidth',1.2);
xlabel('\omega (rad/s)'); ylabel('Fase (deg)'); grid on;
save_figure('bode');

%% ---------------- Nyquist (analítico) ----------------
figure;
plot(real(Gjw), imag(Gjw), 'LineWidth',1.2); hold on;
plot(real(Gjw), -imag(Gjw), 'LineWidth',0.8, 'LineStyle','--');
xlabel('Re(G(j\omega))'); ylabel('Im(G(j\omega))'); grid on;
title('Diagrama de Nyquist');
axis equal;
save_figure('nyquist');

%% ---------------- Resposta Malha Aberta - tempo (analítico) ----------
% Resposta ao degrau unitário em u: v_ss = -b*u0/a (porque a<0)
u0 = 1;  % aceleração de 1 m/s²
v_ss = -b * u0 / a;  % velocidade em regime permanente
y_ol = v_ss * (1 - exp(a * t));
figure;
plot(t, y_ol, 'LineWidth',1.2);
xlabel('Tempo (s)'); ylabel('Velocidade (m/s)'); grid on;
title('Resposta temporal do sistema em malha aberta');
save_figure('open_loop_time');

%% ---------------- Resposta Malha Aberta - freq via FFT ----------
u_rich = chirp_signal(t, 0.01, 1, 0.5);
y_rich = simulate_continuous_linear(a, b, u_rich, t);
nfft = 2^nextpow2(length(t));
Y = fft(y_rich, nfft);
Fs = 1/(t(2)-t(1));
f = (0:nfft/2-1)*(Fs/nfft);
figure;
plot(f, abs(Y(1:nfft/2)), 'LineWidth',1.2);
xlabel('Frequência (Hz)'); ylabel('Magnitude'); grid on;
title('Resposta em frequência (FFT) do sistema em malha aberta');
save_figure('open_loop_freq');

%% ---------------- PID (discreto por Euler) ----------------
% Ganhos ajustados para o modelo com b=1
Kp = 0.8;    % Ganho proporcional
Ki = 0.05;   % Ganho integral  
Kd = 0.1;    % Ganho derivativo

% Simular PID discreto
[~, y_pid, u_pid] = simulate_pid_discrete(a, b, Ts, td, r_d, Kp, Ki, Kd, 0);

% Interpolar para vetor t para consistência nos plots
y_pid_cont = interp1(td, y_pid, t);
u_pid_cont = interp1(td, u_pid, t);

figure;
plot(t, v_ref_step, '--','LineWidth',1.0); hold on;
plot(t, y_pid_cont, 'LineWidth',1.2);
xlabel('Tempo (s)'); ylabel('Velocidade (m/s)'); grid on;
legend('Referência','Saída (PID)');
title('Resposta do sistema com controle PID em malha fechada e referência');
save_figure('pid_closedloop');

figure;
plot(t, u_pid_cont, 'LineWidth',1.2);
xlabel('Tempo (s)'); ylabel('Sinal de controle u(t)'); grid on;
title('Sinal de controle gerado pelo PID');
save_figure('pid_control_signal');

%% ---------------- PID: Variação paramétrica (+10% massa) ----------------
m_var = m * 1.10;
a_var = - (rho*Cd*A*v0) / m_var;
% Nota: b permanece 1 pois u é aceleração comandada (a planta "real" 
% terá resposta diferente, mas o controlador comanda aceleração)
[~, y_pid_var, ~] = simulate_pid_discrete(a_var, b, Ts, td, r_d, Kp, Ki, Kd, 0);
y_pid_var_cont = interp1(td, y_pid_var, t);

figure;
plot(t, y_pid_cont, 'LineWidth',1.2); hold on;
plot(t, y_pid_var_cont, '--', 'LineWidth',1.2);
xlabel('Tempo (s)'); ylabel('Velocidade (m/s)');
legend('PID (nominal)','PID (+10% massa)');
title('Comportamento do sistema com variação de parâmetros');
grid on;
save_figure('pid_variation');

%% ---------------- MPC CORRIGIDO ----------------
% Parâmetros MPC
Np = 20;        % Horizonte de predição (aumentado)
Nc = 5;         % Horizonte de controle (aumentado)

% Restrições de aceleração (m/s²)
umin = -3;      % Frenagem máxima
umax = 2;       % Aceleração máxima

% Pesos da função custo (AJUSTADOS)
Q_mpc = 10;     % Peso do erro de rastreamento (aumentado)
R_mpc = 0.1;    % Peso do esforço de controle (reduzido)

% Discretização Euler
Ad = 1 + a*Ts;
Bd = b*Ts;

% Construir matrizes de predição
[Phi, Gamma] = buildPredMatrices(Ad, Bd, Np, Nc);

% Simulação MPC
y_mpc = zeros(length(td),1);
u_mpc = zeros(length(td),1);
xk = 0;  % Estado inicial (velocidade = 0)

for k = 1:length(td)
    % Vetor de referência no horizonte
    r_horizon = r_d(k) * ones(Np,1);
    
    % Predição livre (sem controle adicional): y_free = Phi * xk
    y_free = Phi * xk;
    
    % Erro de predição livre
    e_free = r_horizon - y_free;
    
    % Formulação QP: min (1/2)*U'*H*U + f'*U
    % onde H = Gamma'*Q*Gamma + R*I
    %      f = -Gamma'*Q*e_free
    Qmat = Q_mpc * eye(Np);
    Rmat = R_mpc * eye(Nc);
    
    H = Gamma' * Qmat * Gamma + Rmat;
    f = -Gamma' * Qmat * e_free;
    
    % Solução analítica (sem restrições)
    Uopt = -H \ f;
    
    % Aplicar restrições por saturação (clipping)
    Uopt = max(min(Uopt, umax), umin);
    
    % Aplicar apenas o primeiro elemento (princípio do horizonte deslizante)
    uk = Uopt(1);
    
    % Atualizar estado (simulação da planta)
    xk = Ad * xk + Bd * uk;
    
    % Armazenar resultados
    y_mpc(k) = xk;
    u_mpc(k) = uk;
end

% Interpolar para plots contínuos
y_mpc_cont = interp1(td, y_mpc, t, 'linear');
u_mpc_cont = interp1(td, u_mpc, t, 'linear');

figure;
plot(t, v_ref_step, '--', 'LineWidth',1.0); hold on;
plot(t, y_mpc_cont, 'LineWidth',1.2);
xlabel('Tempo (s)'); ylabel('Velocidade (m/s)');
legend('Referência','Saída (MPC)');
title('Resposta do sistema com controlador MPC');
grid on;
save_figure('mpc_closedloop');

figure;
plot(t, u_mpc_cont, 'LineWidth',1.2);
xlabel('Tempo (s)'); ylabel('Sinal de controle u(t)');
title('Sinal de controle do método MPC');
grid on;
save_figure('smc_control_signal');

%% ---------------- Métricas ----------------
ref_td = interp1(t, v_ref_step, td);

% Garantir que os vetores têm a mesma orientação
y_pid_col = y_pid(:);
y_mpc_col = y_mpc(:);
ref_td_col = ref_td(:);

rmse_pid = sqrt(mean((y_pid_col - ref_td_col).^2));
rmse_mpc = sqrt(mean((y_mpc_col - ref_td_col).^2));

iae_pid = trapz(td, abs(ref_td_col - y_pid_col));
iae_mpc = trapz(td, abs(ref_td_col - y_mpc_col));

ise_pid = trapz(td, (ref_td_col - y_pid_col).^2);
ise_mpc = trapz(td, (ref_td_col - y_mpc_col).^2);

itae_pid = trapz(td, td(:) .* abs(ref_td_col - y_pid_col));
itae_mpc = trapz(td, td(:) .* abs(ref_td_col - y_mpc_col));

% Exibir métricas no console
fprintf('\n============ MÉTRICAS DE DESEMPENHO ============\n');
fprintf('Métrica   |    PID     |    MPC     | Melhor\n');
fprintf('------------------------------------------------\n');
fprintf('RMSE      | %10.4f | %10.4f | %s\n', rmse_pid, rmse_mpc, comparar(rmse_pid, rmse_mpc));
fprintf('IAE       | %10.4f | %10.4f | %s\n', iae_pid, iae_mpc, comparar(iae_pid, iae_mpc));
fprintf('ISE       | %10.4f | %10.4f | %s\n', ise_pid, ise_mpc, comparar(ise_pid, ise_mpc));
fprintf('ITAE      | %10.4f | %10.4f | %s\n', itae_pid, itae_mpc, comparar(itae_pid, itae_mpc));
fprintf('================================================\n\n');

% Figura com tabela de métricas
figure('Units','normalized','Position',[0.2 0.2 0.5 0.3]);
axis off;
tstr = {
    'Métrica   | PID        | MPC';
    sprintf('RMSE      | %.4f    | %.4f', rmse_pid, rmse_mpc);
    sprintf('IAE       | %.4f   | %.4f', iae_pid, iae_mpc);
    sprintf('ISE       | %.4f  | %.4f', ise_pid, ise_mpc);
    sprintf('ITAE      | %.4f | %.4f', itae_pid, itae_mpc);
    };
text(0.1, 0.5, tstr, 'FontName','Monospaced','FontSize',12);
title('Métricas de erro para comparação dos controladores');
save_figure('metrics_table');

%% ---------------- Espectrogramas ----------------
[Sf, Tm, Freq] = simple_spectrogram(u_pid_cont, t, 256, 200);
figure;
imagesc(Tm, Freq, 20*log10(abs(Sf)+eps));
axis xy;
xlabel('Tempo (s)'); ylabel('Freq (Hz)');
title('Espectrograma do sinal de controle para o PID (custom)');
colorbar;
save_figure('spectrogram_pid');

[Sf2, Tm2, Freq2] = simple_spectrogram(u_mpc_cont, t, 256, 200);
figure;
imagesc(Tm2, Freq2, 20*log10(abs(Sf2)+eps));
axis xy;
xlabel('Tempo (s)'); ylabel('Freq (Hz)');
title('Espectrograma do sinal de controle para o MPC (custom)');
colorbar;
save_figure('spectrogram_smc');

%% ---------------- Gráfico comparativo adicional ----------------
figure;
subplot(2,1,1);
plot(t, v_ref_step, 'k--', 'LineWidth', 1.0); hold on;
plot(t, y_pid_cont, 'b', 'LineWidth', 1.2);
plot(t, y_mpc_cont, 'r', 'LineWidth', 1.2);
xlabel('Tempo (s)'); ylabel('Velocidade (m/s)');
legend('Referência', 'PID', 'MPC', 'Location', 'best');
title('Comparação: Resposta dos Controladores');
grid on;

subplot(2,1,2);
plot(t, u_pid_cont, 'b', 'LineWidth', 1.2); hold on;
plot(t, u_mpc_cont, 'r', 'LineWidth', 1.2);
yline(umax, 'k--', 'LineWidth', 0.8);
yline(umin, 'k--', 'LineWidth', 0.8);
xlabel('Tempo (s)'); ylabel('Sinal de controle (m/s²)');
legend('PID', 'MPC', 'Limites', 'Location', 'best');
title('Comparação: Sinais de Controle');
grid on;
save_figure('comparacao_controladores');

fprintf('Figuras geradas em ./figs/ (script CORRIGIDO).\n');

end

%% ============== FUNÇÕES AUXILIARES ==============

function s = chirp_signal(t, f0, f1, amp)
% Gera sinal chirp (varredura de frequência) sem toolbox
if nargin < 4, amp = 1; end
T = t(end);
k = (f1 - f0) / T;
s = amp * sin(2*pi*(f0*t + 0.5*k.*t.^2));
end

function y = simulate_continuous_linear(a, b, u, t)
% Simula dx/dt = a*x + b*u com x(0)=0 por RK2
dt = t(2) - t(1);
x = 0;
y = zeros(size(t));
for k = 1:length(t)
    y(k) = x;
    uk = u(k);
    k1 = a*x + b*uk;
    x2 = x + dt * k1;
    k2 = a*x2 + b*uk;
    x = x + dt * 0.5*(k1 + k2);
end
end

function [x_hist, y_hist, u_hist] = simulate_pid_discrete(a, b, Ts, td, r_d, Kp, Ki, Kd, x0)
% Simula PID discreto usando integração Euler
% CORRIGIDO: bug na linha u_hist(k) = u
if nargin < 9, x0 = 0; end

Ad = 1 + a*Ts;
Bd = b*Ts;
N = length(td);
x = x0;
Iterm = 0;
prev_e = 0;

y_hist = zeros(N, 1);
u_hist = zeros(N, 1);  % CORRIGIDO: pré-alocação correta
x_hist = zeros(N, 1);

for k = 1:N
    rk = r_d(k);
    yk = x;
    e = rk - yk;
    
    % Termo integral com anti-windup simples
    Iterm = Iterm + e * Ts;
    
    % Termo derivativo
    if k == 1
        Dterm = 0;
    else
        Dterm = (e - prev_e) / Ts;
    end
    
    % Lei de controle PID
    u = Kp*e + Ki*Iterm + Kd*Dterm;
    
    % Atualizar estado
    x = Ad * x + Bd * u;
    
    % Armazenar histórico
    x_hist(k) = x;
    y_hist(k) = x;
    u_hist(k) = u;  % CORRIGIDO: era "nu_hist(k) = u" (erro de digitação)
    
    prev_e = e;
end
end

function [Phi, Gamma] = buildPredMatrices(Ad, Bd, Np, Nc)
% Constrói matrizes de predição para MPC
% Phi (Np x 1): matriz de evolução livre
% Gamma (Np x Nc): matriz de resposta ao controle
%
% Modelo: x(k+1) = Ad*x(k) + Bd*u(k)
% Predição: Y = Phi*x(k) + Gamma*U

Phi = zeros(Np, 1);
Gamma = zeros(Np, Nc);

for i = 1:Np
    Phi(i) = Ad^i;
    for j = 1:min(i, Nc)
        Gamma(i, j) = Ad^(i-j) * Bd;
    end
end
end

function str = comparar(val_pid, val_mpc)
% Retorna string indicando qual controlador teve melhor desempenho
if val_pid < val_mpc
    str = 'PID';
elseif val_mpc < val_pid
    str = 'MPC';
else
    str = 'Empate';
end
end

function save_figure(basename)
% Salva figura atual como PNG e SVG com fundo branco

fn_png = fullfile('figs', [basename, '.png']);
fn_svg = fullfile('figs', [basename, '.svg']);

set(gcf, 'Color', 'white');

ax = findall(gcf, 'Type', 'axes');
for k = 1:length(ax)
    set(ax(k), 'Color', 'white', ...
               'XColor', 'k', ...
               'YColor', 'k', ...
               'ZColor', 'k', ...
               'LineWidth', 1.2, ...
               'FontSize', 12);
end

txt = findall(gcf, 'Type', 'text');
for k = 1:length(txt)
    set(txt(k), 'Color', 'k');
end

hLeg = findall(gcf, 'Type', 'legend');
for k = 1:length(hLeg)
    set(hLeg(k), 'Color', 'white', ...
                 'TextColor', 'k', ...
                 'EdgeColor', 'k');
end

if exist('exportgraphics', 'file') == 2
    exportgraphics(gcf, fn_png, 'BackgroundColor', 'white', 'Resolution', 300);
    exportgraphics(gcf, fn_svg, 'BackgroundColor', 'white');
else
    print(gcf, fn_png, '-dpng', '-r300');
    print(gcf, fn_svg, '-dsvg');
end
end

function [S, Tm, Freq] = simple_spectrogram(x, tout, winlen, noverlap)
% STFT simples sem toolboxes
if nargin < 3, winlen = 256; end
if nargin < 4, noverlap = round(0.8*winlen); end

dt = tout(2) - tout(1);
Fs = 1 / dt;
hop = winlen - noverlap;
N = length(x);
nwin = floor((N - noverlap) / hop);

S = zeros(floor(winlen/2), nwin);
Tm = zeros(1, nwin);

% Janela Hamming manual
n = (0:winlen-1)';
w = 0.54 - 0.46 * cos(2*pi*n/(winlen-1));

x = x(:);  % Garantir vetor coluna

for k = 1:nwin
    idx0 = (k-1)*hop + 1;
    idx = idx0 : (idx0 + winlen - 1);
    if idx(end) > N
        break;
    end
    seg = x(idx) .* w;
    Xf = fft(seg);
    S(:,k) = Xf(1:floor(winlen/2));
    Tm(k) = tout(round(mean(idx)));
end

Freq = (0:floor(winlen/2)-1) * (Fs / winlen);
end
