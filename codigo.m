%% =========================================================================
%% SIMULAÇÃO COMPLETA - CONTROLE DE CRUZEIRO ADAPTATIVO (ACC)
%% Comparação: PID vs MPC
%% Autores: Lucas Catore, João Kawano, Ian Nascimento, Pedro Tiene, Eduardo Martinho
%% UTFPR - Apucarana
%% =========================================================================
clear; clc; close all;

%% =========================================================================
%% PARÂMETROS DO SISTEMA
%% =========================================================================
% Parâmetros do veículo (nominais)
m = 1500;           % Massa do veículo [kg]
Cd = 0.32;          % Coeficiente aerodinâmico
A = 2.2;            % Área frontal [m²]
rho = 1.225;        % Densidade do ar [kg/m³]
mu = 0.015;         % Coeficiente de atrito de rolamento
g = 9.81;           % Aceleração da gravidade [m/s²]

% Parâmetros de simulação
Ts = 0.1;           % Período de amostragem [s]
T_sim = 60;         % Tempo total de simulação [s]
t = 0:Ts:T_sim;     % Vetor de tempo
N = length(t);      % Número de amostras

% Referência de velocidade
v_ref = 22;         % Velocidade de referência [m/s] (≈ 79.2 km/h)

% Modelo linearizado: dv/dt = a*v + b*u
% Linearização em torno de v0 = 20 m/s
v0 = 20;
a_coef = -(rho * Cd * A * v0) / m;  % Coeficiente 'a' do modelo
b_coef = 1;                          % Coeficiente 'b' (normalizado)

%% =========================================================================
%% FUNÇÕES DE SIMULAÇÃO
%% =========================================================================

% Função para simular dinâmica do veículo
simulate_vehicle = @(v, u, m_sim, dt) ...
    v + dt * (u - (0.5*rho*Cd*A*v^2)/m_sim - mu*g);

%% =========================================================================
%% 1. ANÁLISE EM MALHA ABERTA
%% =========================================================================
fprintf('Simulando sistema em malha aberta...\n');

% Entrada constante (força de tração normalizada)
u_ol = 0.8 * ones(1, N);

% Simulação em malha aberta
v_ol = zeros(1, N);
v_ol(1) = 0;  % Velocidade inicial = 0

for k = 1:N-1
    v_ol(k+1) = simulate_vehicle(v_ol(k), u_ol(k), m, Ts);
end

% Figura 1: Resposta temporal em malha aberta
figure('Position', [100, 100, 800, 500]);
plot(t, v_ol, 'b-', 'LineWidth', 1.5);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Velocidade (m/s)', 'FontSize', 12);
title('Resposta temporal do sistema em malha aberta', 'FontSize', 14);
grid on;
saveas(gcf, 'figs/open_loop_time.png');

% Figura 2: Resposta em frequência (FFT)
figure('Position', [100, 100, 800, 500]);
Y = fft(v_ol);
f = (0:N-1)/(N*Ts);
plot(f(1:floor(N/2)), abs(Y(1:floor(N/2))), 'b-', 'LineWidth', 1.5);
xlabel('Frequência (Hz)', 'FontSize', 12);
ylabel('Magnitude', 'FontSize', 12);
title('Resposta em frequência (FFT) do sistema em malha aberta', 'FontSize', 14);
grid on;
xlim([0, 50]);
saveas(gcf, 'figs/open_loop_freq.png');

%% =========================================================================
%% 2. ANÁLISE DE ESTABILIDADE
%% =========================================================================
fprintf('Gerando diagramas de estabilidade...\n');

% Função de transferência do sistema
s = tf('s');
G = b_coef / (s - a_coef);

% Figura 3: Lugar das Raízes
figure('Position', [100, 100, 800, 500]);
rlocus(G);
title('Lugar das raízes - sistema 1a ordem', 'FontSize', 14);
xlabel('Re(s)', 'FontSize', 12);
ylabel('Im(s)', 'FontSize', 12);
grid on;
saveas(gcf, 'figs/root_locus.png');

% Figura 4: Diagrama de Bode
figure('Position', [100, 100, 800, 600]);
bode(G);
title('Diagrama de Bode', 'FontSize', 14);
grid on;
saveas(gcf, 'figs/bode.png');

% Figura 5: Diagrama de Nyquist
figure('Position', [100, 100, 800, 500]);
nyquist(G);
title('Diagrama de Nyquist', 'FontSize', 14);
grid on;
saveas(gcf, 'figs/nyquist.png');

%% =========================================================================
%% 3. CONTROLADOR PID
%% =========================================================================
fprintf('Simulando controlador PID...\n');

% Ganhos do PID
Kp = 0.8;
Ki = 0.05;
Kd = 0.1;

% Simulação PID - Caso Nominal
v_pid = zeros(1, N);
u_pid = zeros(1, N);
e_pid = zeros(1, N);
e_int = 0;
e_prev = 0;

for k = 1:N-1
    % Erro
    e_pid(k) = v_ref - v_pid(k);
    
    % Termos PID
    e_int = e_int + e_pid(k) * Ts;
    e_deriv = (e_pid(k) - e_prev) / Ts;
    
    % Ação de controle
    u_pid(k) = Kp * e_pid(k) + Ki * e_int + Kd * e_deriv;
    
    % Dinâmica do veículo
    v_pid(k+1) = simulate_vehicle(v_pid(k), u_pid(k), m, Ts);
    
    e_prev = e_pid(k);
end
u_pid(N) = u_pid(N-1);
e_pid(N) = v_ref - v_pid(N);

% Figura 6: Resposta PID em malha fechada (ANTIGA - mantida para compatibilidade)
figure('Position', [100, 100, 800, 500]);
plot(t, v_pid, 'Color', [0.85, 0.33, 0.1], 'LineWidth', 1.5);
hold on;
plot(t, v_ref*ones(size(t)), 'b--', 'LineWidth', 1.5);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Velocidade (m/s)', 'FontSize', 12);
title('Resposta do sistema com controle PID em malha fechada e referência', 'FontSize', 14);
legend('Saída (PID)', 'Referência', 'Location', 'best');
grid on;
ylim([0, 25]);
saveas(gcf, 'figs/pid_closedloop.png');

% Figura 7: Sinal de controle PID (ANTIGA - mantida para compatibilidade)
figure('Position', [100, 100, 800, 500]);
plot(t, u_pid, 'b-', 'LineWidth', 1.5);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Sinal de controle u(t)', 'FontSize', 12);
title('Sinal de controle gerado pelo PID', 'FontSize', 14);
grid on;
saveas(gcf, 'figs/pid_control_signal.png');

%% =========================================================================
%% NOVO: Figura PID Combinada (Subplot) - Resposta + Controle
%% =========================================================================
figure('Position', [100, 100, 800, 700]);

% Subplot 1: Velocidade
subplot(2,1,1);
plot(t, v_pid, 'Color', [0.85, 0.33, 0.1], 'LineWidth', 1.5);
hold on;
plot(t, v_ref*ones(size(t)), 'b--', 'LineWidth', 1.5);
xlabel('Tempo (s)', 'FontSize', 11);
ylabel('Velocidade (m/s)', 'FontSize', 11);
title('Resposta do Sistema com Controlador PID', 'FontSize', 12);
legend('Saída (PID)', 'Referência', 'Location', 'best');
grid on;
ylim([0, 25]);

% Subplot 2: Sinal de Controle
subplot(2,1,2);
plot(t, u_pid, 'b-', 'LineWidth', 1.5);
hold on;
% Linhas de limite (para visualização)
plot(t, 2*ones(size(t)), 'k--', 'LineWidth', 1);
plot(t, -3*ones(size(t)), 'k--', 'LineWidth', 1);
xlabel('Tempo (s)', 'FontSize', 11);
ylabel('Aceleração (m/s^2)', 'FontSize', 11);
title('Sinal de Controle - PID', 'FontSize', 12);
legend('u(t) PID', 'Limites físicos', 'Location', 'best');
grid on;
ylim([-5, 20]);

saveas(gcf, 'figs/pid_combined.png');

%% =========================================================================
%% 4. PID COM VARIAÇÃO PARAMÉTRICA (+10% massa)
%% =========================================================================
fprintf('Simulando PID com variação de massa...\n');

m_var = m * 1.10;  % Massa +10%

v_pid_var = zeros(1, N);
u_pid_var = zeros(1, N);
e_int_var = 0;
e_prev_var = 0;

for k = 1:N-1
    e_var = v_ref - v_pid_var(k);
    e_int_var = e_int_var + e_var * Ts;
    e_deriv_var = (e_var - e_prev_var) / Ts;
    
    u_pid_var(k) = Kp * e_var + Ki * e_int_var + Kd * e_deriv_var;
    v_pid_var(k+1) = simulate_vehicle(v_pid_var(k), u_pid_var(k), m_var, Ts);
    
    e_prev_var = e_var;
end
u_pid_var(N) = u_pid_var(N-1);

% Figura 8: Comparação PID nominal vs +10% massa
figure('Position', [100, 100, 800, 500]);
plot(t, v_pid, 'b-', 'LineWidth', 1.5);
hold on;
plot(t, v_pid_var, 'r--', 'LineWidth', 1.5);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Velocidade (m/s)', 'FontSize', 12);
title('Comportamento do sistema com variação de parâmetros', 'FontSize', 14);
legend('PID (nominal)', 'PID (+10% massa)', 'Location', 'best');
grid on;
ylim([0, 25]);
saveas(gcf, 'figs/pid_variation.png');

%% =========================================================================
%% 5. CONTROLADOR MPC
%% =========================================================================
fprintf('Simulando controlador MPC...\n');

% Parâmetros do MPC
Np = 20;           % Horizonte de predição
Nc = 5;            % Horizonte de controle
Q = 10;            % Peso do erro
R = 0.1;           % Peso da ação de controle
u_min = -3;        % Limite inferior de aceleração [m/s²]
u_max = 2;         % Limite superior de aceleração [m/s²]

% Simulação MPC - Caso Nominal
v_mpc = zeros(1, N);
u_mpc = zeros(1, N);

for k = 1:N-1
    % Erro atual
    e_mpc = v_ref - v_mpc(k);
    
    % Cálculo simplificado do MPC (otimização quadrática)
    % Para sistema de primeira ordem: u_opt = saturar(K * e)
    K_mpc = Q / (Q + R);
    u_opt = K_mpc * e_mpc;
    
    % Aplicar restrições de saturação
    u_mpc(k) = max(u_min, min(u_max, u_opt));
    
    % Dinâmica do veículo
    v_mpc(k+1) = simulate_vehicle(v_mpc(k), u_mpc(k), m, Ts);
end
u_mpc(N) = u_mpc(N-1);

% Figura 9: Resposta MPC em malha fechada (ANTIGA - mantida para compatibilidade)
figure('Position', [100, 100, 800, 500]);
plot(t, v_mpc, 'Color', [0.85, 0.33, 0.1], 'LineWidth', 1.5);
hold on;
plot(t, v_ref*ones(size(t)), 'b--', 'LineWidth', 1.5);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Velocidade (m/s)', 'FontSize', 12);
title('Resposta do sistema com controlador MPC', 'FontSize', 14);
legend('Saída (MPC)', 'Referência', 'Location', 'best');
grid on;
ylim([0, 25]);
saveas(gcf, 'figs/mpc_closedloop.png');

% Figura 10: Sinal de controle MPC (ANTIGA - mantida para compatibilidade)
figure('Position', [100, 100, 800, 500]);
plot(t, u_mpc, 'b-', 'LineWidth', 1.5);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Sinal de controle u(t)', 'FontSize', 12);
title('Sinal de controle do método MPC', 'FontSize', 14);
grid on;
ylim([0, 2.5]);
saveas(gcf, 'figs/smc_control_signal.png');

%% =========================================================================
%% NOVO: Figura MPC Combinada (Subplot) - Resposta + Controle
%% =========================================================================
figure('Position', [100, 100, 800, 700]);

% Subplot 1: Velocidade
subplot(2,1,1);
plot(t, v_mpc, 'Color', [0.85, 0.33, 0.1], 'LineWidth', 1.5);
hold on;
plot(t, v_ref*ones(size(t)), 'b--', 'LineWidth', 1.5);
xlabel('Tempo (s)', 'FontSize', 11);
ylabel('Velocidade (m/s)', 'FontSize', 11);
title('Resposta do Sistema com Controlador MPC', 'FontSize', 12);
legend('Saída (MPC)', 'Referência', 'Location', 'best');
grid on;
ylim([0, 25]);

% Subplot 2: Sinal de Controle com limites
subplot(2,1,2);
plot(t, u_mpc, 'r-', 'LineWidth', 1.5);
hold on;
plot(t, u_max*ones(size(t)), 'k--', 'LineWidth', 1);
plot(t, u_min*ones(size(t)), 'k--', 'LineWidth', 1);
xlabel('Tempo (s)', 'FontSize', 11);
ylabel('Aceleração (m/s^2)', 'FontSize', 11);
title('Sinal de Controle - MPC (com limites de saturação)', 'FontSize', 12);
legend('u(t) MPC', 'Limites [-3, 2] m/s^2', 'Location', 'best');
grid on;
ylim([-4, 3]);

saveas(gcf, 'figs/mpc_combined.png');

%% =========================================================================
%% NOVO: MPC COM VARIAÇÃO PARAMÉTRICA (+10% massa)
%% =========================================================================
fprintf('Simulando MPC com variação de massa...\n');

% Simulação MPC com massa +10%
v_mpc_var = zeros(1, N);
u_mpc_var = zeros(1, N);

for k = 1:N-1
    e_mpc_var = v_ref - v_mpc_var(k);
    K_mpc = Q / (Q + R);
    u_opt_var = K_mpc * e_mpc_var;
    u_mpc_var(k) = max(u_min, min(u_max, u_opt_var));
    v_mpc_var(k+1) = simulate_vehicle(v_mpc_var(k), u_mpc_var(k), m_var, Ts);
end
u_mpc_var(N) = u_mpc_var(N-1);

% NOVA Figura: Comparação MPC nominal vs +10% massa
figure('Position', [100, 100, 800, 500]);
plot(t, v_mpc, 'b-', 'LineWidth', 1.5);
hold on;
plot(t, v_mpc_var, 'r--', 'LineWidth', 1.5);
plot(t, v_ref*ones(size(t)), 'k:', 'LineWidth', 1);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Velocidade (m/s)', 'FontSize', 12);
title('Robustez do MPC: variação de +10% na massa', 'FontSize', 14);
legend('MPC (nominal)', 'MPC (+10% massa)', 'Referência', 'Location', 'best');
grid on;
ylim([0, 25]);
saveas(gcf, 'figs/mpc_variation.png');

%% =========================================================================
%% 6. COMPARAÇÃO GERAL PID vs MPC
%% =========================================================================
fprintf('Gerando comparação PID vs MPC...\n');

% Figura 11: Comparação completa
figure('Position', [100, 100, 800, 700]);

% Subplot 1: Respostas
subplot(2,1,1);
plot(t, v_ref*ones(size(t)), 'k--', 'LineWidth', 1.5);
hold on;
plot(t, v_pid, 'b-', 'LineWidth', 1.5);
plot(t, v_mpc, 'r-', 'LineWidth', 1.5);
xlabel('Tempo (s)', 'FontSize', 11);
ylabel('Velocidade (m/s)', 'FontSize', 11);
title('Comparação: Resposta dos Controladores', 'FontSize', 12);
legend('Referência', 'PID', 'MPC', 'Location', 'best');
grid on;
ylim([0, 25]);

% Subplot 2: Sinais de controle
subplot(2,1,2);
plot(t, u_pid, 'b-', 'LineWidth', 1.5);
hold on;
plot(t, u_mpc, 'r-', 'LineWidth', 1.5);
plot(t, u_max*ones(size(t)), 'k--', 'LineWidth', 1);
plot(t, u_min*ones(size(t)), 'k--', 'LineWidth', 1);
xlabel('Tempo (s)', 'FontSize', 11);
ylabel('Sinal de controle (m/s^2)', 'FontSize', 11);
title('Comparação: Sinais de Controle', 'FontSize', 12);
legend('PID', 'MPC', 'Limites', 'Location', 'best');
grid on;
ylim([-5, 20]);

saveas(gcf, 'figs/comparacao_controladores.png');

%% =========================================================================
%% 7. MÉTRICAS DE DESEMPENHO
%% =========================================================================
fprintf('Calculando métricas de desempenho...\n');

% Erro PID
e_pid_full = v_ref - v_pid;
% Erro MPC
e_mpc_full = v_ref - v_mpc;

% RMSE
RMSE_pid = sqrt(mean(e_pid_full.^2));
RMSE_mpc = sqrt(mean(e_mpc_full.^2));

% IAE (Integral Absolute Error)
IAE_pid = sum(abs(e_pid_full)) * Ts;
IAE_mpc = sum(abs(e_mpc_full)) * Ts;

% ISE (Integral Squared Error)
ISE_pid = sum(e_pid_full.^2) * Ts;
ISE_mpc = sum(e_mpc_full.^2) * Ts;

% ITAE (Integral Time-weighted Absolute Error)
ITAE_pid = sum(t .* abs(e_pid_full)) * Ts;
ITAE_mpc = sum(t .* abs(e_mpc_full)) * Ts;

% Exibir métricas
fprintf('\n========== MÉTRICAS DE DESEMPENHO ==========\n');
fprintf('Métrica    |    PID    |    MPC    | Razão\n');
fprintf('----------------------------------------------\n');
fprintf('RMSE       | %8.4f  | %8.4f  | %.1fx\n', RMSE_pid, RMSE_mpc, RMSE_mpc/RMSE_pid);
fprintf('IAE        | %8.4f  | %8.4f  | %.1fx\n', IAE_pid, IAE_mpc, IAE_mpc/IAE_pid);
fprintf('ISE        | %8.4f  | %8.4f  | %.1fx\n', ISE_pid, ISE_mpc, ISE_mpc/ISE_pid);
fprintf('ITAE       | %8.4f  | %8.4f  | %.1fx\n', ITAE_pid, ITAE_mpc, ITAE_mpc/ITAE_pid);
fprintf('==============================================\n');

% Figura: Tabela de métricas
figure('Position', [100, 100, 600, 400]);
axis off;
title('Métricas de erro para comparação dos controladores', 'FontSize', 14);
text(0.5, 0.7, sprintf('Métrica  | PID      | MPC'), 'FontSize', 12, 'FontName', 'Courier', 'HorizontalAlignment', 'center');
text(0.5, 0.6, sprintf('RMSE     | %.4f   | %.4f', RMSE_pid, RMSE_mpc), 'FontSize', 12, 'FontName', 'Courier', 'HorizontalAlignment', 'center');
text(0.5, 0.5, sprintf('IAE      | %.4f  | %.4f', IAE_pid, IAE_mpc), 'FontSize', 12, 'FontName', 'Courier', 'HorizontalAlignment', 'center');
text(0.5, 0.4, sprintf('ISE      | %.4f | %.4f', ISE_pid, ISE_mpc), 'FontSize', 12, 'FontName', 'Courier', 'HorizontalAlignment', 'center');
text(0.5, 0.3, sprintf('ITAE     | %.4f | %.4f', ITAE_pid, ITAE_mpc), 'FontSize', 12, 'FontName', 'Courier', 'HorizontalAlignment', 'center');
saveas(gcf, 'figs/metrics_table.png');

%% =========================================================================
%% 8. ANÁLISE ESPECTRAL (ESPECTROGRAMAS)
%% =========================================================================
fprintf('Gerando espectrogramas...\n');

% Parâmetros STFT
window_size = 64;
overlap = 56;
nfft = 128;
fs = 1/Ts;

% Espectrograma PID
figure('Position', [100, 100, 800, 500]);
spectrogram(u_pid, window_size, overlap, nfft, fs, 'yaxis');
title('Espectrograma do sinal de controle para o PID (custom)', 'FontSize', 14);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Freq (Hz)', 'FontSize', 12);
colorbar;
saveas(gcf, 'figs/spectrogram_pid.png');

% Espectrograma MPC
figure('Position', [100, 100, 800, 500]);
spectrogram(u_mpc, window_size, overlap, nfft, fs, 'yaxis');
title('Espectrograma do sinal de controle para o MPC (custom)', 'FontSize', 14);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Freq (Hz)', 'FontSize', 12);
colorbar;
saveas(gcf, 'figs/spectrogram_smc.png');

%% =========================================================================
%% FIM DA SIMULAÇÃO
%% =========================================================================
fprintf('\n========================================\n');
fprintf('Simulação concluída com sucesso!\n');
fprintf('Todas as figuras salvas na pasta "figs/"\n');
fprintf('========================================\n');

% Lista de arquivos gerados
fprintf('\nArquivos gerados:\n');
fprintf('  - open_loop_time.png\n');
fprintf('  - open_loop_freq.png\n');
fprintf('  - root_locus.png\n');
fprintf('  - bode.png\n');
fprintf('  - nyquist.png\n');
fprintf('  - pid_closedloop.png\n');
fprintf('  - pid_control_signal.png\n');
fprintf('  - pid_combined.png       [NOVO - Subplot PID]\n');
fprintf('  - pid_variation.png\n');
fprintf('  - mpc_closedloop.png\n');
fprintf('  - smc_control_signal.png\n');
fprintf('  - mpc_combined.png       [NOVO - Subplot MPC]\n');
fprintf('  - mpc_variation.png      [NOVO - Robustez MPC]\n');
fprintf('  - comparacao_controladores.png\n');
fprintf('  - metrics_table.png\n');
fprintf('  - spectrogram_pid.png\n');
fprintf('  - spectrogram_smc.png\n');
