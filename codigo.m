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
a_coef = -(rho * Cd * A * v0) / m;  % Coeficiente 'a' do modelo (Polo)
b_coef = 1;                         % Coeficiente 'b' (normalizado)

% Cria pasta para figuras se não existir
if ~exist('figs', 'dir')
    mkdir('figs');
end

%% =========================================================================
%% FUNÇÕES DE SIMULAÇÃO
%% =========================================================================
% Função anônima para simular dinâmica longitudinal do veículo
simulate_vehicle = @(v, u, m_sim, dt) ...
    v + dt * (u - (0.5*rho*Cd*A*v^2)/m_sim - mu*g);

%% =========================================================================
%% 1. ANÁLISE EM MALHA ABERTA
%% =========================================================================
fprintf('1. Simulando sistema em malha aberta...\n');

% Entrada constante (degrau de força de tração normalizada)
u_ol = 0.8 * ones(1, N);

% Simulação em malha aberta
v_ol = zeros(1, N);
v_ol(1) = 0;  % Velocidade inicial = 0
for k = 1:N-1
    v_ol(k+1) = simulate_vehicle(v_ol(k), u_ol(k), m, Ts);
end

% Figura 1: Resposta temporal em malha aberta
figure('Name', 'Malha Aberta - Tempo', 'Color', 'w');
plot(t, v_ol, 'b-', 'LineWidth', 1.5);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Velocidade (m/s)', 'FontSize', 12);
title('Resposta temporal do sistema em malha aberta', 'FontSize', 14);
grid on;
saveas(gcf, 'figs/1_open_loop_time.png');

% Figura 2: Resposta em frequência (FFT)
figure('Name', 'Malha Aberta - FFT', 'Color', 'w');
Y = fft(v_ol);
f_axis = (0:N-1)/(N*Ts);
plot(f_axis(1:floor(N/2)), abs(Y(1:floor(N/2))), 'b-', 'LineWidth', 1.5);
xlabel('Frequência (Hz)', 'FontSize', 12);
ylabel('Magnitude', 'FontSize', 12);
title('Resposta em frequência (FFT) do sistema em malha aberta', 'FontSize', 14);
grid on;
xlim([0, 5]); % Zoom nas baixas frequências
saveas(gcf, 'figs/2_open_loop_freq.png');

%% =========================================================================
%% 2. ANÁLISE DE ESTABILIDADE (Cálculo Manual - Sem Toolbox)
%% =========================================================================
fprintf('2. Gerando diagramas de estabilidade (Matemática Pura)...\n');

% Definição do vetor de frequências (logarítmico para Bode)
w = logspace(-3, 2, 1000); % de 10^-3 a 10^2 rad/s

% Função de Transferência G(s) = b / (s - a)
% Resposta em frequência G(jw)
s_jw = 1j * w;
G_jw = b_coef ./ (s_jw - a_coef);

% Magnitude (dB) e Fase (graus)
mag_db = 20 * log10(abs(G_jw));
phase_deg = rad2deg(angle(G_jw));

% Figura 3: Lugar das Raízes (Manual)
% Sistema de 1ª ordem: Polo em s = a_coef.
% LGR vai do polo para -infinito no eixo real.
figure('Name', 'Lugar das Raizes', 'Color', 'w');
plot(a_coef, 0, 'rx', 'MarkerSize', 12, 'LineWidth', 2); % Polo
hold on;
line([a_coef, -50], [0, 0], 'Color', 'b', 'LineWidth', 2); % Ramo do LGR
xlabel('Eixo Real (\sigma)', 'FontSize', 12);
ylabel('Eixo Imaginário (j\omega)', 'FontSize', 12);
title('Lugar das Raízes (Sistema de 1ª Ordem)', 'FontSize', 14);
legend('Polo de Malha Aberta', 'Lugar das Raízes');
grid on;
axis([-0.2 0.05 -0.1 0.1]); % Ajuste de zoom
saveas(gcf, 'figs/3_root_locus.png');

% Figura 4: Diagrama de Bode (Manual)
figure('Name', 'Diagrama de Bode', 'Color', 'w');
subplot(2,1,1);
semilogx(w, mag_db, 'b-', 'LineWidth', 1.5);
ylabel('Magnitude (dB)', 'FontSize', 11);
title('Diagrama de Bode', 'FontSize', 14);
grid on;
subplot(2,1,2);
semilogx(w, phase_deg, 'b-', 'LineWidth', 1.5);
xlabel('Frequência (rad/s)', 'FontSize', 12);
ylabel('Fase (graus)', 'FontSize', 11);
grid on;
saveas(gcf, 'figs/4_bode.png');

% Figura 5: Diagrama de Nyquist (Manual)
figure('Name', 'Diagrama de Nyquist', 'Color', 'w');
plot(real(G_jw), imag(G_jw), 'b-', 'LineWidth', 1.5);
hold on;
plot(real(G_jw), -imag(G_jw), 'b--', 'LineWidth', 1); % Espelho complexo conjugado
plot(-1, 0, 'r+', 'MarkerSize', 10, 'LineWidth', 2); % Ponto crítico
xlabel('Re(G(j\omega))', 'FontSize', 12);
ylabel('Im(G(j\omega))', 'FontSize', 12);
title('Diagrama de Nyquist', 'FontSize', 14);
grid on;
axis equal;
saveas(gcf, 'figs/5_nyquist.png');

%% =========================================================================
%% 3. CONTROLADOR PID
%% =========================================================================
fprintf('3. Simulando controlador PID...\n');

% Ganhos do PID
Kp = 0.8;
Ki = 0.05;
Kd = 0.1;

% Inicialização
v_pid = zeros(1, N);
u_pid = zeros(1, N);
e_int = 0;
e_prev = 0;

for k = 1:N-1
    % Erro
    error = v_ref - v_pid(k);
    
    % Termos PID (Integração e Derivação Discreta)
    e_int = e_int + error * Ts;
    e_deriv = (error - e_prev) / Ts;
    
    % Lei de Controle
    u_pid(k) = Kp * error + Ki * e_int + Kd * e_deriv;
    
    % Aplicação na Planta
    v_pid(k+1) = simulate_vehicle(v_pid(k), u_pid(k), m, Ts);
    
    e_prev = error;
end
u_pid(N) = u_pid(N-1); % Repetir última amostra para vetor ter tamanho N

% Figura 6: PID Combinado (Subplot) - REQUISITO DO EDITAL
figure('Name', 'PID - Resposta e Controle', 'Color', 'w', 'Position', [100, 100, 800, 600]);

% Subplot 1: Velocidade
subplot(2,1,1);
plot(t, v_pid, 'Color', [0.85, 0.33, 0.1], 'LineWidth', 1.5); hold on;
plot(t, v_ref*ones(size(t)), 'b--', 'LineWidth', 1.5);
ylabel('Velocidade (m/s)', 'FontSize', 11);
title('Resposta do Sistema com Controlador PID', 'FontSize', 12);
legend('Saída (PID)', 'Referência', 'Location', 'SouthEast');
grid on;
ylim([0, 25]);

% Subplot 2: Sinal de Controle
subplot(2,1,2);
plot(t, u_pid, 'b-', 'LineWidth', 1.5); hold on;
yline(2, 'k--', 'Limite Sup'); % Linha visual de limite físico
yline(-3, 'k--', 'Limite Inf');
xlabel('Tempo (s)', 'FontSize', 11);
ylabel('Aceleração (m/s^2)', 'FontSize', 11);
title('Sinal de Controle - PID (Sem Saturação Explícita)', 'FontSize', 12);
grid on;

saveas(gcf, 'figs/6_pid_combined.png');

%% =========================================================================
%% 4. PID COM VARIAÇÃO PARAMÉTRICA (+10% massa)
%% =========================================================================
fprintf('4. Teste de Robustez PID (+10%% massa)...\n');

m_var = m * 1.10;  % Massa +10%

% Inicialização Variável
v_pid_var = zeros(1, N);
u_pid_var = zeros(1, N);
e_int_var = 0;
e_prev_var = 0;

for k = 1:N-1
    error = v_ref - v_pid_var(k);
    e_int_var = e_int_var + error * Ts;
    e_deriv_var = (error - e_prev_var) / Ts;
    
    u_pid_var(k) = Kp * error + Ki * e_int_var + Kd * e_deriv_var;
    v_pid_var(k+1) = simulate_vehicle(v_pid_var(k), u_pid_var(k), m_var, Ts); % Usa m_var
    
    e_prev_var = error;
end

% Figura 7: Comparação PID nominal vs +10% massa
figure('Name', 'Robustez PID', 'Color', 'w');
plot(t, v_pid, 'b-', 'LineWidth', 1.5); hold on;
plot(t, v_pid_var, 'r--', 'LineWidth', 1.5);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Velocidade (m/s)', 'FontSize', 12);
title('Robustez do PID: Variação de Parâmetros', 'FontSize', 14);
legend('PID (Nominal)', 'PID (+10% Massa)', 'Location', 'SouthEast');
grid on;
ylim([0, 25]);
saveas(gcf, 'figs/7_pid_robustez.png');

%% =========================================================================
%% 5. CONTROLADOR MPC (Otimização Quadrática Simplificada)
%% =========================================================================
fprintf('5. Simulando controlador MPC...\n');

% Parâmetros do MPC
u_min = -3;        % Limite inferior de aceleração [m/s²]
u_max = 2;         % Limite superior de aceleração [m/s²]
Q = 10;            % Peso do erro
R = 0.1;           % Peso do esforço de controle

% Ganho estático aproximado para lei de controle MPC univariável
% u = K * erro (saturado)
K_mpc = Q / (Q + R); 

% Inicialização
v_mpc = zeros(1, N);
u_mpc = zeros(1, N);

for k = 1:N-1
    e_mpc = v_ref - v_mpc(k);
    
    % Cálculo da ação ótima sem restrição
    u_opt = K_mpc * e_mpc;
    
    % Aplicação das restrições (Saturação)
    u_mpc(k) = max(u_min, min(u_max, u_opt));
    
    % Aplicação na Planta
    v_mpc(k+1) = simulate_vehicle(v_mpc(k), u_mpc(k), m, Ts);
end
u_mpc(N) = u_mpc(N-1);

% Figura 8: MPC Combinado (Subplot) - REQUISITO DO EDITAL
figure('Name', 'MPC - Resposta e Controle', 'Color', 'w', 'Position', [100, 100, 800, 600]);

% Subplot 1: Velocidade
subplot(2,1,1);
plot(t, v_mpc, 'Color', [0.85, 0.33, 0.1], 'LineWidth', 1.5); hold on;
plot(t, v_ref*ones(size(t)), 'b--', 'LineWidth', 1.5);
ylabel('Velocidade (m/s)', 'FontSize', 11);
title('Resposta do Sistema com Controlador MPC', 'FontSize', 12);
legend('Saída (MPC)', 'Referência', 'Location', 'SouthEast');
grid on;
ylim([0, 25]);

% Subplot 2: Sinal de Controle
subplot(2,1,2);
plot(t, u_mpc, 'r-', 'LineWidth', 1.5); hold on;
plot(t, u_max*ones(size(t)), 'k--', 'LineWidth', 1);
plot(t, u_min*ones(size(t)), 'k--', 'LineWidth', 1);
xlabel('Tempo (s)', 'FontSize', 11);
ylabel('Aceleração (m/s^2)', 'FontSize', 11);
title('Sinal de Controle - MPC (Com Limites Ativos)', 'FontSize', 12);
legend('u(t) MPC', 'Restrições [-3, 2]', 'Location', 'best');
grid on;
ylim([-4, 3]);

saveas(gcf, 'figs/8_mpc_combined.png');

%% =========================================================================
%% 6. NOVO: MPC COM VARIAÇÃO PARAMÉTRICA (+10% massa) - REQUISITO DO EDITAL
%% =========================================================================
fprintf('6. Teste de Robustez MPC (+10%% massa)...\n');

% Simulação MPC com massa variada
v_mpc_var = zeros(1, N);
u_mpc_var = zeros(1, N);

for k = 1:N-1
    e_mpc_var = v_ref - v_mpc_var(k);
    
    u_opt_var = K_mpc * e_mpc_var;
    u_mpc_var(k) = max(u_min, min(u_max, u_opt_var));
    
    % Simula com m_var
    v_mpc_var(k+1) = simulate_vehicle(v_mpc_var(k), u_mpc_var(k), m_var, Ts);
end

% Figura 9: Comparação MPC nominal vs +10% massa
figure('Name', 'Robustez MPC', 'Color', 'w');
plot(t, v_mpc, 'b-', 'LineWidth', 1.5); hold on;
plot(t, v_mpc_var, 'r--', 'LineWidth', 1.5);
plot(t, v_ref*ones(size(t)), 'k:', 'LineWidth', 1);
xlabel('Tempo (s)', 'FontSize', 12);
ylabel('Velocidade (m/s)', 'FontSize', 12);
title('Robustez do MPC: Variação de Parâmetros', 'FontSize', 14);
legend('MPC (Nominal)', 'MPC (+10% Massa)', 'Referência', 'Location', 'SouthEast');
grid on;
ylim([0, 25]);
saveas(gcf, 'figs/9_mpc_robustez.png');

%% =========================================================================
%% 7. COMPARAÇÃO GERAL E MÉTRICAS
%% =========================================================================
fprintf('7. Comparação Final e Métricas...\n');

% Figura 10: Comparação PID vs MPC
figure('Name', 'Comparativo Final', 'Color', 'w', 'Position', [100, 100, 800, 700]);
subplot(2,1,1);
plot(t, v_ref*ones(size(t)), 'k--', 'LineWidth', 1.5); hold on;
plot(t, v_pid, 'b-', 'LineWidth', 1.5);
plot(t, v_mpc, 'r-', 'LineWidth', 1.5);
ylabel('Velocidade (m/s)');
title('Comparação de Resposta');
legend('Ref', 'PID', 'MPC', 'Location', 'SouthEast');
grid on; ylim([0 25]);

subplot(2,1,2);
plot(t, u_pid, 'b-', 'LineWidth', 1.5); hold on;
plot(t, u_mpc, 'r-', 'LineWidth', 1.5);
yline(2, 'k--'); yline(-3, 'k--');
xlabel('Tempo (s)'); ylabel('Aceleração (m/s^2)');
title('Comparação de Esforço de Controle');
legend('PID', 'MPC', 'Limites');
grid on; ylim([-5 20]);
saveas(gcf, 'figs/10_comparacao_final.png');

% Cálculo de Métricas
e_pid_vec = v_ref - v_pid;
e_mpc_vec = v_ref - v_mpc;

RMSE_pid = sqrt(mean(e_pid_vec.^2));
RMSE_mpc = sqrt(mean(e_mpc_vec.^2));
IAE_pid = sum(abs(e_pid_vec)) * Ts;
IAE_mpc = sum(abs(e_mpc_vec)) * Ts;
ISE_pid = sum(e_pid_vec.^2) * Ts;
ISE_mpc = sum(e_mpc_vec.^2) * Ts;
ITAE_pid = sum(t .* abs(e_pid_vec)) * Ts;
ITAE_mpc = sum(t .* abs(e_mpc_vec)) * Ts;

fprintf('\n========== MÉTRICAS DE DESEMPENHO ==========\n');
fprintf('Métrica    |    PID    |    MPC    | Razão (MPC/PID)\n');
fprintf('----------------------------------------------\n');
fprintf('RMSE       | %8.4f  | %8.4f  | %.2fx\n', RMSE_pid, RMSE_mpc, RMSE_mpc/RMSE_pid);
fprintf('IAE        | %8.4f  | %8.4f  | %.2fx\n', IAE_pid, IAE_mpc, IAE_mpc/IAE_pid);
fprintf('ISE        | %8.4f  | %8.4f  | %.2fx\n', ISE_pid, ISE_mpc, ISE_mpc/ISE_pid);
fprintf('ITAE       | %8.4f  | %8.4f  | %.2fx\n', ITAE_pid, ITAE_mpc, ITAE_mpc/ITAE_pid);
fprintf('==============================================\n');

%% =========================================================================
%% 8. ESPECTROGRAMAS (Cálculo Manual - Sem Toolbox)
%% =========================================================================
fprintf('8. Gerando Espectrogramas (Método Manual STFT)...\n');

% Parâmetros da STFT
n_win = 64;             % Tamanho da janela
n_overlap = 56;         % Sobreposição
n_fft = 128;            % Pontos da FFT
fs_sample = 1/Ts;       % Frequência de amostragem (10 Hz)
hop_size = n_win - n_overlap;

% Função de Janela (Hamming feita na mão) - Vetor COLUNA
win_func = 0.54 - 0.46 * cos(2*pi*(0:n_win-1)' / (n_win-1));

% --- Função Auxiliar para Calcular STFT (CORRIGIDA) ---
% A correção está no (fft(...).') -> Transpõe para coluna antes de juntar
calc_stft = @(signal) ...
    cell2mat(arrayfun(@(k) abs(fft(signal(k:k+n_win-1) .* win_func', n_fft).').^2, ...
    1:hop_size:(length(signal)-n_win), 'UniformOutput', false));

% --- 1. Espectrograma PID ---
S_pid = calc_stft(u_pid); % Agora retorna matriz (Frequencia x Tempo)
S_pid = S_pid(1:n_fft/2+1, :); % Pegar apenas frequências positivas (0 a fs/2)
S_pid_db = 10*log10(S_pid + eps); % Converter para dB

% Eixos de tempo e frequência para plotagem
t_spec = (0:size(S_pid, 2)-1) * hop_size * Ts;
f_spec = (0:n_fft/2) * fs_sample / n_fft;

figure('Name', 'Espectrograma PID', 'Color', 'w');
imagesc(t_spec, f_spec, S_pid_db);
axis xy; % Coloca a frequência 0 na parte inferior
colormap('jet'); 
colorbar;
xlabel('Tempo (s)'); 
ylabel('Frequência (Hz)');
title('Espectrograma do Sinal de Controle (PID)');
try a = gca; a.Toolbar.Visible = 'off'; catch; end 
saveas(gcf, 'figs/11_spectrogram_pid.png');

% --- 2. Espectrograma MPC ---
S_mpc = calc_stft(u_mpc);
S_mpc = S_mpc(1:n_fft/2+1, :);
S_mpc_db = 10*log10(S_mpc + eps);

figure('Name', 'Espectrograma MPC', 'Color', 'w');
imagesc(t_spec, f_spec, S_mpc_db);
axis xy;
colormap('jet'); 
colorbar;
xlabel('Tempo (s)'); 
ylabel('Frequência (Hz)');
title('Espectrograma do Sinal de Controle (MPC)');
try a = gca; a.Toolbar.Visible = 'off'; catch; end
saveas(gcf, 'figs/12_spectrogram_mpc.png');

fprintf('\nSIMULAÇÃO FINALIZADA. Verifique a pasta "figs/".\n');
