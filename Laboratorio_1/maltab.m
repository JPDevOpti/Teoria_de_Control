clear all; clc; close all;
%%
F1 = tf([1 4], [1 7 13]); %Función de transferencia 1 dada
F2 = tf([10 5], [1 2 5 3]);%Función de transferencia 2 dada
[p,z] = pzmap(F2); %Obtener polos y ceros
disp(p) % Imprimir los polos complejos conjugados que generan la respuesta subamortiguada del sistema :o

%Realimentaciones negativas unitarias para evaluar en lazo cerrado
y1=feedback(F1,1);
y2=feedback(F2,1);

t= 0:0.1:10; %Vector de tiempo para graficar

%Valores de las constantes para el punto de la bicicleta
m1=8; m2=5; 
k1=7; k2=7; b=18;

%Extraer datos de Simulink
out=sim('Informe1.slx'); 

%% Entrada escalón G1

%LAZO ABIERTO

%Obtención datos Simulink
StepF1 = out.F1_Step.signals.values();
tiempo = out.F1_Step.time();

figure(1)
sgtitle('Respuesta ante el escalón (G1)')

subplot(2,1,1)

[escalonF1LA, time] = step(F1,4); %MatLab 
plot(time, escalonF1LA, LineWidth=1)
hold on;

plot(tiempo, StepF1,'--y', LineWidth=1); %Simulink

title('Lazo abierto');
grid on ;
legend('MatLab', 'Simulink');
xlabel('Tiempo (s)')
ylabel('Amplitud')
xlim([0,4])
ylim([0,0.35])


% LAZO CERRADO

StepF1_LC = out.F1_Step_LC.signals.values();
tiempo = out.F1_Step_LC.time();

subplot(2,1,2)

[escalonF1LC, time] = step(y1,4); %MatLab
plot(time, escalonF1LC,LineWidth=1);
hold on;

plot(tiempo, StepF1_LC,'--y',LineWidth=1); %Simulink

title('Lazo cerrado');
grid on ;
legend('MatLab', 'Simulink');
xlabel('Tiempo (s)')
ylabel('Amplitud')
xlim([0,4])
ylim([0,0.35])


%% Entrada impulso G1

ImpulseF1 = out.F1_Impulse.signals.values(:,:);
tiempo = out.F1_Impulse.time();

%LAZO ABIERTO
figure(2)

sgtitle('Respuesta ante el impulso (G1)')

subplot(2,1,1)

[impulseF1LA, time] = impulse(F1,4); %MatLab
plot(time,impulseF1LA,LineWidth=1)
hold on;
plot(tiempo, ImpulseF1,'--y',LineWidth=1); %Simulink

title('Lazo abierto');
grid on ;
legend('MatLab', 'Simulink');
xlabel('Tiempo (s)')
ylabel('Amplitud')
xlim([0,4])


%LAZO CERRADO

ImpulseF1_LC = out.F1_Impulse_LC.signals.values(:,:);
tiempo = out.F1_Impulse_LC.time();

subplot(2,1,2)

[impulseF1LC, time] = impulse(y1,4); %MatLab
plot(time,impulseF1LC,LineWidth=1)
hold on;
plot(tiempo, ImpulseF1_LC,'--y',LineWidth=1); %Simulink

title('Lazo cerrado');
grid on ;
legend('MatLab', 'Simulink');
xlabel('Tiempo (s)')
ylabel('Amplitud')
xlim([0,4])

%% Entrada rampa G1

%LAZO ABIERTO

RampF1 = out.F1_Ramp.signals.values();
tiempo = out.F1_Ramp.time();

figure(3)
sgtitle('Respuesta ante la función rampa (G1)')

subplot(2,1,1)

[rampaF1LA, time]= lsim(F1, t, t); % MatLab
plot(time,rampaF1LA,LineWidth=1)
hold on;
plot(tiempo, RampF1,'--y',LineWidth=1); %Simulink

title('Lazo abierto');
grid on ;
legend('MatLab', 'Simulink');
xlabel('Tiempo (s)')
ylabel('Amplitud')
xlim([0,6])
ylim([0,2])

% LAZO CERRADO

RampF1_LC = out.F1_Ramp_LC.signals.values();
tiempo = out.F1_Ramp_LC.time();

subplot(2,1,2)

[rampaF1LC, time]=lsim(y1, t, t); % MatLab
plot(time,rampaF1LC,LineWidth=1)
hold on;
plot(tiempo, RampF1_LC,'--y',LineWidth=1); % Simulink

title('Lazo cerrado');
grid on ;
legend('MatLab', 'Simulink');
xlabel('Tiempo (s)')
ylabel('Amplitud')
xlim([0,6])
ylim([0,2])


%% Respuesta al Escalón (G2)

% LAZO ABIERTO

StepF2 = out.F2_Step.signals.values();
tiempo = out.F2_Step.time();

figure(4)
sgtitle('Respuesta ante el escalón (G2)')

subplot(2,1,1)

[escalonF2LA, time]= step(F2,10); %Matlab
plot(time,escalonF2LA, LineWidth=1)
hold on;

plot(tiempo, StepF2,'--y', LineWidth=1); %Simulink

title('Lazo abierto');
grid on ;
legend('MatLab', 'Simulink');
xlabel('Tiempo (s)')
ylabel('Amplitud')


% LAZO CERRADO

StepF2_LC = out.F2_Step_LC.signals.values();
tiempo = out.F2_Step_LC.time();

subplot(2,1,2)

[escalonF2LC, time]=step(y2,10); %MatLab
plot(time,escalonF2LC, LineWidth=1)
hold on;

plot(tiempo, StepF2_LC,'--y', LineWidth=1); %Simulink

title('Lazo cerrado');
grid on ;
legend('MatLab', 'Simulink');
xlabel('Tiempo (s)')
ylabel('Amplitud')


%% Respuesta al impulso (G2)

% LAZO ABIERTO

ImpulseF2 = out.F2_Impulse.signals.values(:,:);
tiempo = out.F2_Impulse.time();

figure(5)
sgtitle('Respuesta ante el impulso (G2)')

subplot(2,1,1)

[impulseF2LA, time]=impulse(F2,10); %MatLab
plot(time,impulseF2LA,LineWidth=1)
hold on;
plot(tiempo, ImpulseF2,'--y',LineWidth=1); %Simulink

title('Lazo abierto');
grid on ;
legend('MatLab', 'Simulink');
xlabel('Tiempo (s)')
ylabel('Amplitud')

% LAZO CERRADO

ImpulseF2_LC = out.F2_Impulse_LC.signals.values(:,:);
tiempo = out.F2_Impulse_LC.time();

subplot(2,1,2)

[impulseF2LC, time]=impulse(y2,10); %MatLab
plot(time,impulseF2LC,LineWidth=1)
hold on;
plot(tiempo, ImpulseF2_LC,'--y',LineWidth=1); %Simulink

title('Lazo cerrado');
grid on ;
legend('MatLab', 'Simulink');
xlabel('Tiempo (s)')
ylabel('Amplitud')

%% Respuesta ante fn Rampa (G2) 

% LAZO ABIERTO
RampF2 = out.F2_Ramp.signals.values();
tiempo = out.F2_Ramp.time();

figure(6)
sgtitle('Respuesta ante la función rampa (G2)')

subplot(2,1,1)

[rampaF2LA, time] = lsim(F2, t, t); %MatLab
plot(time,rampaF2LA,LineWidth=1)
hold on;
plot(tiempo, RampF2,'--y',LineWidth=1); %Simulink

title('Lazo abierto');
grid on ;
legend('MatLab', 'Simulink');
xlabel('Tiempo (s)')
ylabel('Amplitud')

% LAZO CERRADO 

RampF2_LC = out.F2_Ramp_LC.signals.values();
tiempo = out.F2_Ramp_LC.time();

subplot(2,1,2)

[rampaF2LC, time] = lsim(y2, t, t); %MatLab
plot(time,rampaF2LC,LineWidth=1)
hold on;
plot(tiempo, RampF2_LC,'--y',LineWidth=1); %Simulink

title('Lazo cerrado');
grid on ;
legend('MatLab', 'Simulink');
xlabel('Tiempo (s)')
ylabel('Amplitud')

%% SUSPENSIÓN DE UNA BICICLETA

FTXU = tf([k1*m2 k1*b k2*k1], [m2*m1 (m2*b)+(b*m1) (m2*k1)+(m2*k2)+(k2*m1) b*k1 k2*k1]); %FT = X/U
FTYU = tf([b*k1 k2*k1], [m2*m1 (m2*b)+(b*m1) (m2*k1)+(m2*k2)+(k2*m1) b*k1 k2*k1]); %FT = Y/U
tiempo = 0:0.1:80; %Vector de tiempo para graficar

[p,z] = pzmap(FTXU); %Obtener polos y ceros
disp(p) % Imprimir los polos para ver los complejos conjugados que generan 
% la respuesta subamortiguada del sistema :o

%Realimentaciones para realizar el análisis en lazo cerrado con
%realimentación negativa unitaria
XU_feedback = feedback(FTXU,1); 
YU_feedback = feedback(FTYU,1);

%% Respuesta al escalón

%LAZO ABIERTO

% Para G = X/U
[escalonXU_LA,tiempobici] = step(FTXU,tiempo);

% Para G = Y/U
[escalonYU_LA,tiempobici] = step(FTYU,tiempo);


% LAZO CERRADO

% Para G = X/U
[escalonXU_LC,tiempobici] = step(XU_feedback,tiempo);

% Para G = Y/U
[escalonYU_LC,tiempobici] = step(YU_feedback,tiempo);

%% Respuesta al impulso

% LAZO ABIERTO

% Para G = X/U
[impulsoXU_LA,tiempobici] = impulse(FTXU,tiempo);

% Para G = Y/U
[impulsoYU_LA,tiempobici] = impulse(FTYU,tiempo);


% LAZO CERRADO

% Para G = X/U
[impulsoXU_LC,tiempobici] = impulse(XU_feedback,tiempo);

% Para G = Y/U
[impulsoYU_LC,tiempobici] = impulse(YU_feedback,tiempo);

%% Respuesta a la fn rampa

% LAZO ABIERTO

rampa_bici = (tiempo >=0).*tiempo;

% Para G = X/U
rampaXU_LA = lsim(FTXU, rampa_bici, tiempo);

% Para G = Y/U
rampaYU_LA = lsim(FTYU, rampa_bici, tiempo);


% LAZO CERRADO

% Para G = X/U
rampaXU_LC = lsim(XU_feedback, rampa_bici, tiempo);

% Para G = Y/U
rampaYU_LC = lsim(YU_feedback, rampa_bici, tiempo);

%% Graficación
figure(7);
% Para la primera función de transferencia G=X/U

sgtitle("Respuestas TF1=X/U");

subplot(3,1,1)

plot(tiempobici,escalonXU_LA,LineWidth=1)
hold on
plot(tiempobici,escalonXU_LC,LineWidth=1)
title("Respuesta al escalón")
legend("Lazo abierto", "Lazo cerrado")
ylabel("Amplitud")
xlabel("Tiempo (s)")

subplot(3,1,2)

plot(tiempobici,impulsoXU_LA,LineWidth=1)
hold on
plot(tiempobici,impulsoXU_LC,LineWidth=1)
title("Respuesta al impulso")
legend("Lazo abierto", "Lazo cerrado")
ylabel("Amplitud")
xlabel("Tiempo (s)")

subplot(3,1,3)

plot(tiempobici,rampaXU_LA,LineWidth=1)
hold on
plot(tiempobici,rampaXU_LC,LineWidth=1)
title("Respuesta a la rampa")
legend("Lazo abierto", "Lazo cerrado")
ylabel("Amplitud")
xlabel("Tiempo (s)")

%%
% Para la segunda función de transferencia G=Y/U
figure(8);

sgtitle("Para TF2=Y/U");

subplot(3,1,1)
plot(tiempo,escalonYU_LA,LineWidth=1)
hold on
plot(tiempo,escalonYU_LC,LineWidth=1)
title("Respuesta al escalón")
legend("Lazo abierto", "Lazo cerrado")
ylabel("Amplitud")
xlabel("Tiempo (s)")

subplot(3,1,2)
plot(tiempo,impulsoYU_LA,LineWidth=1)
hold on
plot(tiempo,impulsoYU_LC,LineWidth=1)
title("Respuesta al impulso")
legend("Lazo abierto", "Lazo cerrado")
ylabel("Amplitud")
xlabel("Tiempo (s)")

subplot(3,1,3)
plot(tiempo,rampaYU_LA,LineWidth=1)
hold on
plot(tiempo,rampaYU_LC,LineWidth=1)
title("Respuesta a la rampa")
legend("Lazo abierto", "Lazo cerrado")
ylabel("Amplitud")
xlabel("Tiempo (s)")





