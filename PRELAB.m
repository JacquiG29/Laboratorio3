%% Tarea 1 Control 2 Problema 7
% Jacqui y Gaby
% 16009
% 21/02/2020
%% Parte 1
% Sugerencia
R1 = 1*10^3;
R2 = 10*10^3;
R3 = R2;
C1 = 1*10^(-6);
C2 = 0.1*10^(-6);
C3 = 10*10^(-6);

% Numerador
b0=1/(R1*R2*R3*C1*C2*C3);
b1 = 0;
b2 = 0;

% Denominador
a0 = b0;
a1=(1/(R2*R3*C2*C3))+((R3+R2+R1)/(R1*R2*R3*C1*C3));
a2=((C3*R2*R3+C3*R1*R3+C1*R1*R2+C1*R1*R3)/(R2*R3*R1*C3*C1));
G0=tf(b0,[1,a2,a1,a0]);

%% Parte 2
% Realización canónica controlable
A_controlable = [0    1   0;...
                 0    0   1;...
                 -a0 -a1 -a2];
B_controlable = [0;0;1];
C_controlable = [b0,b1,b2];

% Realización canónica observable
A_observable = [-a2    1   0;...
                -a1    0   1;...
                -a0    0   0];
B_observable = [b2; b1; b0];
C_observable = [1,0,0];

% A partir de las ecuaciones diferenciales
A_ed =[-1*(R1+R2)/(R1*R2*C1)   1/(R2*C1)      -1/(R2*C1)     ;...
         0                 0              -1/(R3*C2)     ;...
       -1/(R2*C3)          1/(R2*C3)   -(R2+R3)/(R2*R3*C3)];
B_ed = [1/(R1*C1); 0; 0];
C_ed = [0 1 0];

% A partir de la función de transferencia G0
[A_g0,B_g0,C_g0,D_g0] = tf2ss(b0,[1,a2,a1,a0]);

% A partir del Linear Analysis Tool
load('linsys.mat');
A_LAT = linsys1.A;
B_LAT = linsys1.B;
C_LAT = linsys1.C;

%% Función de transferencia a partir de las matrices
% G(s) = C*inv(sI-A)*B+D
% Primero definimos la s como transfer function para no utilizar el paquete
% simbólico de Matlab.
s = tf('s');

% TF Función de transferencia G0
G_g0 = C_g0*inv(s*eye(3)-A_g0)*B_g0;

% TF Canónica Controlable
G_controlable = C_controlable*inv(s*eye(3)-A_controlable)*B_controlable;

% TF Canónica observable
G_observable = C_observable*inv(s*eye(3)-A_observable)*B_observable;

% TF Ecuaciones diferenciales
G_ed = C_ed*inv(s*eye(3)-A_ed)*B_ed;

% TF Linear Analysis Tool
G_LAT = C_LAT*inv(s*eye(3)-A_LAT)*B_LAT;

%% Comprobación de las funciones de transferencia

G0
G_controlable
G_observable
G_ed
G_LAT


