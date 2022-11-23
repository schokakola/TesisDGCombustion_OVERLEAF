clc; clear all; close all;
%% Velocity multiplier which defines the strain
velMult = 3; 
%% Problem geometry description
C.L = 0.02; % Domain length, m
C.cp = 1.3; % Heat capacity;
C.initialCellNumber = 100;
C.chemActive = true;
C.variableKineticParameters = true;
%% Boundary conditions
C.vLeft = 0.0243 * velMult; % Velocity left (fuel)
C.vRight = -0.0243 * velMult * 3; % Velocity right (oxidizer)
C.TF0 = 300; % Temperature left (fuel)
C.TO0 = 300; % Temperature right (oxidizer)
C.CompositionFuel = [0.2, 0.0, 0.0, 0.0, 0.8];
C.CompositionOx = [0.0, 0.23, 0.0, 0.0, 0.77];
%% Physical parameters
C.p0 = 101325; % Pressure, Pa
C.viscosity0 = 1.716e-5; % kg/( m s) ==> viscosity at T = 273.15 for air
C.Tref = 273; % Reference temperature from powerlaw
C.a = 2 / 3; % PowerLaw exponent
C.Coef_Stoic = [-1, -2, 1, 2, 0]; % Chemical reaction 1CH4 + 2O2 -> 1CO2 + 2H2O
C.MM = [16, 32, 44, 18, 28];
C.s = (C.Coef_Stoic (2) * C.MM(2)) / (C.Coef_Stoic(1) * C.MM(1));
C.phi = C.s * C.CompositionFuel(1) / C.CompositionOx(2);
C.zst = 1.0 / (1.0 + C.phi);
C.Q = 50100;
C.R = 8.314 * 1000; % gas constant, kg m^2 / (s^2 K mol)
%% Calculate flame sheet (infinite reaction rate)
x = zeros(1);
FlameSheetSolution = zeros(1);
lambdaOut = 100;
[sol, lambdaOut] = CounterFlowFlame_MixtureFraction(C);
x = sol(1, :); % x coordinates
dummyZero = zeros(1, length(x)); 
FlameSheetSolution = [sol(2, :); sol(3, :); sol(4, :); sol(5, :); ...
    dummyZero; sol(6, :); dummyZero; sol(7, :); dummyZero; sol(8, :); ...
    dummyZero; sol(9, :); dummyZero; sol(10, :); dummyZero];
%% Finite reaction rate calculation
[xint, Sxint, calculatedStrain, calculatedMaxTemperature, lambdaOut] = CounterFlowFlame_FullChemistry(C, FlameSheetSolution, x, lambdaOut, true);

% Export results for comparison with two dimensional formulation