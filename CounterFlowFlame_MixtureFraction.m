function [solutionArray, lambdaOut] = CounterFlowFlame_MixtureFraction(myconfig)
options = bvpset('stats', 'off', 'NMax', 10000, 'AbsTol', 1e-8);

%% Solver configuration
lambda = -100; %Initial estimation for dpdz
TF0 = myconfig.TF0;
TO0 = myconfig.TO0;
Q = myconfig.Q;
cp = myconfig.cp;
zst = myconfig.zst;
L = myconfig.L;
YF0 = myconfig.fuelInletConcentration(1);
YO0 = myconfig.oxidizerInletConcentration(2);
Coef_Stoic = myconfig.Coef_Stoic;
MM = myconfig.MM;
Tref = myconfig.Tref;
a = myconfig.a;
Prandtl = 1;

%% Solve system
sol = bvpinit(linspace(0, L, myconfig.initialCellNumber), @mat4init, lambda);
sol = bvp4c(@mat4ode, @mat4bc, sol, options);
Sxint = deval(sol, linspace(0, L, 200));
lambdaOut = sol.parameters;
Sxint = deval(sol, sol.x);

%% Recover primitive variables from Z
ZArray = Sxint(4, :);
for i = 1:length(ZArray)
    z_T(i) = getTemperatureFromZ(ZArray(i));
    z_Y1(i) = getY1FromZ(ZArray(i));
    z_Y2(i) = getY2FromZ(ZArray(i));
    z_Y3(i) = getY3FromZ(ZArray(i));
    z_Y4(i) = getY4FromZ(ZArray(i));
    z_Y5(i) = getY5FromZ(ZArray(i));
end
solutionArray = [sol.x; Sxint(1, :); Sxint(2, :); Sxint(3, :); z_T; z_Y1; z_Y2; z_Y3; z_Y4; z_Y5; Sxint(4, :)];

%% =======================
% Note: the system to be solved is of first order on v, second order on U,
% and second order on the scalars (T, Y1,Y2...)
% Second order equations are brought to a first order by introducing a
% transformation
function dydx = mat4ode(~, y, lambda) % equation being solved
    v = y(1);
    U1 = y(2);
    U2 = y(3);
    Z1 = y(4); % Z
    Z2 = y(5); % dZ/dy
    if (Z1 > 1) %Limit values
        Z1 = 1;
    end
    if (Z1 < 0)
        Z1 = 0;%Limit values
    end

    rho_ = GetRhoFromZ(Z1);
    mu_ = GetMuFromZ(Z1) / Prandtl;
    % Calculation of density and viscosity derivatives is problematic,
    % because the derivative at zSt is not defined.
    %Is simply ignored (seems to still give an adequate starting solution )
    rho_p = 0.0;
    mu_p = 0.0;

    dydx = [(-rho_ * U1 - rho_p * v) / rho_; ... %Conti
        U2; ... %Mom
        (lambda + rho_ * v * U2 + rho_ * U1^2 - mu_p * U2) * 1.0 / mu_; ... %Mom
        Z2; ... .%MassFraction
        (rho_ * v * Z2 - mu_p * Z2) / mu_; ... %MassFraction
        ];
end


function res = mat4bc(ya, yb, lambda) % boundary conditions
    va = ya(1);
    U1a = ya(2);
    Za = ya(4);

    vb = yb(1);
    U1b = yb(2);
    Zb = yb(4);

    res = [va - myconfig.vLeft; ... % v(0) = v0
        vb - myconfig.vRight; ... % v(L) = vl
        U1a - 0; ... % U(0) = 0
        U1b - 0; ... % U(L) = 0
        Za - 1.0; ... % Z(0) = 1.0
        Zb - 0.0; ... %Z(L) = TL
        ];
end

function yinit = mat4init(x) % initial guess function.
    yinit = [0.0; ... %v
        0.0; ... %U1
        0.0; ... %U2
        0.5; ... % Z1
        0.0; ... % Z2
        ];
end
% Definition of helper functions for density and transport parameters
function viscosity = mu(T)
    viscosity = myconfig.viscosity0 * (T / Tref)^(a);
end
function density = rho(T, Yk)
    mult = 0.0;
    for c = 1:(length(Yk))
        mult = mult + Yk(c) / MM(c);
    end
    density = myconfig.p0 / (myconfig.R * T * mult); % kg / m^3
end
%==============================================
%Transformations from Z

function Temperature_fromZ = getTemperatureFromZ(Z)
    if (Z >= zst)
        Temperature_fromZ = Z * TF0 + (1 - Z) * TO0 + Q * YF0 / cp * zst * (1 - Z) / (1 - zst);
    else
        Temperature_fromZ = Z * TF0 + (1 - Z) * TO0 + Q * YF0 / cp * Z;
    end
end

function Y1_fromZ = getY1FromZ(Z)
    if (Z >= zst)
        Y1_fromZ = YF0 * (Z - zst) / (1 - zst);
    else
        Y1_fromZ = 0.0;
    end
end

function Y2_fromZ = getY2FromZ(Z)
    if (Z >= zst)
        Y2_fromZ = 0;
    else
        Y2_fromZ = YO0 * (1 - Z / zst);
    end
end

function Y3_fromZ = getY3FromZ(Z)
    nu_CO2 = Coef_Stoic(3);
    nu_O2 = Coef_Stoic(2);
    nu_CH4 = Coef_Stoic(1);
    MW_CO2 = MM(3);
    MW_O2 = MM(2);
    MW_CH4 = MM(1);
    if (Z >= zst)
        Y3_fromZ = -YO0 * (nu_CO2 * MW_CO2) / (nu_O2 * MW_O2) * (1 - Z);
    else
        Y3_fromZ = -YF0 * (nu_CO2 * MW_CO2) / (nu_CH4 * MW_CH4) * Z;
    end
end

function Y4_fromZ = getY4FromZ(Z)
    nu_H2O = Coef_Stoic(4);
    nu_O2 = Coef_Stoic(2);
    nu_CH4 = Coef_Stoic(1);
    MW_H2O = MM(4);
    MW_O2 = MM(2);
    MW_CH4 = MM(1);
    if (Z >= zst)
        Y4_fromZ = -YO0 * (nu_H2O * MW_H2O) / (nu_O2 * MW_O2) * (1 - Z);
    else
        Y4_fromZ = -YF0 * (nu_H2O * MW_H2O) / (nu_CH4 * MW_CH4) * Z;
    end
end

function Y5_fromZ = getY5FromZ(Z)
    if (Z >= zst)
        YNOxi0 = 1.0 - YO0;
        YNFuel0 = 1.0 - YF0;
        Y5_fromZ = YNOxi0 * (1 - Z) + YNFuel0 * Z;
    else
        YNOxi0 = 1.0 - YO0;
        YNFuel0 = 1.0 - YF0;
        Y5_fromZ = YNOxi0 * (1 - Z) + YNFuel0 * Z;
    end
end

function MuFromZ = GetMuFromZ(z)
    MuFromZ = mu(getTemperatureFromZ(z));
end

function rhoFromZ = GetRhoFromZ(Z)
    rhoFromZ = rho(getTemperatureFromZ(Z), [getY1FromZ(Z), getY2FromZ(Z), getY3FromZ(Z), getY4FromZ(Z), getY5FromZ(Z)]);
end
end
