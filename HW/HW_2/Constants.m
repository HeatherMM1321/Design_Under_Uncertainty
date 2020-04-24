%% Motor Design Variables x
D   = ?   % rotor diameter (cm)
L   = ?;  % rotor axial length (cm)
Lwa = 14.12;                % armature wire length (m)
Lwf = 309.45;               % field wire length (m)
ds  = 0.612;                % slot depth (cm)

%% Coupling Variables b, shared with control problem
n   = 122;                   % rotational speed (rev/s)
V   = 40;                   % design voltage (V)
Pmin= 3.94;                 % minimum required power (kW)
Tmin= 5.12e-3;         	% minimum required torque (kNm)

%% Parameter Vector a (constants)
dcu = ?;               % copper density at 25C (kg/m^3)
dfe = ?;               % iron density at 25C (kg/m^3)
fi  = 0.7;                  % pole arc to pole pitch ratio
p   = 2;                    % number of poles
S   = 27;                   % number of slots (teeth on rotor)
rho = 1.8e-8;               % resistivity (ohm-m) of copper at 25C


example of the b vector:

b = [n V Pmin Tmin];