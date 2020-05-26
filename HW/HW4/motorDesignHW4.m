function varargout = motorDesign(D, L, dcu, dfe, nargout)
% DC Motor design function
%nargout=1;

% Copyright 2010 The MathWorks, Inc.
%% Design Vector X
% Convert to consistent units
D   = D*1e-2;            % rotor diameter (m)
L   = L*1e-2;            % rotor axial length (m)
Lwa = 14.12;                 % armature wire length (m)
Lwf = 309.45;                 % field wire length (m)
ds  = 0.612*1e-2;            % slot depth (m)

%% Parameter Vector a (constants)
dcu = dcu;                 % copper density density at 25C (kg/m^3)
dfe = dfe;                 % iron density density at 25C (kg/m^3)
fi  = 0.7;                  % pole arc to pole pitch ratio
p   = 2;                    % number of poles
S   = 27;                   % number of slots (teeth on rotor)
rho = 1.8e-8;               % resistivity (ohm-m) of copper at 25C

%% Coupling Variables b, shared with control problem
n   = 122;                   % rotational speed (rev/s)
V   = 40;                   % design voltage (V)
Pmin= 3.94;                 % minimum required power (kW)
Tmin= 5.12e-3;                % minimum required torque (kNm)

%% Derived parameters and constants
mu0 = 4*pi*1e-7;            % magnetic constant
ap  = p;                    % parallel circuit paths (equals poles)
eff = 0.85;                 % efficiency
bfc  = 40e-3;               % pole depth (m)
fcf  = 0.55;                % field coil spacing factor
Awa  = 2.0074e-006;         % cross sectional area of armature winding (m^2)
Awf  = 0.2749e-006;%3.2419e-005;         % cross sectional area of field coil winding (m^2)

%% Design Calculations
ac   = 2*Pmin*Lwa/( pi*eff*ap*V*( 2*L + 2.3*pi*D/p + 5*ds)*D ); % specific electric loading (kA/m)
Bg   = Pmin/(pi^2*fi*ac*D^2*L*n);   % Maximum Flux Density (T)
Bav  = fi*Bg;                       % average gap density (T)
tau  = pi*D/p;                      % pole pitch
fa   = Bav*tau*L;                   % average flux per pole in armature (Wb)
hf   = 0.76*fa/L;                   % height of field winding (m)
Lmtf = (L+2*bfc+fa/L);              % mean turn length of field coil (m)
I    = Awf^2*V/( hf*Lmtf*bfc*fcf*p*rho );       % field current (A)

%% Control Problem Variables
Ra = rho*Lwa/4/Awa;                 % armature resitance (Ohms)
La = pi*fi*mu0*Lwa^2*D*L/(0.025*pi*D*(2*L+2.3*pi*D/p+5*ds)^2); % armature inductance (H)
Ja = 0.5*0.0020*pi*dfe*L*(D/2-ds)^2 + 0.5*dcu*Lwa*Awa*((D/2)^2+(D/2-ds)^2); % armature inertia (kg m^2)
Km = Tmin*V/Pmin;                   % motor torque constant

v = [La,Ra,Km,Ja];

%% Design Requirements (Equality Constraints) Ceq <= 0
ceq = [];

%% Design Requirements (Inequality Constraints) C <= 0
c(1)  = pi*D/p - 0.38;                  % pole pitch upper bound (m)
c(2)  = 2*D*Bg/(D-2*ds) - 2;          % magnetization margin??
c(3)  = pi*D*n - 25;                    % armature periferal speed UB (m/s)
c(4)  = 8 - pi*D*n;                     % armature periferal speed LB (m/s)
c(5)  = L*p /(pi*D) - 1.0;              % length to pole pitch ratio UB
c(6)  = 0.6 - L*p/(pi*D);               % length to pole pitch ratio LB
c(7)  = 2*D/(D-2*ds) - 3.5;             %
c(8)  = 2.4 - 2*D/(D-2*ds);
c(9)  = Bg - 0.8;                       % maximum flux density UB (T)
c(10) = 0.3 - Bg;                       % maximum flux density LB (T)
c(11) = ds - 2*pi*(D-2*ds)/S;           % geometry limit (m)
c(12) = ds - 0.5*D;                     % slot depth UB (m)
c(13) = ac - 21;                        % electric loading UB (kA)
c(14) = 6 - ac;                         % electric loading LB (kA)
c(15) = Pmin - pi^2*fi*Bg*ac*D^2*L*n;   % minimum power requirement (kW)
c(16) = Tmin - pi/2*fi*Bg*ac*D^2*L;     % minimum torque requirement (Nm)
c(17) = 1.2*Lwa*Awa/( 2*L+2.3*pi*D/p +5*ds ) - pi*( (D/2)^2 - (D/2-ds)^2 );
c(18) = 400 - I^2*Lwf*rho/( Awf*(Lmtf+bfc)*hf );
c(19) = I^2*Lwf*rho/( Awf*(Lmtf+bfc)*hf ) - 850;

%% Define Use (Objective or Constraint)
switch nargout
    case 1 % objective function
        varargout{1} = dcu*(Awa*Lwa + Awf*Lwf)+dfe*L*pi*(D+ds)^2; % weight
    case 2
        varargout{1} = 22.171 - (dcu*(Awa*Lwa + Awf*Lwf)+dfe*L*pi*(0.9*D+ds)^2); % weight
    case 3 % constraints
        varargout{1} = c;
        varargout{2} = ceq;
    case {4,5} % objective, constraints, and control outputs
        varargout{1} = dcu*(Awa*Lwa + Awf*Lwf)+dfe*L*pi*(D+ds)^2;
        varargout{2} = c;
        varargout{3} = ceq;
        varargout{4} = v;
    case 6
        varargout{1}=(D)+0.02*cos(dcu);
    
end






















% make it more expensive
%eigs(magic(200));
