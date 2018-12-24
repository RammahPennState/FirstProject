% This function is to solve delay differential equations (DDEs).
% If there are d equations, y is a d×1 vector that approximates the soltion y(t).
% Z is a d × k array. Z(:,j) approximates y(t ? tuaj).
function Dydt = comp1DDE(t,y,Z,p,ts,Uin)
% The parameters
V1 = p.V1; V2 = p.V2; V3 = p.V3; T1 = p.T1; T2 = p.T2; E = p.E; td = p.td;
a1 = p.a1; C1 = p.C1; beta = p.beta; C4 = p.C4; U0 = p.U0; Um = p.Um; Rg = p.Rg;
alpha = p.alpha; C5 = p.C5; Rm = p.Rm;
% The states
Ip = y(1);
Ii = y(2);
G  = y(3);
%
% time lags
tua = Z(:,1); % this is a (3x1) vector represents: row {Ip,Ii,G} (states), and column delay 'tua' (min).
% Ip(t-tua) will be tua(1,1), which is the first state in the DDEs.
% There are 3 ddes questions.
dydt = zeros(3,1);
% ****************
% The feeding function is converted to continuous-time.
u_external = interp1(ts,Uin,t);
% *****************************
% the 5 rate functions, which use the concentrations.
cx = Ip/V1; % MU/ml
cy = Ii/V2; % MU/ml
cz = 0.1*G/V3; % Glucose concentration [mg/dl]
% Ip with the delay tua.
cx_tua = tua(1)/V1;
% the function f1
a = C1/a1;
b = -0.03; % comes form (1/30).
% Note that in f1: exp([C1/a1 - G/Vg*a])= exp(a + b*cz), where cz = 0.1*G/Vg , b = 1/30. 
f1 = Rm/(1+ exp(a+b*cz)); % [mU/min]
% Note in f2: I*(E*ti + V/Vi*E*ti) = (I/Vi)(E*ti + Vi/E*ti)= cy*[1 + Vi/E*ti] where
% cy = I/Vi.
f2 = ((Um-U0)/(1+ exp(beta*log(C4) - beta*log(cy*(1+V2/(E*T2)))))+ U0);% [mg/min dl]

f3 = Rg/(1+ exp(alpha*cx_tua-alpha*C5));
% *****************************
% DDEs 
dydt(1) = f1 - ((E/V1)+(1/T1))*Ip + (E/V2)*Ii; % [mU/min]
dydt(2) = (E/V1)*Ip -((E/V2)+(1/T2))*Ii; % [mU/min]
dydt(3) = f3 - cz *f2 + u_external;% [mg/min]
Dydt = [dydt(1); dydt(2); dydt(3)];
end
