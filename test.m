% f = 2


% c = 1
% rho = 1
% u0 = 0
% 'DD', 0, 0

%fem1d (20, 'muniform', 'DD', [0 0], 'f', 'c', 'rho', 0.1, 1, 'trapezoid', 'u0', 'eulerImplicit')
%fem1d (20, 'mquadratic', 'DD', [0 0], 'f', 'c', 'rho', 0.1, 1, 'trapezoid', 'u0', 'eulerImplicit')
%fem1d (20, 'mrandom', 'DD', [0 0], 'f', 'c', 'rho', 0.1, 1, 'trapezoid', 'u0', 'eulerImplicit')


% f = 2
% c = 1
% rho = 1
% u0 = 0
% 'DN', 0, 0

%fem1d (20, 'muniform', 'DN', [0 0], 'f', 'c', 'rho', 0.1, 1, 'trapezoid', 'u0', 'eulerImplicit')
%fem1d (20, 'mquadratic', 'DN', [0 0], 'f', 'c', 'rho', 0.1, 1, 'trapezoid', 'u0', 'eulerImplicit')
%fem1d (20, 'mrandom', 'DN', [0 0], 'f', 'c', 'rho', 0.1, 1, 'trapezoid', 'u0', 'eulerImplicit')



% f = 10
% c = 1
% rho = 1
% u0 = -x
% 'DN', 0, -1

fem1d (20, 'muniform', 'DN', [0 -1], 'f', 'c', 'rho', 0.1, 1, 'trapezoid', 'u0', 'eulerImplicit')
fem1d (20, 'mquadratic', 'DN', [0 -1], 'f', 'c', 'rho', 0.1, 1, 'trapezoid', 'u0', 'eulerImplicit')
fem1d (20, 'mrandom', 'DN', [0 -1], 'f', 'c', 'rho', 0.1, 1, 'trapezoid', 'u0', 'eulerImplicit')