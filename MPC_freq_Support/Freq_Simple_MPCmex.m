N = 50;
Ts = 0.02;

M = 4.0;
D = 1.5;
Rp = 0.05;
Tg = 0.2;

import casadi.*;
opti = casadi.Opti();

% State variable
w = opti.variable(N + 1);
dotw = opti.variable(N + 1);

% Control input
u = opti.variable(N);

% Initial states
w0p = opti.parameter();
dotw0p = opti.parameter();
Pdp = opti.parameter();
Q11 = opti.parameter();
Q22 = opti.parameter();
R11 = opti.parameter();

% Cost function
J = sum(Q11 * (w - 0).^2) + sum(Q22 * dotw.^2) + sum(R11 * u.^2);
opti.minimize(J);

for k = 1:N
    opti.subject_to(w(k + 1) == w(k) + Ts * (dotw(k + 1)));
    opti.subject_to(dotw(k + 1) == dotw(k) + Ts * (- (D / (M * Tg) + 1 / (Rp * M * Tg)) * w(k + 1) - (D / M + 1 / Tg) * dotw(k + 1) - 1 / (M * Tg) * (Pdp - u(k))));
end

opti.subject_to(w(1) == w0p);
opti.subject_to(dotw(1) == dotw0p);

for k=1:N
    opti.subject_to(-1 <= u(k) <= 1);
end
% for k=2:N
%     opti.subject_to(-10*Ts <= u(k) - u(k-1) <= 10*Ts);
% end

opts.qpsol = 'qrqp';
opti.solver('sqpmethod', opts);

% opti.solver('bonmin', struct('expand',false), struct('max_iter',100));
M = opti.to_function('M', {w0p, dotw0p, Pdp, Q11, Q22, R11}, {u(1), J});

M.generate('Mc', struct('mex', true));
mex Mc.c -largeArrayDims
