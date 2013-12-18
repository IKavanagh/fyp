function v = orderof(s)
% Determines the order of the scalar passed in
% orderof(1e3) will return 3, orderof(1e-3) will return -3
v = floor(log(abs(s))./log(10));