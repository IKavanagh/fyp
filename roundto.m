function r = roundto(n, d)
% Rounds n to d decimal places
r = round(n*(10^d)) / (10^d);