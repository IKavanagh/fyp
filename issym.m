function v = issym(m)
% Tests if the matrix passed into the function is symmetrical, 1 == true
test = @(x) isequal(x,x');
v = test(m);