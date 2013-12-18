function v = isposdef(m)
% Tests if the matrix passed into the function is postive definite,
% 1 == true
[~, p] = chol(m);
v = (p == 0);