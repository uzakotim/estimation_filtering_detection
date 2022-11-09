% matric analysis
%
% derivatives of x'Ay, x'Ax
%
syms x1 x2 y1 y2 y3 A
%
x = [x1 x2].'
y = [y1 y2 y3].'
A = [1 2 3; 4 5 6]
input('Hit any key to continue ...')
%
xAy = x.'*A*y
%
disp('Compare d./dy')
d_dy = [diff(xAy,y1);
        diff(xAy,y2);
        diff(xAy,y3)]
ATx = A.'*x
input('Hit any key to continue ...')
%
disp('Compare d./dx')
d_dx = [diff(xAy,x1);
        diff(xAy,x2)]
Ay = A*y
input('Hit any key to continue ...')
