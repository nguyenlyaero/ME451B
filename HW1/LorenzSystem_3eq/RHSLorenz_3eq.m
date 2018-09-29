function f = RHSLorenz_3eq(t,XYZ)
global sigma r s
%f = zeros(3,1);
f(1,1) = -sigma*XYZ(1) + sigma*XYZ(2);
f(2,1) = r*XYZ(1) - XYZ(2) - XYZ(1)*XYZ(3);
f(3,1) = -s*XYZ(3) + XYZ(1)*XYZ(2);