function [x0,y0,z0] = FixedPointsLorenz_3eq(s,r)
x0(1) = 0; y0(1) = 0; z0(1) = 0;
if (r>1)
    z0(2:3) = r - 1;
    x0(2) = sqrt(s)*sqrt(z0(2));
    x0(3) = -sqrt(s)*sqrt(z0(3));
    y0(2) = sqrt(s)*sqrt(z0(2));
    y0(3) = -sqrt(s)*sqrt(z0(3));
end