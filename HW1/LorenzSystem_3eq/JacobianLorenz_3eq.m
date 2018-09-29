function e = JacobianLorenz_3eq(sigma,r,s,x0,y0,z0)
dfidxj = [-sigma,sigma,0;r-z0,-1,-x0;y0,x0,-s];
e = eig(dfidxj);