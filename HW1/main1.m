sigma=10; w=8/3; tau=0.4;
syms r s A B C D E l
J=[ -sigma sigma*r 0 -s*sigma 0;
    1-C -1 -A 0 0;
    w*B w*A -w 0 0;
    1-E 0 0 -tau -A;
    w*D 0 0 w*A -w*tau ];
J=subs(J,[A B C D E], [0 0 0 0 0]);
character=det(J-eye(5).*l)==0;
l1=-57/5;
l2S=10*(s-r)+72/5;
l1l2S=-(-4*r+10*s+4);
HopfCond1=solve(l1*l2S==l1l2S,r);
HopfCond2=solve(l2S>0,r);

% Test
sTest=10;
rTest=subs(HopfCond1,s,sTest);
eigen=eig(subs(J, [r s], [rTest sTest]));