sigma=10; w=8/3; tau=0.4;
s=0;
syms r A B C D E
SteadyEqn=[-A + r*B - s*D==0;
           -B + A*(1-C)==0;
           -C + A*B==0;
           -tau*D + A*(1-E)==0;
           -tau*E + A*D==0];

rs=[linspace(1e-10, 1, 20)'; linspace(1, 3, 100)'];
% Steady-States

SteadySol=solve(SteadyEqn, A, B, C, D, E);

nsol=length(SteadySol.A);
SteadyAs=zeros(120,nsol);
SteadyBs=zeros(120,nsol);
SteadyCs=zeros(120,nsol);
SteadyDs=zeros(120,nsol);
SteadyEs=zeros(120,nsol);
for i=1:120
    SteadyAs(i,:)=subs(SteadySol.A, r, rs(i));
    SteadyBs(i,:)=subs(SteadySol.B, r, rs(i));
    SteadyCs(i,:)=subs(SteadySol.C, r, rs(i));
    SteadyDs(i,:)=subs(SteadySol.D, r, rs(i));
    SteadyEs(i,:)=subs(SteadySol.E, r, rs(i));
end

AllRealFilter=(imag(SteadyAs)==0).*(imag(SteadyBs)==0)...
    .*(imag(SteadyCs)==0).*(imag(SteadyDs)==0).*(imag(SteadyEs)==0);
SteadyAs=AllRealFilter.*SteadyAs;
SteadyBs=AllRealFilter.*SteadyBs;
SteadyCs=AllRealFilter.*SteadyCs;
SteadyDs=AllRealFilter.*SteadyDs;
SteadyEs=AllRealFilter.*SteadyEs;


figure;
hold on;
plot(rs, AllRealFilter.*SteadyAs,'k');
plot(rs, AllRealFilter.*SteadyBs,'b');
plot(rs, AllRealFilter.*SteadyCs,'r');
plot(rs, AllRealFilter.*SteadyDs, 'y');
plot(rs, AllRealFilter.*SteadyEs, 'm');
title('Steady States');
hold off;

% Stability
J=[ -1 r 0 -s 0;
    1-C -1 -A 0 0;
    B A -1 0 0;
    1-E 0 0 -tau -A;
    D 0 0 A -tau ];
Jeig=zeros(120,5,nsol);

for i=1:nsol
    for j=1:120
        Jeig(j,:,i)=  eig(subs(J, [A, B, C, D, E, r],...
            [SteadyAs(j,i),SteadyBs(j,i),SteadyCs(j,i),SteadyDs(j,i),SteadyEs(j,i),rs(j)]));
    end
end

% Plot eigenvalues
figure;
hold on;
for i=1:nsol
plot(rs, max(real(Jeig(:,:,i)).'));
end
hold off;