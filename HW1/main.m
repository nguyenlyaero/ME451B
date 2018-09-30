sigma=10; w=8/3; tau=0.4;
s=0;
syms r A B C D E
SteadyEqn=[-A + r*B - s*D==0;
           -B + A*(1-C)==0;
           -C + A*B==0;
           -tau*D + A*(1-E)==0;
           -tau*E + A*D==0];

rs=[linspace(1e-10, 1, 20)'; linspace(1,1.245,50)'; linspace(1.245,1.255,100)'; linspace(1.255,1.5,50)'; linspace(1.5, 3, 50)'];
% Steady-States

SteadySol=solve(SteadyEqn, A, B, C, D, E);

nsol=length(SteadySol.A);
SteadyAs=zeros(270,nsol);
SteadyBs=zeros(270,nsol);
SteadyCs=zeros(270,nsol);
SteadyDs=zeros(270,nsol);
SteadyEs=zeros(270,nsol);
for i=1:270
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

% Stability
J=[ -1 r 0 -s 0;
    1-C -1 -A 0 0;
    B A -1 0 0;
    1-E 0 0 -tau -A;
    D 0 0 A -tau ];
Jeig=zeros(270,5,nsol); % r,eig,sol

for i=1:nsol
    for j=1:270
        Jeig(j,:,i)=  eig(subs(J, [A, B, C, D, E, r],...
            [SteadyAs(j,i),SteadyBs(j,i),SteadyCs(j,i),SteadyDs(j,i),SteadyEs(j,i),rs(j)]));
    end
end

%% s=0
% Plot Steady States
figure;
hold on;
plot(rs, SteadyAs,'k');
plot(rs, SteadyBs,'b');
plot(rs, SteadyCs,'r');
plot(rs, SteadyDs,'y');
plot(rs, SteadyEs,'m');
title('Steady States');
hold off;


% Plot eigenvalues
figure;
hold on;
for i=1:nsol
plot(rs, max(real(Jeig(:,:,i)).'));
end
hold off;

%% s=0.1 and s=0.2
% Range of subcritical Steady States
range=find(SteadyAs(:,2)~=0);
i1=range(1); i2=range(end)+1;

rsFiltered={rs; rs(i1:i2); rs(i1:i2); [rs(i1); rs(i2:end)]; [rs(i1);rs(i2:end)]};
SteadyAsFiltered={SteadyAs(:,1); SteadyAs(i1:i2,2);...
    SteadyAs(i1:i2,3); [SteadyAs(i1,2); SteadyAs(i2:end,4)]; [SteadyAs(i1,3); SteadyAs(i2:end,5)]};
SteadyBsFiltered={SteadyBs(:,1); SteadyBs(i1:i2,2);...
    SteadyBs(i1:i2,3); [SteadyBs(i1,2); SteadyBs(i2:end,4)]; [SteadyBs(i1,3); SteadyBs(i2:end,5)]};
SteadyCsFiltered={SteadyCs(:,1); SteadyCs(i1:i2,2);...
    SteadyCs(i1:i2,3); [SteadyCs(i1,2); SteadyCs(i2:end,4)]; [SteadyCs(i1,3); SteadyCs(i2:end,5)]};
SteadyDsFiltered={SteadyDs(:,1); SteadyDs(i1:i2,2);...
    SteadyDs(i1:i2,3); [SteadyDs(i1,2); SteadyDs(i2:end,4)]; [SteadyDs(i1,3); SteadyDs(i2:end,5)]};
SteadyEsFiltered={SteadyEs(:,1); SteadyEs(i1:i2,2);...
    SteadyEs(i1:i2,3); [SteadyEs(i1,2); SteadyEs(i2:end,4)]; [SteadyEs(i1,3); SteadyEs(i2:end,5)]};
% Plot Steady States
figure;
hold on;
for i=1:5
plot(cell2mat(rsFiltered(i)), cell2mat(SteadyCsFiltered(i)));
end
hold off;

% Plot eigenvalues
JeigFiltered={Jeig(:,:,1); Jeig(i1:i2,:,2); Jeig(i1:i2,:,3); ...
    [Jeig(i1,:,2); Jeig(i2:end,:,4)]; [Jeig(i1,:,3); Jeig(i2:end,:,4)]};

figure;
hold on;
for i=1:nsol
plot(cell2mat(rsFiltered(i)), max(real(cell2mat(JeigFiltered(i))).'));
end
hold off;
