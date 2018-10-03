sigma=10; w=8/3; tau=0.4;
s=10;
syms r A B C D E
SteadyEqn=[-A + r*B - s*D==0;
           -B + A*(1-C)==0;
           -C + A*B==0;
           -tau*D + A*(1-E)==0;
           -tau*E + A*D==0];

%rs=[linspace(1e-10, 1, 20)'; linspace(1,1.245,50)'; linspace(1.245,1.255,100)';...
%   linspace(1.255,1.5,50)'; linspace(1.5, 3, 50)'];
rs=[linspace(1e-10, 8.45, 48)'; linspace(8.45, 8.51, 10)'; linspace(8.51, 8.517,50)';...
    linspace(8.517, 8.533, 10)'; linspace(8.533, 8.85, 20)'; linspace(8.85, 50, 132)'];
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
J=[ -sigma sigma*r 0 -s*sigma 0;
    1-C -1 -A 0 0;
    w*B w*A -w 0 0;
    1-E 0 0 -tau -A;
    w*D 0 0 w*A -w*tau ];
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
subplot(3,2,1);
plot(rs, SteadyAs,'k');
ylabel('A');
xlabel('r');
subplot(3,2,2);
plot(rs, SteadyBs,'k');
ylabel('B');
xlabel('r');
subplot(3,2,3);
plot(rs, SteadyCs,'k');
ylabel('C');
xlabel('r');
subplot(3,2,4);
plot(rs, SteadyDs,'k');
ylabel('D');
xlabel('r');
subplot(3,2,5);
plot(rs, SteadyEs,'k');
ylabel('E');
xlabel('r');


range=find(rs==1);
% Plot eigenvalues
figure;
hold on;
plot(rs,max(real(Jeig(:,:,1)).'), 'k');
for i=2:nsol
plot(rs(range(1):end), max(real(Jeig(range(1):end,:,i)).'), 'k--');
legend('Trivial Solution', 'Branch 1 and 2');
xlabel('r'); ylabel('Maximum real part of eivenvalues')
end
hold off;

% Plot bifurcation
figure;
plot(rs, SteadyAs,'k');
xlabel('r'); ylabel('A');

%% s=0.1 and s=0.2 and s=10
% Range of subcritical Steady States
range=find(SteadyAs(:,2)~=0);
i1=range(1); i2=range(end)+1;

rsFiltered={rs; [rs(i1+1); rs(i1:i2)]; [rs(i1+1); rs(i1:i2)]; rs(i1+1:end); rs(i1+1:end)};
SteadyAsFiltered={SteadyAs(:,1); [SteadyAs(i1+1,4); SteadyAs(i1:i2,2)];...
    [SteadyAs(i1+1,5); SteadyAs(i1:i2,3)]; SteadyAs(i1+1:end,4); SteadyAs(i1+1:end,5)};
SteadyBsFiltered={SteadyBs(:,1); [SteadyBs(i1+1,4); SteadyBs(i1:i2,2)];...
    [SteadyBs(i1+1,5); SteadyBs(i1:i2,3)]; SteadyBs(i1+1:end,4); SteadyBs(i1+1:end,5)};
SteadyCsFiltered={SteadyCs(:,1); [SteadyCs(i1+1,4); SteadyCs(i1:i2,2)];...
    [SteadyCs(i1+1,5); SteadyCs(i1:i2,3)]; SteadyCs(i1+1:end,4); SteadyCs(i1+1:end,5)};
SteadyEsFiltered={SteadyEs(:,1); [SteadyEs(i1+1,4); SteadyEs(i1:i2,2)];...
    [SteadyEs(i1+1,5); SteadyEs(i1:i2,3)]; SteadyEs(i1+1:end,4); SteadyEs(i1+1:end,5)};
SteadyDsFiltered={SteadyDs(:,1); [SteadyDs(i1+1,4); SteadyDs(i1:i2,2)];...
    [SteadyDs(i1+1,5); SteadyDs(i1:i2,3)]; SteadyDs(i1+1:end,4); SteadyDs(i1+1:end,5)};
% Plot Steady States
figure;
subplot(3,2,1);
hold on;
for i=1:nsol
plot(cell2mat(rsFiltered(i)), cell2mat(SteadyAsFiltered(i)));
end
hold off;
xlabel('r'); ylabel('A');
subplot(3,2,2);
hold on;
for i=1:nsol
plot(cell2mat(rsFiltered(i)), cell2mat(SteadyBsFiltered(i)));
end
hold off;
xlabel('r'); ylabel('B');
subplot(3,2,3);
hold on;
for i=1:nsol
plot(cell2mat(rsFiltered(i)), cell2mat(SteadyCsFiltered(i)));
end
hold off;
xlabel('r'); ylabel('C');
subplot(3,2,4);
hold on;
for i=1:nsol
plot(cell2mat(rsFiltered(i)), cell2mat(SteadyDsFiltered(i)));
end
hold off;
xlabel('r'); ylabel('D');
subplot(3,2,5);
hold on;
for i=1:nsol
plot(cell2mat(rsFiltered(i)), cell2mat(SteadyEsFiltered(i)));
end
hold off;
xlabel('r'); ylabel('E');

% Plot eigenvalues
JeigFiltered={Jeig(:,:,1); cat(1, Jeig(i1+1,:,4), Jeig(i1:i2,:,2)); cat(1, Jeig(i+1,:,5), Jeig(i1:i2,:,3)); ...
    Jeig(i1+1:end,:,4); Jeig(i1+1:end,:,5)};

figure;
hold on;
for i=1:nsol
plot(cell2mat(rsFiltered(i)), max(real(cell2mat(JeigFiltered(i))).'));
xlabel('r'); ylabel('Maximum real part of eivenvalues');
end
hold off;

% Plot bifurcation
figure;
hold on;
for i=1:nsol
    plot(cell2mat(rsFiltered(i)), cell2mat(SteadyAsFiltered(i)));
end

xlabel('r'), ylabel('A');

