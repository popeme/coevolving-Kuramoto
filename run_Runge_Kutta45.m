%% new Runge-Kutta integrator-- trying it

load sch200_SC.mat
master_Kij = Kij;

% Initialize Parameters
N = size(Kij,1);                            % number of nodes
frequency_mean = 40;                        % mean intrinsic frequency (Hz)
f_std = 0.1;                                % freq std

% Time Parameters
runlength = 792+20;                         % runlength used for all analyses (20=length of HRF in sec, equals transient later discarded)
transient = 20;
dt = 0.001;                                 % size of integration step (do not change), in seconds
h = dt;
tspan = single([0:dt:runlength]);                  % solve diff equations at these time points
tspan_transient = [0:dt:transient];        % solve diff equations at these time points (transient)
step = 720;                                 % TR in msec, used for binning Ybold and R
lts = 1100;                                 % length of run in TRs

% Delay/Frustration Matrix
vel = 12;                                   % conduction velocity in mm/ms
D2 = D.*(1/vel);                            % delay in milliseconds
Aij_master = D2./((1/dt)/frequency_mean).*(2*pi);  % phase frustration (expresses delay in fractions of a cycle of the average oscillator)

%Scaled Coupling
kval = 280;

%Nearest Neighbors
m = 20;

% initial condition (phase) and intrinsic phases
ths_ic = (2*pi).*rand(N,1);
omegas = (2*pi).*(frequency_mean+f_std.*randn(N,1));

%set up variables to store
ths = zeros(numel(tspan),N);
ths_test = zeros(numel(tspan_transient),N);
edges_all = logical(zeros(numel(tspan),6040));

% Integration of transient
ths_test(1,:) = ths_ic;
for ii = 1:numel(tspan_transient)-1
    [Kij,~] = edge_selection(ths_test(ii,:)',master_Kij,m,'Most');
    Aij = zeros(200);
    Aij(Kij>0) = Aij_master(Kij>0);
    [y_f, out, h_next] = Runge_Kutta_Fehlberg_4_5_KC(@KS_fcn_RK,ths_test(ii,:)',tspan_transient(ii),h,tspan_transient(ii+1),1e-6, Kij, N, kval, Aij, omegas);
	h = h_next;
    ths_test(ii+1,:) = y_f;
    disp(ii)
end

%profile on
%Integration of whole run, saving
ths(1,:) = ths_test(end,:);
[Kij,edges_all(1,:)] = edge_selection(ths(1,:)',master_Kij,m,'Most');
Aij = zeros(200);
Aij(Kij>0) = Aij(Kij>0);
for ii = 1:numel(tspan)-1
    [y_f, out, h_next] = Runge_Kutta_Fehlberg_4_5_KC(@KS_fcn_RK,ths(ii,:)',tspan(ii),h,tspan(ii+1),1e-6, Kij, N, kval, Aij, omegas);
	h = h_next;
    ths(ii+1,:) = y_f;
    [Kij,edges_all(ii+1,:)] = edge_selection(ths(ii,:)',master_Kij,m,'Most');
    Aij = zeros(200);
    Aij(Kij>0) = Aij_master(Kij>0);
    if mod(ii,1000) == 0
        disp(ii/(numel(tspan)-1)*100)
    end
%profile info
end

save('test2.mat','ths','edges_all','-v7.3')
%%
% Convert to BOLD
load test2.mat
load sch200_FC.mat
Ybold_reg = convert_to_BOLD(ths(2:end,:)',lts,step,N);

yboldz = zscore(Ybold_reg);
FC = corr(yboldz);
FCfish = fcn_fisher(FC);
nFCfish = fcn_fisher(nFCavg);

inds = logical(triu(ones(200),1));

r = corr(nFCfish(inds),FCfish(inds));
r2 = corr(FCnorm_avg(inds),FCfish(inds));
