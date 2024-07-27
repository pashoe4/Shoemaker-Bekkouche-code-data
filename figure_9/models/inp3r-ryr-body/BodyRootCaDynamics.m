function out=BodyRootCaDynamics(mpfile, varargin)
%%%% BodyRootCaDynamics.m %%%%
% Compute equilibrium state of a compartment with ER + RYRs, followed by
% time-domain response of an annular cell body model with ER,
% with 2-D (circumferential & radial) diffusion,
% and junctions with a set of cellular processes,
% with defined 'root' regions that have different ER-related variables.
% 'Excitation' takes form of incoming calcium waves from the processes
% InP3 model is Siekmann-Cao-Sneyd model with re-defined CI3 state;
% RyR model is from Breit & Quiesser / Keizer & Levine.
% RyRs may be discretely distributed along dendrite, < 1/compartment
%
% User is prompted for a file containing desired parameter values.
% This file contains the name of a parameter file for the
% associated cellular processes, including roots
p = inputParser;
addOptional(p,'ker',0);
addOptional(p,'tuning_only',false);
addOptional(p,'tend_sim',1.0);% Simulation time
parse(p,varargin{:});
ker=p.Results.ker;%Set virtual calcium source for ER
tuning_only=p.Results.tuning_only; % Enabels early return for ker-tuning
tend_sim = p.Results.tend_sim;
out='';
dt = 1E-5;   % time step for temporal integration (units s)
%dt = 4E-6;   % small time step for temporal integration (units s)

iplot = 200; % plot interval in units of dt
tout = iplot*dt; % plot interval

%% Load parameters

% load model parameters
%mpfile = input('Specify name of model parameter file: ','s');
load(mpfile); % cell body model parameters

% load rate constants / function tables for InP3 receptors
load('InP3R_rates.mat');

% load rate constants / function tables for Ryanodine receptors
load('RyR_rates.mat');

%% Set virtual calcium source for ER
%ker = 0; % for default value RrI=10; applies to processes elements
%ker = -2.9; % for RrI=5; applies to processes elements
%ker = 6.9; %for RrI=20; applies to processes elements
kerb = 2*ker; % applies to body & root elements

%% Set distal Ca & InP3 influx parameters for activated processes

% influx rate
Cain = 200; % units uM.s^-1
InP3in = 150; % units uM.s^-1

BVCain = BCc * Cain;
BVInP3in = BInP3 * InP3in;

% influx pulse duration
tin = 0.02; % units s

%% Scale rate constants by dt for purposes of temporal integration

%%% Universal calcium parameters %%%

% Virtual calcium sources for ER
DCer = dt*ker; % processes (note virtual source assumed constant)
DCerb = dt*kerb; % body (note virtual source assumed constant)

% Calcium molar leakage rates
%   (note extracellular leakage flux assumed constant)
DClx = dt*klx;
DClxbf = dt*klxbf;
DClxbw = dt*klxbw;
%   ER membrane leakage
kle = dt*kle;
kler = dt*kler;
kleb = dt*kleb;

% Cytosol Ca diffusion coefficient
DCac = dt*DCac;
DCacr = dt*DCacr;
% ER Ca diffusion coefficient
DCae = dt*DCae;
DCaer = dt*DCaer;

% calcium buffering rates
kbCf = dt*kbCf;
kbCb = dt*kbCb;

%%% Calcium pump parameters %%%

% pump clearance rate constants
kpN = dt*kpN;
kpNbf = dt*kpNbf;
kpNbw = dt*kpNbw;
kpP = dt*kpP;
kpPbf = dt*kpPbf;
kpPbw = dt*kpPbw;
kpS = dt*kpS;
kpSr = dt*kpSr;
kpSb =  dt*kpSb;

%%% InP3 receptor and related parameters %%%

% InP3 diffusion coefficient
DInP3 = dt*DInP3;
DInP3r = dt*DInP3r;

% Calcium molar influx rate for open InP3R channels
krI = dt*krI;
krIr = dt*krIr;
krIb = dt*krIb;

% InP3 generation rate
kmIf = dt*kmIf;
kmIfbf = dt*kmIfbf;
kmIfbw = dt*kmIfbw;
% InP3 degradation rate
kmIb = dt*kmIb;

% InP3 receptor state rates
rI12 = dt*rI12;  
rI21 = dt*rI21;
rI26 = dt*rI26;
rI62 = dt*rI62;
rI45 = dt*rI45;
rI54 = dt*rI54;
rI34 = dt*rI34;

rI43 = dt*rI43;
rI42 = dt*rI42;
rI24 = dt*rI24;

% rate parameter for InP3 receptor gating
rIg = rIg*dt;

% PKC activation rate
kiPf = dt*kiPf;
% PKC deactivation rate
kiPb = dt*kiPb;

%%% Ryanodine receptor parameters %%%

% Calcium molar influx rate for open RyR channels
krR = dt*krR;

% RyR receptor state rates
rR23 = dt*rR23;
rR32 = dt*rR32;
rR31 = dt*rR31;
rR43 = dt*rR43;

kR13 = dt*kR13;
kR34 = dt*kR34;

rsat13 = dt*rsat13;
rsat34 = dt*rsat34;

%% Set up initial conditions / initial states for equilibrium calculation

% Small calcium concentrations used to ensure convergence to equilibrium
Cc = 0.01; % cytosolic Ca; units uM
Ce = 200; % ER Ca; units uM

CaCalB = Cc*kbCf*CalB0/kbCb; % bound calcium value

% approximate quiescent InP3-related states
InP3 = kmIf*Cc/kmIb; % cytosolic InP3;
CI1 = 0;
CI2 = 0;
CI3 = 0;
CI4 = 1;
OI5 = 0;
OI6 = 0;
PKCv = 0;
% gating variables
mIt = 0;
hIt = 1;

% allocate InP3R state increments
DI21 = 0;
DI26 = 0;
DI42 = 0;
DI45 = 0;
DI43 = 0;

% approximate quiescent RyR states
CR1 = 1;
CR2 = 0;
OR3 = 0;
OR4 = 0;

% allocate RyR state increments
DR13 = 0;
DR32 = 0;
DR34 = 0;

% allocate InP3R calcium increment
DCI = 0;

% allocate RyR calcium increments
DCR = 0;

%% Relaxation to observe (unexcited) equilibrium states, or instability

% output variables to examine history
Cct=[]; Cet = []; time = [];

DClx0 = volRoc*DClx; % rescaled external leakage for single compartment
tend = 280; % end time
for t = 0:dt:tend % loop on time
    
    % Compute state increments
    
    % Ca leakage from ER
    DCle = volRic*kle.*(Ce-Cc);
    
    % Ca buffering in cytosol
    DCfb = kbCf.*Cc.*(CalB0-CaCalB);
    DCbf = kbCb.*CaCalB;
    
    % Calcium pumps
    
        % NCX pumps
    DCN = volRoc*kpN.*Cc./(kpNC+ Cc);
    
    % PMCA pumps
    Cc2pw =  Cc.* Cc;
    DCP = volRoc*kpP.*Cc2pw./(kpPC2pw+Cc2pw);
    
    % SERCA pumps
    DCS = volRic*kpS.* Cc./((kpSC+ Cc).*Ce);
    
    % InP3R and related states
    
    % Gating function dynamics
    Cc6pw = Cc.^6;
    mI =  Cc6pw./(Cc6pw+km6pw);
    hI =  kh6pw./(Cc6pw+kh6pw);
    
    In3pw = InP3.^3;
    mIp = In3pw./(In3pw+kbp3pw);
    
%     % No dynamics; instantaneous gating
%     mIt = mI;
%     hIt = hI;

    % single-pole delay
        % (done here since these play role in InPrR state updates)
    mIt = mIt + rIg.*( mI - mIt );
    hIt = hIt + rIg.*( hI - hIt );

    % InP3 states increments
    phI43 = rI43.*(1 - hIt);
    phGate = mIp.*mIt.*hIt;
    phI42 = rI42.*(    phGate);
    phI24 = rI24.*(1 - phGate);
    DI21 = rI21.*CI2 - rI12.*CI1;
    DI26 = rI26.*CI2 - rI62.*OI6;
    DI42 = phI42.*CI4 - phI24.*CI2;
    DI45 = rI45.*CI4 - rI54.*OI5;
    DI43 = phI43.*CI4 - rI34.*CI3;
    
    % probability of open state for InP3Rs
    OpenI = OI5 + OI6;
    % Ca increment through open InP3R channels
    DCI = volRic*krI*OpenI.*( Ce - Cc );
    
    % cytosolic InP3 increment
    DInP3c = (volRoc*kmIf* Cc)./(1+kiPI*PKCv) - kmIb*InP3;
    
    % PKCv update
    DPKCv = kiPf* Cc.*(PKC0-PKCv) - kiPb*PKCv;
    
    % RyR states
    
    rR13 = kR13.*Cc.^4;
    rR34 = kR34.*Cc.^3;
    phR13 = rR13.*rsat13./(rR13+rsat13);
    phR34 = rR34.*rsat34./(rR34+rsat34);

    % RyR states increments
    DR13 = phR13.*CR1 - rR31.*OR3;
    DR32 = rR32.*OR3 - rR23.*CR2;
    DR34 = phR34.*OR3 - rR43.*OR4;
    
    % probability of open state for RyR
    OpenR = OR3 + OR4;
    
    % Ca increment through open RyR channels
    if Drflag % if there is more than 1 RyR per compartment
        DCR = volRic*krR*OpenR.*( Ce - Cc );
    else % if there are multiple compartments for each RyR
        DCR = volRic*krR*OpenR.*( Ce - Cc )/MrR;
    end
        
    % Apply increments to states
    
    % sum of ER calcium increments
    DCe = ( DCS - DCR-DCI - DCle )/volR + DCer;
    % update Ce
    Ce = Ce + DCe;
    
    % sum of cytosolic calcium increment
    DCc = DCbf-DCfb - DCN-DCP-DCS + DCI+DCR + DClx0+DCle;
    % update Cc
    Cc = Cc + DCc;
    
    % InP3 concentration increment
    InP3 = InP3+DInP3c;
    
    % Calbindin state
    CaCalB = CaCalB + DCfb - DCbf;
    
    % InP3R and related states
    CI1 = CI1 + DI21;
    CI2 = CI2 + DI42 - DI21 - DI26;
    CI3 = CI3 + DI43;
    CI4 = CI4 - DI45 - DI42 - DI43;
    OI5 = OI5 + DI45;
    OI6 = OI6 + DI26;
    
    PKCv = PKCv + DPKCv;
    
    % RyR states
    CR1 = CR1 - DR13;
    CR2 = CR2 + DR32;
    OR3 = OR3 + DR13 - DR32 - DR34;
    OR4 = OR4 + DR34;
    
    if ~mod(t,tout)
        time = [time;t];
        Cct = [Cct,Cc];
        Cet = [Cet,Ce];
    end
    
end

if tuning_only
    out=Ce;
    disp("Tuning iteration completed. Returning.")
    return
end

figure(11)
plot(time, Cct,'b-')
grid on
xlabel('Time (s)');
ylabel('Concentration (uM)');
set(gca,'fontsize',15);
figure(12)
plot(time, Cet,'-','color', [0,0.3,0.6])
grid on
xlabel('Time (s)');
ylabel('Concentration (uM)');
set(gca,'fontsize',15);

Cc_Cmprt = Cc;
Ce_Cmprt = Ce;
InP3_Cmprt = InP3;

%% Set up initial conditions / initial states, cellular processes

ngt = ng + ngr; % total # grid points, process proper + root segment

for k=1:np % loop over number of processes joining cell body
    
    Ccp{k} = Cc*ones(ngt+2,3); % cytosolic Ca
    Cep{k} = Ce*ones(ngt+2,1); % ER Ca
    InP3p{k} = InP3*ones(ngt+2,3); % cytosolic InP3;
    CaCalBp{k} = CaCalB*ones(ngt,3); % bound calcium value
    
    % quiescent InP3R states
    CI1p{k} = CI1*ones(ngt,1);
    CI2p{k} = CI2*ones(ngt,1);
    CI3p{k} = CI3*ones(ngt,1);
    CI4p{k} = CI4*ones(ngt,1);
    OI5p{k} = OI5*ones(ngt,1);
    OI6p{k} = OI6*ones(ngt,1);
    % gating variables
    mItp{k} = mIt*ones(ngt,1);
    hItp{k} = hIt*ones(ngt,1);
    
    % activated PKC
    PKCvp{k} = PKCv*ones(ngt,1);
    
    % quiescent RyR  states -- defined only in process proper
    if Drflag % if >= 1 RyR/compartment, initiate states @ each compartment
        CR1p{k} = CR1*ones(ng,1);
        CR2p{k} = CR2*ones(ng,1);
        OR3p{k} = OR3*ones(ng,1);
        OR4p{k} = OR4*ones(ng,1);
    else % if < 1 RyR/compartment, initiate @ compartments where present
        CR1p{k} = CR1*IRyR;
        CR2p{k} = CR2*IRyR;
        OR3p{k} = OR3*IRyR;
        OR4p{k} = OR4*IRyR;
    end
    
end

% allocate array for ER leak increments
DCle = zeros(ngt,1);

% allocate arrays for pump calcium increments
DCS = zeros(ngt,1);
DCP = zeros(ngt,1);
DCN = zeros(ngt,1);

% allocate arrays for InP3R state increments
DI21 = zeros(ngt,1);
DI26 = zeros(ngt,1);
DI42 = zeros(ngt,1);
DI45 = zeros(ngt,1);
DI43 = zeros(ngt,1);

% allocate arrays for RyR state increments
DR13 = zeros(ng,1);
DR32 = zeros(ng,1);
DR34 = zeros(ng,1);

% allocate arrays for receptor calcium increments
DCI = zeros(ngt,1);
DCR = zeros(ngt,1);

% divisors for joint boundary value calculation at root-body junctions
denCc = -1 / ( DCac*JPb + DCacr*JBb );
denCe = -1 / ( DCae*JPb + DCaer*JBb );
denInP3 = -1 / ( DInP3*JPb + DInP3r*JBb );

%% Set up initial conditions / initial states, cell body

Ccb = Cc*ones(nc,nr,nt); % cytosolic Ca
Ceb = Ce*ones(nc,nr,nt); % ER Ca
InP3b = InP3*ones(nc,nr,nt); % cytosolic InP3;
CaCalBb = CaCalB*ones(nc,nr,nt); % bound calcium value

% quiescent InP3R states
CI1b = CI1*ones(nc,nr,nt);
CI2b = CI2*ones(nc,nr,nt);
CI3b = CI3*ones(nc,nr,nt);
CI4b = CI4*ones(nc,nr,nt);
OI5b = OI5*ones(nc,nr,nt);
OI6b = OI6*ones(nc,nr,nt);
% gating variables
mItb = mIt*ones(nc,nr,nt);
hItb = hIt*ones(nc,nr,nt);

% activated PKC
PKCvb = PKCv*ones(nc,nr,nt);

% allocate array for ER leak increments
DCleb = zeros(nc,nr,nt);

% allocate arrays for pump calcium increments
DCSb = zeros(nc,nr,nt);
DCPb = zeros(nc,nr,nt);
DCNb = zeros(nc,nr,nt);
DInP3Mib = zeros(nc,nr,nt);

% allocate arrays for InP3R state increments
DI21b = zeros(nc,nr,nt);
DI26b = zeros(nc,nr,nt);
DI42b = zeros(nc,nr,nt);
DI45b = zeros(nc,nr,nt);
DI43b = zeros(nc,nr,nt);

% allocate arrays for receptor calcium increments
DCIb = zeros(nc,nr,nt);
% DCRb = zeros(nc,nr,nt);

% % allocate arrays for circumferential diffusion calcium increments
% DCcbC = zeros(nc,nr,nt);
% DCebC = zeros(nc,nr,nt);
% DInP3bC = zeros(nc,nr,nt);

%% Set up initial boundary values, cell body / process junctions

% joint boundary values
BVCc = Cc*ones(1,np);
BVCe = Ce*ones(1,np);
BVInP3 = InP3*ones(1,np);

% %% Optional graphing parameters for processes
% 
% xgrid1 = [0, 0.5*dx:dx:(ngt-0.5)*dx, ngt*dx]';
% %xgrid2 = [0.5*dx:dx:(ngt-0.5)*dx]';
% xlims = [0,(ld+lr)];
% xtic = [0:(ld+lr)/10:(ld+lr)];

%% Loop on time

%tend = 1.0; % end time
tend=tend_sim;
iplot = 100; % plot interval in units of dt
tout = iplot*dt; % plot interval
it = 0;
data_time=[];
for t = 0:dt:tend % loop on time
    
    %% CELLULAR PROCESS EVOLUTION
    
    for k=1:np % loop over processes attached to cell body
        
        % Get calcium and InP3 values
        Cc = Ccp{k};
        Ce = Cep{k};
        InP3 = InP3p{k};
        
        Cc(Cc<0) = 0; % added to prevent numerical undershoot in [Ca]
        InP3(InP3<0) = 0; % added to prevent numerical undershoot in [InP3]
        
        CcM = Cc(2:end-1,:);
        CeM = Ce(2:end-1,1);
        InP3M = InP3(2:end-1,:);
        
        % Compute increments for all states, cellular processes
        
        % Ca leakage from ER
        
        DCle(kpr) = kle.*(CeM(kpr)-CcM(kpr,1));
        DCle(krt) = kler.*(CeM(krt)-CcM(krt,1));
        
        % Ca buffering in cytosol
        
        DCfb = kbCf.*CcM.*(CalB0-CaCalBp{k});
        DCbf = kbCb.*CaCalBp{k};
        
        % Calcium pumps
        
        % NCX pumps
        DCN = kpN.*CcM(:,3)./(kpNC+CcM(:,3));
        
        % PMCA pumps
        Cc2pw = CcM(:,3).*CcM(:,3);
        DCP = kpP.*Cc2pw./(kpPC2pw+Cc2pw);
        
        % SERCA pumps
        %   process proper
        DCS(kpr) = kpS.*CcM(kpr,1)./((kpSC+CcM(kpr,1)).*CeM(kpr));
        %   root segment
        DCS(krt) = kpSr.*CcM(krt,1)./((kpSC+CcM(krt,1)).*CeM(krt));
        
        % InP3R and related states
        
        % Gating function dynamics
        Cc6pw = CcM(:,1).^6;
        mI =  Cc6pw./(Cc6pw+km6pw);
        hI =  kh6pw./(Cc6pw+kh6pw);
        
        In3pw = InP3M(:,1).^3;
        mIp = In3pw./(In3pw+kbp3pw);
        
        %     % No dynamics; instantaneous gating
        %     mItp{k} = mI;
        %     hItp{k} = hI;
        
        % single-pole delay
        % (done here since these play role in InPrR state updates)
        mItp{k} = mItp{k} + rIg.*( mI - mItp{k} );
        hItp{k} = hItp{k} + rIg.*( hI - hItp{k} );
        
        % InP3R states increments
        phI43 = rI43.*(1 - hItp{k});
        phGate = mIp.*mItp{k}.*hItp{k};
        phI42 = rI42.*(    phGate);
        phI24 = rI24.*(1 - phGate);
        DI21 = rI21.*CI2p{k} - rI12.*CI1p{k};
        DI26 = rI26.*CI2p{k} - rI62.*OI6p{k};
        DI42 = phI42.*CI4p{k} - phI24.*CI2p{k};
        DI45 = rI45.*CI4p{k} - rI54.*OI5p{k};
        DI43 = phI43.*CI4p{k} - rI34.*CI3p{k};
        
        % probability of open state for InP3Rs
        OpenIp = OI5p{k} + OI6p{k};
        % Ca increment through open InP3R channels
            % process proper
        DCI(kpr) = krI*OpenIp(kpr).*(CeM(kpr)-CcM(kpr,1));
            % root segment
        DCI(krt) = krIr*OpenIp(krt).*(CeM(krt)-CcM(krt,1));

        % InP3 increment from generation at plasma membrane
        DInP3Mi = kmIf*CcM(:,3)./(1+kiPI*PKCvp{k});
        % InP3 decrement, bulk degradation
        DInP3Md = kmIb*InP3M;
        
        % PKCv update
        DPKCv = kiPf*CcM(:,3).*(PKC0-PKCvp{k}) - kiPb*PKCvp{k};
        
        % RyR states
        
        rR34 = kR34.*CcM(kpr,1).^3;
        rR13 = kR13.*CcM(kpr,1).^4;
        phR13 = rR13.*rsat13./(rR13+rsat13);
        phR34 = rR34.*rsat34./(rR34+rsat34);
        
        % RyR states increments
        DR13 = phR13.*CR1p{k} - rR31.*OR3p{k};
        DR32 = rR32.*OR3p{k} - rR23.*CR2p{k};
        DR34 = phR34.*OR3p{k} - rR43.*OR4p{k};
        
        % probability of open state for RyR
        OpenR = OR3p{k} + OR4p{k};
        
        % Ca increment through open RyR channels
        DCR(kpr) = krR*OpenR.*(CeM(kpr)-CcM(kpr,1));
        
        % Apply updates to all states %
        
        % Source/sink contributions to diffusive states (at grid centers)
        
        % sum of ER calcium increments
        DCeM = volRie*( DCS - DCI-DCR - DCle );
        DCeM(kpr) = DCeM(kpr) + DCer; % virtual source for process proper
        DCeM(krt) = DCeM(krt) + DCerb; % virtual source for process root               
        
        CeM = CeM+DCeM;
        
        % sum of cytosolic calcium increments
        DCcMb = DCbf-DCfb;
        CcM = CcM+DCcMb;
        DCcMe = -DCS + DCI+DCR + DCle;
        CcM(:,1) = CcM(:,1)+DCcMe;
        DCcMx = -DCN-DCP + DClx;
        CcM(:,3) = CcM(:,3)+DCcMx;
        
        % InP3 increments/decrements
        InP3M(:,3) = InP3M(:,3) + DInP3Mi;
        InP3M = InP3M - DInP3Md;
        
        % diffusion processes
        %   NOTE: proximal boundary conditions applied in final section
        
        % cytosolic calcium diffusion
        
        %   radial diffusion
        DCcR = DCac*(RD*CcM')';
        CcM = CcM + DCcR;
        
        % update Cc array
        Cc = [Cc(1,:); CcM; Cc(end,:)];
        %   axial diffusion
        DCcD = DCac*LP*Cc;
        Cc = Cc + DCcD;
        % apply distal sealed-end boundary condition
        Cc(1,:) = BL*Cc; % grad(Cc) (i.e. flux) = 0, distal end
        
       % ER calcium diffusion
        % update Ce array
        Ce = [Ce(1,1); CeM; Ce(end,1)];
        DCeD = DCae*LP*Ce;
        Ce = Ce + DCeD;
        % apply distal sealed-end boundary condition
        Ce(1,1) = BL*Ce; % grad(Ce) (i.e. flux) = 0, distal end
        
        % InP3 diffusion
        
        %   radial diffusion
        DInP3R = DInP3*(RD*InP3M')';
        InP3M = InP3M + DInP3R;
        
        % update InP3 array
        InP3 = [InP3(1,:); InP3M; InP3(end,:)];
        % axial diffusion
        DInP3D = DInP3*LP*InP3;
        InP3 = InP3 + DInP3D;
        % apply distal sealed-end boundary condition
        InP3(1,:) = BL*InP3; % grad(InP3) (i.e. flux) = 0, distal end

        % during input period, modify Ca & InP3 B.V.s to reflect in. fluxes
        if t < tin && ( k==1 )
            Cc(1,:) = Cc(1,:) + BVCain; % Ca pulse, distal end
            InP3(1,:) = InP3(1,:)+ BVInP3in;  % InP3 pulse, distal end
        end
        
        % replace Ca & Inp3 distributions in cell arrays
        Ccp{k} = Cc;
        Cep{k} = Ce;
        InP3p{k} = InP3;
        
        % Calbindin state
        
        CaCalBp{k} = CaCalBp{k} + DCfb - DCbf;
        
        % InP3R and related states
        CI1p{k} = CI1p{k} + DI21;
        CI2p{k} = CI2p{k} + DI42 - DI21 - DI26;
        CI3p{k} = CI3p{k} + DI43;
        CI4p{k} = CI4p{k} - DI45 - DI42 - DI43;
        OI5p{k} = OI5p{k} + DI45;
        OI6p{k} = OI6p{k} + DI26;
        
        PKCvp{k} = PKCvp{k} + DPKCv;
        
        % RyR states
        CR1p{k} = CR1p{k} - DR13;
        CR2p{k} = CR2p{k} + DR32;
        OR3p{k} = OR3p{k} + DR13 - DR32 - DR34;
        OR4p{k} = OR4p{k} + DR34;
        
    end
    
    %% CELL BODY EVOLUTION
    
    % Ca leakage from ER
    
    DCleb = kleb.*(Ceb-Ccb);

    % Ca buffering in cytosol (bulk)
    
    DCfbb = kbCf.*Ccb.*(CalB0-CaCalBb);
    DCbfb = kbCb.*CaCalBb;
 
    % Calcium pumps
    
    % SERCA pumps
    DCSb = kpSb.*Ccb./((kpSC+Ccb).*Ceb);

    % NCX pumps
    DCNb(ksw,nr,:) = 0; % initialize sidewall terms
    %   face terms
    for k=1:nr
        DCNb(:,k,1) = kpNbf(k).*Ccb(:,k,1)./(kpNC+Ccb(:,k,1));
        DCNb(:,k,end) = kpNbf(k).*Ccb(:,k,end)./(kpNC+Ccb(:,k,end));
    end
    %   increment sidewall terms
    DCNb(ksw,nr,:) = DCNb(ksw,nr,:) +...
                    kpNbw.*Ccb(ksw,nr,:)./(kpNC+Ccb(ksw,nr,:));
         
    % PMCA pumps
    Cc2pw = Ccb.*Ccb;
    DCPb(ksw,nr,:) = 0; % initialize sidewall terms
    %   face terms
    for k=1:nr
        DCPb(:,k,1) = kpPbf(k).*Cc2pw(:,k,1)./(kpPC2pw+Cc2pw(:,k,1));
        DCPb(:,k,end) = kpPbf(k).*Cc2pw(:,k,end)./(kpPC2pw+Cc2pw(:,k,end));
    end
    %   increment sidewall terms
    DCPb(ksw,nr,:) = DCPb(ksw,nr,:) +...
                    kpPbw.*Cc2pw(ksw,nr,:)./(kpPC2pw+Cc2pw(ksw,nr,:));
                     
    % InP3R and related states
        
    % Gating function dynamics
    Cc6pw = Ccb.^6;
    mI =  Cc6pw./(Cc6pw+km6pw);
    hI =  kh6pw./(Cc6pw+kh6pw);
    
    In3pw = InP3b.^3;
    mIp = In3pw./(In3pw+kbp3pw);
    
%     % No dynamics; instantaneous gating
%     mItb = mI;
%     hItb = hI;

    % single-pole delay
        % (done here since these play role in InPrR state updates)
    mItb = mItb + rIg.*( mI - mItb );
    hItb = hItb + rIg.*( hI - hItb );
    
    % InP3R states increments
    phI43 = rI43.*(1 - hItb);
    phGate = mIp.*mItb.*hItb;
    phI42 = rI42.*(    phGate);
    phI24 = rI24.*(1 - phGate);
    DI21b = rI21.*CI2b - rI12.*CI1b;
    DI26b = rI26.*CI2b - rI62.*OI6b;
    DI42b = phI42.*CI4b - phI24.*CI2b;
    DI45b = rI45.*CI4b - rI54.*OI5b;
    DI43b = phI43.*CI4b - rI34.*CI3b;

    % probability of open state for InP3Rs
    OpenIb = OI5b + OI6b;
    % Ca increment through open InP3R channels
    for k=1:nr
        DCIb(:,k,:) = krIb(k)*OpenIb(:,k,:).*(Ceb(:,k,:)-Ccb(:,k,:));
    end
    
    % InP3 increment from generation at plasma membrane
    DInP3Mib(ksw,nr,:) = 0;  % initialize sidewall terms
    %   face terms
    for k=1:nr
        DInP3Mib(:,k,1) = kmIfbf(k)*Ccb(:,k,1)./(1+kiPI*PKCvb(:,k,1));
        DInP3Mib(:,k,end) =kmIfbf(k)*Ccb(:,k,end)./(1+kiPI*PKCvb(:,k,end));
    end
    %   increment sidewall terms
    DInP3Mib(ksw,nr,:) = DInP3Mib(ksw,nr,:) +...
                       kmIfbw*Ccb(ksw,nr,:)./(1+kiPI*PKCvb(ksw,nr,:));         
    % InP3 decrement, bulk degradation
    DInP3Mdb = kmIb*InP3b;
        
    % PKCv update
    DPKCvb = kiPf*Ccb.*(PKC0-PKCvb) - kiPb*PKCvb;
    
    
    % Apply updates to all states
    
    % Source/sink contributions to diffusive states (at grid centers)

    % sum of ER calcium increments/decrements
    DCeb = ( DCSb - DCIb - DCleb )./volRb + DCerb;
    Ceb = Ceb+DCeb;
    
    % sum of cytosolic calcium increments/decrements
    DCcbv = DCbfb-DCfbb -DCSb + DCIb + DCleb -DCNb-DCPb;
    Ccb = Ccb + DCcbv;
 
    % PM leakage
    for k=1:nr
        Ccb(:,k,1) = Ccb(:,k,1) + DClxbf(k);
        Ccb(:,k,end) = Ccb(:,k,end) + DClxbf(k);
    end
    %   sidewall terms
    Ccb(:,nr,:) =   Ccb(:,nr,:) + DClxbw; % sidewall PM leakage fluxes

    % sum of InP3 increments/decrements
    InP3b = InP3b + DInP3Mib - DInP3Mdb;
    
    % Calbindin state
    CaCalBb = CaCalBb + DCfbb - DCbfb;
    
    % InP3R and related states
    CI1b = CI1b + DI21b;
    CI2b = CI2b + DI42b - DI21b - DI26b;
    CI3b = CI3b + DI43b;
    CI4b = CI4b - DI45b - DI42b - DI43b;
    OI5b = OI5b + DI45b;
    OI6b = OI6b + DI26b;
    
    PKCvb = PKCvb + DPKCvb;
    
    % diffusion processes
    
    % circumferential diffusion
    for k=1:nc % loop on circumferential index
        indx = mod( (k-kD):(k+kD-2), nc ) + 1; % contributing indices
        for kk=1:nr % loop on radial index
            LC = LB{kk}; % circumferential diffusion ceff's @ this radius
            DCcbC(k,kk,:) = DCac * pagemtimes(LC,Ccb(indx,kk,:)); % cyto Ca
            DCebC(k,kk,:) = DCae * pagemtimes(LC,Ceb(indx,kk,:)); % ER Ca
            DInP3bC(k,kk,:) = DInP3 * pagemtimes(LC,InP3b(indx,kk,:)); %InP3
        end
    end
    
    % radial diffusion
    %   temporary permuted arrays: put radial index first
    Ccbx = permute(Ccb,[2 1 3]);
    Cebx = permute(Ceb,[2 1 3]);
    InP3bx = permute(InP3b,[2 1 3]);
    %   non-junction elements: use RBW for sealed BC
    DCcbRx(:,ksw,:) = DCacr * pagemtimes(RBW,Ccbx(:,ksw,:)); % cyto Ca
    DCebRx(:,ksw,:) = DCaer * pagemtimes(RBW,Cebx(:,ksw,:)); % ER Ca
    DInP3bRx(:,ksw,:) = DInP3r * pagemtimes(RBW,InP3bx(:,ksw,:)); % InP3
    %   junction elements: use RBJ for joint BC
    DCcbRx(:,kjn,:) = DCacr * pagemtimes(RBJ,Ccbx(:,kjn,:)); % cyto Ca
    DCebRx(:,kjn,:) = DCaer * pagemtimes(RBJ,Cebx(:,kjn,:)); % ER Ca
    DInP3bRx(:,kjn,:) = DInP3r * pagemtimes(RBJ,InP3bx(:,kjn,:)); % InP3
    %   permute indices back to correct order
    DCcbR = permute(DCcbRx,[2 1 3]); % cyto Ca
    DCebR = permute(DCebRx,[2 1 3]); % ER Ca
    DInP3bR = permute(DInP3bRx,[2 1 3]); % InP3
    % contributions of boundary values for junction elements
    k0 = 0;
    for k=1:np % loop over junctions
        indx = k0 + (1:ncp);
        % cyto Ca
        DCcbR(indx,end,:) = DCcbR(indx,end,:) + DCacr*RBJb*BVCc(k);
        % ER Ca
        DCebR(indx,end,:) = DCebR(indx,end,:) + DCaer*RBJb*BVCe(k);
        % cyto InP3
        DInP3bR(indx,end,:) = DInP3bR(indx,end,:) + DInP3r*RBJb*BVInP3(k);
        k0 = k0 + nct;
    end

    % z-direction diffusion
    %   temporary permuted arrays: put z-index first
    Ccbx = permute(Ccb,[3 1 2]);
    Cebx = permute(Ceb,[3 1 2]);
    InP3bx = permute(InP3b,[3 1 2]);    
    for k=1:nr % loop on radial index
        TC = TB{k}; % z-direction diffusion ceff's @ this radius
        DCcbTx(:,:,k) = DCac * TC*Ccbx(:,:,k); % cyto Ca
        DCebTx(:,:,k) = DCae * TC*Cebx(:,:,k); % ER Ca
        DInP3bTx(:,:,k) = DInP3 * TC*InP3bx(:,:,k); % InP3
    end
    %   permute indices back to correct order
    DCcbT = permute(DCcbTx,[2 3 1]); % cyto Ca
    DCebT = permute(DCebTx,[2 3 1]); % ER Ca
    DInP3bT = permute(DInP3bTx,[2 3 1]); % InP3
    
    % update Ceb array
    Ceb = Ceb + DCebC + DCebR + DCebT;
    % update Cec array
    Ccb = Ccb + DCcbC + DCcbR + DCcbT;
    % update InP3b array
    InP3b = InP3b + DInP3bC + DInP3bR + DInP3bT;

    %% JOINT BOUNDARY CONDITIONS
    
    BVCc = zeros(1,np);
    BVCe = zeros(1,np);
    BVInP3 = zeros(1,np);
    k0 = 0;
    for k=1:np % loop over junctions
        
        % Contributions from process
        BVCc(k) = DCac.*fxp*(JP*Ccp{k})'; % cyto Ca
        BVCe(k) = DCae.*JP*Cep{k}; % ER elements
        BVInP3(k) = DInP3.*fxp*(JP*InP3p{k})';  % InP3
        
        % contributions from cell body elements
        %   compute mean values over thickness elements
        Ccbx = zeros(nc,nr);
        Cebx = zeros(nc,nr);
        InP3bx = zeros(nc,nr);
        for kk=1:nt
            Ccbx = Ccbx + Ccb(:,:,kk);
            Cebx = Cebx + Ceb(:,:,kk);
            InP3bx = InP3bx + InP3b(:,:,kk);
        end
        Ccbx = Ccbx./nt;
        Cebx = Cebx./nt;
        InP3bx = InP3bx./nt;
        %   compute cell body contributions
        indx = k0 + (1:ncp); % circumf. indices of junction elements
        BVCc(k) = BVCc(k) + DCacr.*sum(JB*Ccbx(indx,:)')./ncp;
        BVCe(k) = BVCe(k) + DCaer.*sum(JB*Cebx(indx,:)')./ncp;
        BVInP3(k) = BVInP3(k) + DInP3r.*sum(JB*InP3bx(indx,:)')./ncp;
        k0 = k0 + nct;
        
        % normalize
        BVCc(k) = denCc.*BVCc(k);
        BVCe(k) = denCe.*BVCe(k);
        BVInP3(k) = denInP3.*BVInP3(k);
        
        Ccp{k}(end,:) = BVCc(k);
        Cep{k}(end) = BVCe(k);
        InP3p{k}(end,:) = BVInP3(k);
        
    end

    %% Store plot data at intervals of tout
    if ~mod(t,tout)
        it = it+1;
        data_time = [data_time;t];
        % data for cell body
        Ccbt{it,1} = Ccb;
        Cebt{it,1} = Ceb;
        InP3bt{it,1} = InP3b;
        
        % data for processes
        for k=1:np
            Ccpt{it,k} = Ccp{k};
            Cept{it,k} = Cep{k};
            InP3pt{it,k} = InP3p{k};
            
        end
        
% graphing
    
%     % set plot figure for cytosolic Ca in process #1
% 
%         figure(1)
%         plot(xgrid1,Ccp{1}(:,1),'g.-');
%         hold on
%         grid on
%         plot(xgrid1,Ccp{1}(:,2),'.-','color',[0,0.75,0.75]);
%         plot(xgrid1,Ccp{1}(:,3),'b.-');
%         plot(xgrid1,InP3p{1}(:,1),'r-')
%         xlabel('Distance (um)');
%         ylabel('Concentration (uM)');
%         set(gca,'fontsize',15);
%         xlim(xlims);
%         xticks(xtic);
%         ylim([0 5]);
%         hold off

%         % plot other selected values
%         figure(2)
%         plot(rc,Ccb(2,:),'g.-');
%         hold on
%         grid on
%         plot(rc,Ccb(16,:),'.-','color',[0,0.75,0.75]);
%         plot(rc,Ccb(28,:),'b.-');
%         plot(rc,InP3b(4,:),'r-')
%         xlabel('Distance (um)');
%         ylabel('Concentration (uM)');
%         set(gca,'fontsize',15);
%         xlim([0,rc(end)+0.5*dr]);
%         ylim([0 9]);
%         hold off
        
    end
        
end

% filename for outputs
%outfile = ['OUT',mpfile(2:end)];
outfile = ['OUT','_',mpfile];
% save data for plotting
save(outfile, 'Ccbt','Ccpt','Cebt','Cept','InP3bt','InP3pt',...
               'ng','ngr','np','nc','nct','ncp','nr','nt','rb','rc','ld','data_time');
