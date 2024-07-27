%%%% ProcCaDynamicsMultiSeg.m %%%%
% Compute equilibrium state of a compartment with ER + RYRs, followed by
% time-domain response of a multi-segment (possibly branching) process,
% with 2-D (axial & radial) diffusion.
% 'Excitation' takes form of pulsed InP3 flux at distal end.
% Initial values computed with average parameters
% over all segments w/o multiple RyRs per compartment.
% > SETS BOUNDARY VALUES BASED ON WEIGHTED B.C. VECTORS FROM MOLE <
% InP3 model is Siekmann-Cao-Sneyd model with re-defined CI3 state;
% RyR model is from Breit & Quiesser / Keizer & Levine.
% RyRs may be discretely distributed along dendrite, < 1/compartment
%
% User is prompted for a file containing desired parameter values.
% Script uses CURRENT VERSION of input files.
% Renamed from 'Cmprt_DendrCaDynamicsMultiSeg_v1c.m'

clear;

%% set up time parameters
%dt = 1E-5;   % time step for temporal integration (units s)
dt = 5E-6;   % small time step for temporal integration (units s)

iplot = 100; % plot interval in units of dt
tout = iplot*dt; % plot interval

%% Load parameters

% load model parameters
filename = input('Specify name of model parameter file: ','s');
load(filename);

% load rate constants / function tables for InP3 receptors
load('InP3R_rates.mat');

% load rate constants / function tables for Ryanodine receptors
load('RyR_rates.mat');

%% Set distal Ca & InP3 influx parameters for activation of process

% influx rate
Cain = 200; % units uM.s^-1
InP3in = 150; % units uM.s^-1

BVCain = BCc * Cain;
BVInP3in = BInP3 * InP3in;

% influx pulse duration
tin = 0.02; % units s

%% Scale rate constants by dt for purposes of temporal integration

%%% Universal calcium parameters %%%

% Calcium molar leakage rates
DClx = dt.*klx; % (note extracellular leakage assumed constant)
kle = dt.*kle;

% Cytosol Ca diffusion coefficient
DCac = dt*DCac;
% ER Ca diffusion coefficient
DCae = dt*DCae;

% calcium buffering rates
kbCf = dt*kbCf;
kbCb = dt*kbCb;

%%% Calcium pump parameters %%%

% pump clearance rate constants
kpN = dt.*kpN;
kpP = dt.*kpP;
kpS = dt.*kpS;

%%% InP3 receptor and related parameters %%%

% InP3 diffusion coefficient
DInP3 = dt*DInP3;

% Calcium molar influx rate for open InP3R channels
krI = dt.*krI;

% InP3 generation rate
kmIf = dt.*kmIf;
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
krR = dt.*krR;

% RyR receptor state rates
rR23 = dt*rR23;
rR32 = dt*rR32;
rR31 = dt*rR31;
rR43 = dt*rR43;

kR13 = dt*kR13;
kR34 = dt*kR34;

rsat13 = dt*rsat13;
rsat34 = dt*rsat34;

%% Take means of parameters for equilibrium calculation

Avflag = ~Drflag; % identify segments with MrR <= 1
volRocX = mean(volRoc(Avflag));
volRicX = mean(volRic(Avflag));
DClxX = mean(DClx(Avflag));
kleX = mean(kle(Avflag));
kpNX = mean(kpN(Avflag));
kpPX = mean(kpP(Avflag));
kpSX = mean(kpS(Avflag));
kmIfX = mean(kmIf(Avflag));
krIX = mean(krI(Avflag));
krRX = mean(krR(Avflag));
MrRX = mean(MrR(Avflag));

%% Set up initial conditions / initial states for equilibrium calculation

% Small calcium concentrations used to ensure convergence to equilibrium
CcX = 0.01; % cytosolic Ca; units uM
CeX = 200; % ER Ca; units uM

CaCalBX = CcX*kbCf*CalB0/kbCb; % bound calcium value

% approximate quiescent InP3-related states
InP3X = kmIfX*CcX/kmIb; % cytosolic InP3;
CI1X = 0;
CI2X = 0;
CI3X = 0;
CI4X = 1;
OI5X = 0;
OI6X = 0;
PKCvX = 0;
% gating variables
mItX = 0;
hItX = 1;

% approximate quiescent RyR states
CR1X = 1;
CR2X = 0;
OR3X = 0;
OR4X = 0;

%% Relaxation to observe (unexcited) equilibrium states, or instability

% output variables to examine history
Cct=[]; Cet = []; time = [];

DClx0 = volRocX*DClxX; % rescaled external leakage for single compartment
tend = 100; % end time
for t = 0:dt:tend % loop on time
    
    % Compute state increments
    
    % Ca leakage from ER
    DCle = volRicX*kleX.*(CeX-CcX);
    
    % Ca buffering in cytosol
    DCfb = kbCf.*CcX.*(CalB0-CaCalBX);
    DCbf = kbCb.*CaCalBX;
    
    % Calcium pumps
    
        % NCX pumps
    DCN = volRocX*kpNX.*CcX./(kpNC+ CcX);
    
    % PMCA pumps
    Cc2pw =  CcX.* CcX;
    DCP = volRocX*kpPX.*Cc2pw./(kpPC2pw+Cc2pw);
    
    % SERCA pumps
    DCS = volRicX*kpSX.* CcX./((kpSC+ CcX).*CeX);
    
    % InP3R and related states
    
    % Gating function dynamics
    Cc6pw = CcX.^6;
    mI =  Cc6pw./(Cc6pw+km6pw);
    hI =  kh6pw./(Cc6pw+kh6pw);
    
    In3pw = InP3X.^3;
    mIp = In3pw./(In3pw+kbp3pw);
    
%     % No dynamics; instantaneous gating
%     mItX = mI;
%     hItX = hI;

    % single-pole delay
        % (done here since these play role in InPrR state updates)
    mItX = mItX + rIg.*( mI - mItX );
    hItX = hItX + rIg.*( hI - hItX );

    % InP3 states increments
    phI43 = rI43.*(1 - hItX);
    phGate = mIp.*mItX.*hItX;
    phI42 = rI42.*(    phGate);
    phI24 = rI24.*(1 - phGate);
    DI21 = rI21.*CI2X - rI12.*CI1X;
    DI26 = rI26.*CI2X - rI62.*OI6X;
    DI42 = phI42.*CI4X - phI24.*CI2X;
    DI45 = rI45.*CI4X - rI54.*OI5X;
    DI43 = phI43.*CI4X - rI34.*CI3X;
    
    % probability of open state for InP3Rs
    OpenI = OI5X + OI6X;
    % Ca increment through open InP3R channels
    DCI = volRicX*krIX*OpenI.*( CeX - CcX );
    
    % cytosolic InP3 increment
    DInP3c = (volRocX*kmIfX* CcX)./(1+kiPI*PKCvX) - kmIb*InP3X;
    
    % PKCv update
    DPKCv = kiPf* CcX.*(PKC0-PKCvX) - kiPb*PKCvX;
    
    % RyR states
    
    rR13 = kR13.*CcX.^4;
    rR34 = kR34.*CcX.^3;
    phR13 = rR13.*rsat13./(rR13+rsat13);
    phR34 = rR34.*rsat34./(rR34+rsat34);

    % RyR states increments
    DR13 = phR13.*CR1X - rR31.*OR3X;
    DR32 = rR32.*OR3X - rR23.*CR2X;
    DR34 = phR34.*OR3X - rR43.*OR4X;
    
    % probability of open state for RyR
    OpenR = OR3X + OR4X;
    
    DCR = volRicX*krRX*OpenR.*( CeX - CcX )./MrRX;
        
    % Apply increments to states
    
    % sum of ER calcium increments
    DCe = ( DCS - DCR-DCI - DCle )/volR;
    % update Ce
    CeX = CeX + DCe;
    
    % sum of cytosolic calcium increment
    DCc = DCbf-DCfb - DCN-DCP-DCS + DCI+DCR + DClx0+DCle;
    % update Cc
    CcX = CcX + DCc;
    
    % InP3 concentration increment
    InP3X = InP3X+DInP3c;
    
    % Calbindin state
    CaCalBX = CaCalBX + DCfb - DCbf;
    
    % InP3R and related states
    CI1X = CI1X + DI21;
    CI2X = CI2X + DI42 - DI21 - DI26;
    CI3X = CI3X + DI43;
    CI4X = CI4X - DI45 - DI42 - DI43;
    OI5X = OI5X + DI45;
    OI6X = OI6X + DI26;
    
    PKCvX = PKCvX + DPKCv;
    
    % RyR states
    CR1X = CR1X - DR13;
    CR2X = CR2X + DR32;
    OR3X = OR3X + DR13 - DR32 - DR34;
    OR4X = OR4X + DR34;
    
    if ~mod(t,tout)
        time = [time;t];
        Cct = [Cct,CcX];
        Cet = [Cet,CeX];
    end
    
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

Cc_Cmprt = CcX;
Ce_Cmprt = CeX;
InP3_Cmprt = InP3X;

%% Set up initial conditions / initial states for dendrite

for k=1:ns % loop over segments
    
    % load cell arrays for dendrite with approximate quiescent parameter values
    Cc{k} = CcX*ones(ng(k)+2,3); % cytosolic Ca
    Ce{k} = CeX*ones(ng(k)+2,1); % ER Ca
    InP3{k} = InP3X*ones(ng(k)+2,3); % cytosolic InP3;
    CaCalB{k} = CaCalBX*ones(ng(k),3); % bound calcium value
    
    % quiescent InP3R states
    CI1{k} = CI1X*ones(ng(k),1);
    CI2{k} = CI2X*ones(ng(k),1);
    CI3{k} = CI3X*ones(ng(k),1);
    CI4{k} = CI4X*ones(ng(k),1);
    OI5{k} = OI5X*ones(ng(k),1);
    OI6{k} = OI6X*ones(ng(k),1);
    % gating variables
    mIt{k} = mItX*ones(ng(k),1);
    hIt{k} = hItX*ones(ng(k),1);
    
    % activated PKC
    PKCv{k} = PKCvX*ones(ng(k),1);
    
    % quiescent RyR  states
    if Drflag % if >= 1 RyR/compartment, initiate states @ each compartment
        CR1{k} = CR1X*ones(ng(k),1);
        CR2{k} = CR2X*ones(ng(k),1);
        OR3{k} = OR3X*ones(ng(k),1);
        OR4{k} = OR4X*ones(ng(k),1);
    else % if < 1 RyR/compartment, only initiate @ compartments where present
        CR1{k} = CR1X*IRyR{k};
        CR2{k} = CR2X*IRyR{k};
        OR3{k} = OR3X*IRyR{k};
        OR4{k} = OR4X*IRyR{k};
    end
    
end
    
% divisor for joint boundary value calculation at segment junctions
denBV = -1 / ( 2*JPb );

%% Set up for time loop

nin = round(4/dx); % number of compartments receiving input (4um length)

% define series of segments for plotting
kp = cnx(2,1);
xgrid{1} = [0, 0.5*dx:dx:(ng(kp)-0.5)*dx, ng(kp)*dx]';
x0 = xgrid{1}(end);
kplot = kp;
for k=1:nx
    kp = cnx(1,k);
    xg = [0, 0.5*dx:dx:(ng(kp)-0.5)*dx, ng(kp)*dx]';
    xgrid{k+1} = xg + x0;
    kplot = [kplot,kp];
    x0 = xgrid{k+1}(end);
end
nplot = length(kplot);
xlims = [0,xgrid{end}(end)];
xtic = [0:10*dx:xgrid{end}(end)];

%% Loop on time

tend = 0.6; % end time
is = 0; % counter for data storage

% output variables to examine history
Cct={}; Cet = {}; InP3t = {};

for t = 0:dt:tend % loop on time
    
    for k=1:ns
        
        Cc{k}(Cc{k}<0) = 0; % added to prevent numerical undershoot in [Ca]
        InP3{k}(InP3{k}<0) = 0; % added to prevent numerical undershoot in [InP3]
        
        %% internal Cc, Ce, and InP3R for use in computing states @ grid points
        CcM = Cc{k}(2:end-1,:);
        CeM = Ce{k}(2:end-1,1);
        InP3M = InP3{k}(2:end-1,:);
               
        %% Compute increments for all states %%
        
        %% Ca leakage from ER
        
        DCle = kle(k).*(CeM-CcM(:,1));
        
        %% Ca buffering in cytosol
        
        DCfb = kbCf.*CcM.*(CalB0-CaCalB{k});
        DCbf = kbCb.*CaCalB{k};
        
        %% Calcium pumps
        
        % NCX pumps
        DCN = kpN(k).*CcM(:,3)./(kpNC+CcM(:,3));
        
        % PMCA pumps
        Cc2pw = CcM(:,3).*CcM(:,3);
        DCP = kpP(k).*Cc2pw./(kpPC2pw+Cc2pw);
        
        % SERCA pumps
        DCS = kpS(k).*CcM(:,1)./((kpSC+CcM(:,1)).*CeM);
        
        %% InP3R and related states
        
        % Gating function dynamics
        Cc6pw = CcM(:,1).^6;
        mI =  Cc6pw./(Cc6pw+km6pw);
        hI =  kh6pw./(Cc6pw+kh6pw);
        
        In3pw = InP3M(:,1).^3;
        mIp = In3pw./(In3pw+kbp3pw);
        
        %     % No dynamics; instantaneous gating
        %     mIt = mI;
        %     hIt = hI;
        
        % single-pole delay
        % (done here since these play role in InPrR state updates)
        mIt{k} = mIt{k} + rIg.*( mI - mIt{k} );
        hIt{k} = hIt{k} + rIg.*( hI - hIt{k} );
        
        % InP3R states increments
        phI43 = rI43.*(1 - hIt{k});
        phGate = mIp.*mIt{k}.*hIt{k};
        phI42 = rI42.*(    phGate);
        phI24 = rI24.*(1 - phGate);
        DI21 = rI21.*CI2{k} - rI12.*CI1{k};
        DI26 = rI26.*CI2{k} - rI62.*OI6{k};
        DI42 = phI42.*CI4{k} - phI24.*CI2{k};
        DI45 = rI45.*CI4{k} - rI54.*OI5{k};
        DI43 = phI43.*CI4{k} - rI34.*CI3{k};
        
        % probability of open state for InP3Rs
        OpenI = OI5{k} + OI6{k};
        % Ca increment through open InP3R channels
        DCI = krI(k).*OpenI.*(CeM-CcM(:,1));
        
        % InP3 increment from generation at plasma membrane
        DInP3Mi = kmIf(k)*CcM(:,3)./(1+kiPI*PKCv{k});
        % InP3 decrement, bulk degradation
        DInP3Md = kmIb*InP3M;
        
        % PKCv update
        DPKCv = kiPf*CcM(:,3).*(PKC0-PKCv{k}) - kiPb*PKCv{k};
        
        %% RyR states
        
        rR34 = kR34.*CcM(:,1).^3;
        rR13 = kR13.*CcM(:,1).^4;
        phR13 = rR13.*rsat13./(rR13+rsat13);
        phR34 = rR34.*rsat34./(rR34+rsat34);
        
        % RyR states increments
        DR13 = phR13.*CR1{k} - rR31.*OR3{k};
        DR32 = rR32.*OR3{k} - rR23.*CR2{k};
        DR34 = phR34.*OR3{k} - rR43.*OR4{k};
        
        % probability of open state for RyR
        OpenR = OR3{k} + OR4{k};
        
        % Ca increment through open RyR channels
        DCR = krR(k).*OpenR.*(CeM-CcM(:,1));
        
        %% Apply updates to all states %%
        
        %% Source/sink contributions to diffusive states (at grid centers)
        
        % sum of ER calcium increments
        DCeM = volRie(k)*( DCS - DCI-DCR - DCle );
        CeM = CeM+DCeM;
        
        % sum of cytosolic calcium increments
        DCcMb = DCbf-DCfb;
        CcM = CcM+DCcMb;
        DCcMe = -DCS + DCI+DCR + DCle;
        CcM(:,1) = CcM(:,1)+DCcMe;
        DCcMx = -DCN-DCP + DClx(k);
        CcM(:,3) = CcM(:,3)+DCcMx;
        
        % InP3 increments/decrements
        InP3M(:,3) = InP3M(:,3) + DInP3Mi;
        InP3M = InP3M - DInP3Md;
        
        %% Calbindin state
        
        CaCalB{k} = CaCalB{k} + DCfb - DCbf;
        
        %% InP3R and related states
        CI1{k} = CI1{k} + DI21;
        CI2{k} = CI2{k} + DI42 - DI21 - DI26;
        CI3{k} = CI3{k} + DI43;
        CI4{k} = CI4{k} - DI45 - DI42 - DI43;
        OI5{k} = OI5{k} + DI45;
        OI6{k} = OI6{k} + DI26;
        
        PKCv{k} = PKCv{k} + DPKCv;
        
        %% RyR states
        CR1{k} = CR1{k} - DR13;
        CR2{k} = CR2{k} + DR32;
        OR3{k} = OR3{k} + DR13 - DR32 - DR34;
        OR4{k} = OR4{k} + DR34;
        
        %% diffusion processes
        
        % ER calcium diffusion
        % update Ce array
        Ce{k} = [Ce{k}(1,1); CeM; Ce{k}(end,1)];
        DCeD = DCae*LP{k}*Ce{k};
        Ce{k} = Ce{k} + DCeD;
        
        % InP3 diffusion
        
        %   radial diffusion
        DInP3R = DInP3*(RD{k}*InP3M')';
        InP3M = InP3M + DInP3R;
        
        % update InP3 array
        InP3{k} = [InP3{k}(1,:); InP3M; InP3{k}(end,:)];
        % axial diffusion
        DInP3D = DInP3*LP{k}*InP3{k};
        InP3{k} = InP3{k} + DInP3D;
        
        % cytosolic calcium diffusion
        
        %   radial diffusion
        DCcR = DCac*(RD{k}*CcM')';
        CcM = CcM + DCcR;
        
        % update Cc array
        Cc{k} = [Cc{k}(1,:); CcM; Cc{k}(end,:)];
        %   axial diffusion
        DCcD = DCac*LP{k}*Cc{k};
        Cc{k} = Cc{k} + DCcD;
        
        % apply sealed-end boundary conditions
        
        if sbcl(k)
            Ce{k}(1,1) = BL{k}*Ce{k}; % grad(Ce) (flux) = 0, distal end
            InP3{k}(1,:) = BL{k}*InP3{k}; % grad(InP3) = 0, distal end
            Cc{k}(1,:) = BL{k}*Cc{k}; % grad(Cc) = 0, distal end
        end
        
        if sbcr(k)
            Ce{k}(end,1) = BR{k}*Ce{k}; % grad(Ce) = 0, proximal end
            InP3{k}(end,:) = BR{k}*InP3{k}; % grad(InP3) = 0, prox. end
            Cc{k}(end,:) = BR{k}*Cc{k}; % grad(Cc) = 0, prox. end
        end
        
    end
    
    % during input pulse, modify Ca & InP3 B.V.s @ distal end of segment 1
    if t < tin
        Cc{1}(1,:) = Cc{1}(1,:) + BVCain; % Ca pulse, distal end
        InP3{1}(1,:) = InP3{1}(1,:) + BVInP3in;  % InP3 pulse, distal end
        Cc{2}(1,:) = Cc{2}(1,:) + BVCain; % Ca pulse, distal end
        InP3{2}(1,:) = InP3{2}(1,:) + BVInP3in;  % InP3 pulse, distal end
    end

    %% axial boundary values for diffusants
    
    for k=1:nx % loop over junctions
        
        % Compute nodal boundary values for Cc, Ce, and InP3
        k2 = cnx(2,k); % index of first distal segment
        CeB = JPR{k2}*Ce{k2};
        InB = JPR{k2}*InP3{k2};
        CcB = JPR{k2}*Cc{k2};
        k3 = cnx(3,k);
        if k3 % if 2nd distal  segment exists
            CeB = 0.5*( CeB + JPR{k3}*Ce{k3} );
            InB = 0.5*( InB + JPR{k3}*InP3{k3} );
            CcB = 0.5*( CcB + JPR{k3}*Cc{k3} );
        end
        k1 = cnx(1,k); % index of proximal segment
        CeB = CeB + JPL{k1}*Ce{k1};
        InB = InB + JPL{k1}*InP3{k1};
        CcB = CcB + JPL{k1}*Cc{k1};

        CcB = denBV .* CcB; % boundary values of Cc for this node
        CeB = denBV .* CeB; % boundary value  of Ce for this node
        InB = denBV .* InB; % boundary values of InP3 for this node
        
        % set the boundary values in the arrays
        Ce{k1}(1) = CeB;
        InP3{k1}(1,:) = InB;
        Cc{k1}(1,:) = CcB;
        Ce{k2}(end) = CeB;
        InP3{k2}(end,:) = InB;
        Cc{k2}(end,:) = CcB;
        if k3 % if 2nd distal  segment exists
            Ce{k3}(end) = CeB;
            InP3{k3}(end,:) = InB;
            Cc{k3}(end,:) = CcB;
        end
        
    end
    
        % set up figure for cytosolic Ca
         if ~mod(t,tout)
    
             for np = 1:nplot
                kp = kplot(np);
                
                figure(1)
                plot(xgrid{np},Cc{kp}(:,1),'g.-');
                xlim(xlims);
                xticks(xtic);
                ylim([0 8]);
                hold on
                grid on
                plot(xgrid{np},Cc{kp}(:,2),'.-','color',[0,0.75,0.75]);
                plot(xgrid{np},Cc{kp}(:,3),'b.-');
                plot(xgrid{np},InP3{kp}(:,1),'r-')
                xlabel('Distance (um)');
                ylabel('Concentration (uM)');
                set(gca,'fontsize',15);
               
%                 % set up figure for ER Ca
%                 figure(2)
%                 plot(xgrid{np},Ce{kp}(:,1),'.-','color', [0.3,0.3,0.8]);
%                 xlim(xlims);
%                 xticks(xtic);
%                 ylim([0 500]);
%                 hold on
%                 grid on
%                 xlabel('Distance (uM)');
%                 ylabel('Concentration (uM)');
%                 set(gca,'fontsize',15);
             end
             
             figure(1)
             hold off
             
%              figure(2)
%              hold off
             
             % store data at plot intervals
             is = is+1;
             
             % data for segment 1
             kp = kplot(1);
             CcX = Cc{kp};
             InP3X = InP3{kp};
             CeX = Ce{kp};
             % append data for subsequent segments to be plotted
             for np = 2:nplot
                 kp = kplot(np);
                 CcX = [ CcX; Cc{kp}(2:end,:) ];
                 InP3X = [ InP3X; InP3{kp}(2:end,:) ];
                 CeX = [ CeX; Ce{kp}(2:end,:) ];
             end
             
             % store data
             Cct{is} = CcX;
             InP3t{is} = InP3X;
             Cet{is} = CeX;
             
             
         end
        
end

% filename for outputs
filename = ['OUT',filename(2:end)];
% save data for plotting
save(filename, 'xgrid', 'tout', 'Cct', 'InP3t', 'Cet');

