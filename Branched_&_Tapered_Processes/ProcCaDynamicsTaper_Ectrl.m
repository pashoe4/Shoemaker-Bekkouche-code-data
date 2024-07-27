%%%% ProcCaDynamicsTaper.m %%%%
% Compute equilibrium state of a compartment with ER + RYRs, followed by
% time-domain response of a tapered cylindrical process with ER + InP3Rs,
% with 2-D (axial & radial) diffusion.
% 'Excitation' takes form of pulsed Ca & InP3 fluxes at distal end.
% InP3 model is Siekmann-Cao-Sneyd model with re-defined CI3 state;
% RyR model is currently not implemented.
%
% User is prompted for a file containing desired parameter values.
% Script uses CURRENT VERSION of input files, is intended for release.

clear;

%dt = 1E-5;   % time step for temporal integration (units s)
dt = 2E-6;   % small time step for temporal integration (units s)

iplot = 200; % plot interval in units of dt
tout = iplot*dt; % plot interval

%% Load parameters

% load model parameters
mpfile = input('Specify name of model parameter file: ','s');
load(mpfile); % tapered process model parameters

% load rate constants / function tables for InP3 receptors
load('InP3R_rates.mat');

% load rate constants / function tables for Ryanodine receptors
load('RyR_rates.mat');

%% Set virtual calcium source for ER
ker = -1.6;
%ker = -1.56E-16./vole * kvol;
%ker = -1.24E-16./vole * kvol;
%ker = 0;

%% Set distal Ca & InP3 influx parameters for activated processes

% influx rate
% Cain = 100; % units uM.s^-1
% InP3in = 100; % units uM.s^-1
Cain = 200; % units uM.s^-1
InP3in = 200; % units uM.s^-1

BVCain = BCc * Cain;
BVInP3in = BInP3 * InP3in;

% influx pulse duration
tin = 0.02; % units s

%% Scale rate constants by dt for purposes of temporal integration

%%% Universal calcium parameters %%%

% Virtual calcium source for ER
DCer = dt*ker; % (note virtual source assumed constant)

% Calcium molar leakage rates
DClx = dt*klx; % (note extracellular leakage assumed constant)
kle = dt*kle;

% Cytosol Ca diffusion coefficient
DCac = dt*DCac;
% ER Ca diffusion coefficient
DCae = dt*DCae;

% calcium buffering rates
kbCf = dt*kbCf;
kbCb = dt*kbCb;

%%% Calcium pump parameters %%%

% pump clearance rate constants
kpN = dt*kpN;
kpP = dt*kpP;
kpS = dt*kpS;

%%% InP3 receptor and related parameters %%%

% InP3 diffusion coefficient
DInP3 = dt*DInP3;

% Calcium molar influx rate for open InP3R channels
krI = dt*krI;

% InP3 generation rate
kmIf = dt*kmIf;
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

% Input flux increments
Cain = dt*Cain;
InP3in = dt*InP3in;

%% Take means of parameters for equilibrium calculation

Avflag = ~Drflag; % identify segments with MrR <= 1
DClxX = mean(DClx);
kleX = mean(kle);
kpNX = mean(kpN);
kpPX = mean(kpP);
kpSX = mean(kpS);
kmIfX = mean(kmIf);
krIX = mean(krI);
krRX = mean(krR);
MrRX = mean(MrR(Avflag));
DCerX = mean(DCer);

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

it = 0;

DClx0 = volRoc*DClxX; % rescaled external leakage for single compartment
toutr = 20*tout; % downsample outputs

tend = 140; % end time
for t = 0:dt:tend % loop on time
    
    % Compute state increments
    
    % Ca leakage from ER
    DCle = volRic*kleX.*(CeX-CcX);
    
    % Ca buffering in cytosol
    DCfb = kbCf.*CcX.*(CalB0-CaCalBX);
    DCbf = kbCb.*CaCalBX;
    
    % Calcium pumps
    
        % NCX pumps
    DCN = volRoc*kpNX.*CcX./(kpNC+ CcX);
    
    % PMCA pumps
    Cc2pw =  CcX.* CcX;
    DCP = volRoc*kpPX.*Cc2pw./(kpPC2pw+Cc2pw);
    
    % SERCA pumps
    DCS = volRic*kpSX.* CcX./((kpSC+ CcX).*CeX);
    
    % InP3R and related states
    
    % Gating function dynamics
    Cc6pw = CcX.^6;
    mI =  Cc6pw./(Cc6pw+km6pw);
    hI =  kh6pw./(Cc6pw+kh6pw);
    
    In3pw = InP3X.^3;
    mIp = In3pw./(In3pw+kbp3pw);
    
%     % No dynamics; instantaneous gating
%     mIt = mI;
%     hIt = hI;

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
    DCI = volRic*krIX*OpenI.*( CeX - CcX );
    
    % cytosolic InP3 increment
    DInP3c = (volRoc*kmIfX* CcX)./(1+kiPI*PKCvX) - kmIb*InP3X;
    
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
    
    DCR = volRic*krRX*OpenR.*( CeX - CcX )./MrRX;
        
    % Apply increments to states
    
    % sum of ER calcium increments
    DCe = ( DCS - DCR-DCI - DCle )/volR + DCerX;
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
    
    if ~mod(t,toutr)
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

% load arrays for dendrite with approximate quiescent parameter values

Cc = CcX*ones(ngt+2,3); % cytosolic Ca
Ce = CeX*ones(ngt+2,1); % ER Ca
InP3 = InP3X*ones(ngt+2,3); % cytosolic InP3;
CaCalB = CaCalBX*ones(ngt,3); % bound calcium value

%%%%%%%%
Cc(1:10,:) = 2; % kickstart wave in early segment!
InP3(1:10,:) = 1; % kickstart wave in early segment!
%%%%%%%%

% quiescent InP3R states
CI1 = CI1X*ones(ngt,1);
CI2 = CI2X*ones(ngt,1);
CI3 = CI3X*ones(ngt,1);
CI4 = CI4X*ones(ngt,1);
OI5 = OI5X*ones(ngt,1);
OI6 = OI6X*ones(ngt,1);
% gating variables
mIt = mItX*ones(ngt,1);
hIt = hItX*ones(ngt,1);

% activated PKC
PKCv = PKCvX*ones(ngt,1);

% quiescent RyR  states
CR1 = CR1X*IRyR;
CR2 = CR2X*IRyR;
OR3 = OR3X*IRyR;
OR4 = OR4X*IRyR;

%%%% debug arrays %%%%
irp = [11,31,39,45];
CR1X = [];
CR2X = [];
OR3X = [];
OR4X = [];
DR13X = [];
DR32X = [];
DR34X = [];
DR13P = [];
DR32P = [];
DR34P = [];
DR13N = [];
DR32N = [];
DR34N = [];

%%%%%%%%%%%%%%%%%%%%%%

%% Set up for time loop

% grids for plots
xgrid1 = [0, 0.5*dx:dx:(ngt-0.5)*dx, ngt*dx]';
xlims = [-1,ldt+1];
xtic = [0:2:ldt];
xtap = [ld0,ld+ld0];

Ccmax = 0; % to keep track of spatial max Cc

%% Loop on time

Cct = {};  % cell array for calcium
it = 0;
ngstart = ng0+1;  % index of CcM at start of taper
ngend = ng0+ng+2; % index of Cc at end of taper

tend = 1.6;
%tend = 1.2; % end time
for t = 0:dt:tend % loop on time
    
    Cc(Cc<0) = 0; % added to prevent numerical undershoot in [Ca]
    InP3(InP3<0) = 0; % added to prevent numerical undershoot in [InP3]

%% internal Cc, Ce, and InP3R for use in computing states @ grid points
    CcM = Cc(2:end-1,:);    
    CeM = Ce(2:end-1,1);
    InP3M = InP3(2:end-1,:);
    
%% Compute increments for all states %%

    %% Ca leakage from ER

    DCle = kle.*(CeM-CcM(:,1));

    %% Ca buffering in cytosol
    
    DCfb = kbCf.*CcM.*(CalB0-CaCalB);
    DCbf = kbCb.*CaCalB;
 
    %% Calcium pumps
    
    % NCX pumps
    DCN = kpN.*CcM(:,3)./(kpNC+CcM(:,3));
    
    % PMCA pumps
    Cc2pw = CcM(:,3).*CcM(:,3);
    DCP = kpP.*Cc2pw./(kpPC2pw+Cc2pw);

    % SERCA pumps
    DCS = kpS.*CcM(:,1)./((kpSC+CcM(:,1)).*CeM);

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
    mIt = mIt + rIg.*( mI - mIt );
    hIt = hIt + rIg.*( hI - hIt );
    
    % InP3R states increments
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
    DCI = krI.*OpenI.*(CeM-CcM(:,1));
    
    % InP3 increment from generation at plasma membrane
    DInP3Mi = kmIf.*CcM(:,3)./(1+kiPI*PKCv);
    % InP3 decrement, bulk degradation
    DInP3Md = kmIb.*InP3M;
    
    % PKCv update
    DPKCv = kiPf.*CcM(:,3).*(PKC0-PKCv) - kiPb*PKCv;
    
    %% RyR states
    
    rR34 = kR34.*CcM(:,1).^3;
    rR13 = kR13.*CcM(:,1).^4;
    phR13 = rR13.*rsat13./(rR13+rsat13);
    phR34 = rR34.*rsat34./(rR34+rsat34);
    
    % RyR states increments
    DR13 = phR13.*CR1 - rR31.*OR3;
    DR32 = rR32.*OR3 - rR23.*CR2;
    DR34 = phR34.*OR3 - rR43.*OR4;

    % probability of open state for RyR
    OpenR = OR3 + OR4;
    
    % Ca increment through open RyR channels
    DCR = krR.*OpenR.*(CeM-CcM(:,1));

%% Apply updates to all states %%

    %% Source/sink contributions to diffusive states (at grid centers)

    % sum of ER calcium increments 
    DCeM = volRie*( DCS - DCI-DCR - DCle ) + DCer;
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

    
    %% diffusion processes
    
    % ER calcium diffusion
    % update Ce array
    Ce = [Ce(1,1); CeM; Ce(end,1)];
    DCeD = DCae*LP*Ce;
    Ce = Ce + DCeD;
    % apply boundary conditions:
    Ce(1,1) = BL*Ce; % grad(Ce) (i.e. flux) = 0, distal end
    Ce(end,1) = BR*Ce; % grad(Ce) = 0, proximal end

    % InP3 diffusion
    
        %   radial diffusion
    DInP3R = DInP3*drm2pw.*(R0*InP3M')';
    InP3M = InP3M + DInP3R;
    
    % update InP3 array
    InP3 = [InP3(1,:); InP3M; InP3(end,:)];
    % axial diffusion
    DInP3D = DInP3*LP*InP3;
    InP3 = InP3 + DInP3D;
    % apply boundary conditions:
    InP3(1,:) = BL*InP3; % flux = 0, distal end
    InP3(end,:) = BR*InP3; % grad(InP3) = 0, prox. end

    % cytosolic calcium diffusion
    
    %   radial diffusion
    DCcR = DCac*drm2pw.*(R0*CcM')';
    CcM = CcM + DCcR;
    
    % update Cc array
    Cc = [Cc(1,:); CcM; Cc(end,:)];
    %   axial diffusion
    DCcD = DCac*LP*Cc;
    Cc = Cc + DCcD;
    % apply boundary conditions:
    Cc(1,:) = BL*Cc; % grad(Cc) (i.e. flux) = 0, distal end
    Cc(end,:) = BR*Cc; % grad(Cc) = 0, prox. end   

    % during input period, modify Ca & InP3 B.V.s to reflect input fluxes
    if t < tin
        Cc(1,:) = Cc(1,:) + BVCain; % Ca pulse, distal end
        InP3(1,:) = InP3(1,:) + BVInP3in;  % InP3 pulse, distal end
%         Cc(end,:) = Cc(end,:) + BVCain; % Ca pulse, proximal end
%         InP3(end,:) = InP3(end,:) + BVInP3in;  % InP3 pulse, proximal end
    end

    %% Calbindin state
    
    CaCalB = CaCalB + DCfb - DCbf;
    
    %% InP3R and related states
    CI1 = CI1 + DI21;
    CI2 = CI2 + DI42 - DI21 - DI26;
    CI3 = CI3 + DI43;
    CI4 = CI4 - DI45 - DI42 - DI43;
    OI5 = OI5 + DI45;
    OI6 = OI6 + DI26;
    
    PKCv = PKCv + DPKCv;
    
    %% RyR states
    CR1 = CR1 - DR13;
    CR2 = CR2 + DR32;
    OR3 = OR3 + DR13 - DR32 - DR34;
    OR4 = OR4 + DR34;
    
%% data archiving & graphing
    
    % set up figure for cytosolic Ca
    if ~mod(t,tout)
        
        % archive Cc data
        it = it+1;
        
        Cct{it} = Cc(ngstart:ngend,:);
        
%         CR1X = [ CR1X; CR1(irp)' ];
%         CR2X = [ CR2X; CR2(irp)' ];
%         OR3X = [ OR3X; OR3(irp)' ];
%         OR4X = [ OR4X; OR4(irp)' ];
%         DR13P = [DR13P; (phR13(irp).*CR1(irp))'];
%         DR32P = [DR32P;  (rR32.*OR3(irp))'];
%         DR34P = [DR34P; (phR34(irp).*OR3(irp))'];
%         DR13N = [DR13N; (-rR31.*OR3(irp))'];
%         DR32N = [DR32N; (-rR23.*CR2(irp))'];
%         DR34N = [DR34N; (-rR43.*OR4(irp))'];

        
        % graph Cc and InP3
        figure(3)
        plot(xgrid1,Cc(:,1),'g.-');
        hold on
        grid on
        plot(xgrid1,Cc(:,2),'.-','color',[0,0.75,0.75]);
        plot(xgrid1,Cc(:,3),'b.-');
        plot(xgrid1,InP3(:,1),'r-');
        plot(xtap,[0,0],'r^');
        xlabel('Distance (um)');
        ylabel('Concentration (uM)');
        set(gca,'fontsize',15);
        xlim(xlims);
        xticks(xtic);
        ylim([0 5]);
        hold off
        
%         pause
        
        % graph ER Ca
        figure(4)
        plot(xgrid1,Ce(:,1),'.-','color', [0.3,0.3,0.8]);
        hold on
        grid on
        plot(xtap,[0,0],'r^');
        xlabel('Distance (uM)');
        ylabel('Concentration (uM)');
        set(gca,'fontsize',15);
        xlim(xlims);
        xticks(xtic);
        ylim([0 450]);
        hold off
        
    end
    
end

rdt = rd(ngstart-1:ngend-1);
% filename for outputs
outfile = ['OUT',mpfile];
% save data for plotting
save(outfile, 'Cct','rdt','dx','tout');
