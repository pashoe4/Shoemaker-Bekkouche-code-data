function out=ProcCaDynamics3D_ranRyR(mpfile, varargin)
%%%% ProcCaDynamics3D_ranRyR.m %%%%
% Compute equilibrium state of a compartment with ER + RYRs, followed by
% time-domain response of a cylindrical process with ER + RYRs,
% with 2-D (axial & radial) diffusion.
% 'Excitation' takes form of elevated Ca at distal end
% InP3 model is Siekmann-Cao-Sneyd model with re-defined CI3 state;
% RyR model is from Breit & Quiesser / Keizer & Levine.
% RyRs may be discretely distributed along dendrite, < 1/compartment,
% randomly distributed around the circumference of the ER
%
% User is prompted for a file containing desired parameter values.

dt = 1E-5;   % time step for temporal integration (units s)
%dt = 4E-6;   % small time step for temporal integration (units s)

iplot = 200; % plot interval in units of dt
tout = iplot*dt; % plot interval
it = 0;%index for saving time and data
data_time=[];%saving time for data processing after simulation
%% Load parameters

% load model parameters
load(mpfile); 

% load rate constants / function tables for InP3 receptors
load('InP3R_rates.mat');

% load rate constants / function tables for Ryanodine receptors
load('RyR_rates.mat');

%% Scale rate constants by dt for purposes of temporal integration

%%% Universal calcium parameters %%%

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

%% Set up initial conditions / initial states for equilibrium calculation

% Small calcium concentrations used to ensure convergence to equilibrium
% Moved to setParametersProc3D.m to enable setting as input with other vars
%Cc = 0.01; % cytosolic Ca; units uM
%Ce = 200; % ER Ca; units uM

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
%tend = 225; % end time
tend = 150; % end time
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
        % note scaling by 1/nc for circumferential divisions
    DCR = volRic*krR*OpenR.*( Ce - Cc )/(nc*MrR);
        
    % Apply increments to states
    
    % sum of ER calcium increments
    DCe = ( DCS - DCR-DCI - DCle )/volR;
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

doPlot=0;
if doPlot
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
end

Cc_Cmprt = Cc;
Ce_Cmprt = Ce;
InP3_Cmprt = InP3;

%% Set up initial conditions / initial states for dendrite

% load arrays for dendrite with approximate quiescent parameter values

Cc = Cc*ones(ng+2,3,nc); % cytosolic Ca
Ce = Ce*ones(ng+2,1); % ER Ca
InP3 = InP3*ones(ng+2,3,nc); % cytosolic InP3;
CaCalB = CaCalB*ones(ng,3,nc); % bound calcium value

% quiescent InP3R states
CI1 = CI1*ones(ng,1,nc);
CI2 = CI2*ones(ng,1,nc);
CI3 = CI3*ones(ng,1,nc);
CI4 = CI4*ones(ng,1,nc);
OI5 = OI5*ones(ng,1,nc);
OI6 = OI6*ones(ng,1,nc);
% gating variables
mIt = mIt*ones(ng,1,nc);
hIt = hIt*ones(ng,1,nc);

% activated PKC
PKCv = PKCv*ones(ng,1,nc);

% quiescent RyR  states
CR1 = CR1*ones(nRyR,1);
CR2 = CR2*ones(nRyR,1);
OR3 = OR3*ones(nRyR,1);
OR4 = OR4*ones(nRyR,1);

% pre-allocate arrays for 3D increments
DCcCx = zeros(nc,ng,3);
DInP3Cx = zeros(nc,ng,3);
DCcA = zeros(ng+2,3,nc);
DInP3A = zeros(ng+2,3,nc);

% allocate array for ER leak increments
DCle = zeros(ng,1,nc);

% allocate arrays for pump calcium increments
DCS = zeros(ng,1,nc);
DCP = zeros(ng,1,nc);
DCN = zeros(ng,1,nc);

% allocate arrays for InP3R state increments
DI21 = zeros(ng,1,nc);
DI26 = zeros(ng,1,nc);
DI42 = zeros(ng,1,nc);
DI45 = zeros(ng,1,nc);
DI43 = zeros(ng,1,nc);

% allocate arrays for RyR state increments
DR13 = zeros(nRyR,1);
DR32 = zeros(nRyR,1);
DR34 = zeros(nRyR,1);

% allocate arrays for receptor calcium increments
DCI = zeros(ng,1,nc);
DCR = zeros(ng,1,nc);

%% Set up for time loop

nin = round(4/dx); % number of compartments receiving input (4um length)

% impose initial condition trigger @ distal end to start response
    Ccin6pw = 2.0^6;
    mI =  Ccin6pw./(Ccin6pw+km6pw);
    hI =  kh6pw./(Ccin6pw+kh6pw);
for k = 1:nin+1 % set initial conditions in distal compartments
    Cc(k,:,:) = 2.0; % [Ca], units uM
    mIt(k,1,:) = mI;
    hIt(k,1,:) = hI;
    InP3(k,:,:) = 1.5; % [InP3], units uM
end

% grids for plots
xgrid1 = [0, 0.5*dx:dx:(ng-0.5)*dx, ng*dx]';
xgrid2 = [0.5*dx:dx:(ng-0.5)*dx]';
xlims = [-1,ld];
xtic = [0:ld/10:ld];

Ccmax = 0; % to keep track of spatial max Cc

% location markers for Ca wave computations
Nm = MrR;
m0 = ceil(ng/2);
k = Mr0;
while k<m0
    k = k+MrR;
end
m2 = k;
m1 = m2-Nm;
while k<ng
    k = k+MrR;
end
m3 = k-round(4*MrR0/MrR)*Nm;
IP = m1;

%% Loop on time

tend = 1.5; % end time

Cc_mid = []; % for minimal dendrite history

for t = 0:dt:tend % loop on time
    
    Cc(Cc<0) = 0; % added to prevent numerical undershoot in [Ca]
    InP3(InP3<0) = 0; % added to prevent numerical undershoot in [InP3]

%% internal Cc, Ce, and InP3R for use in computing states @ grid points
    CcM = Cc(2:end-1,:,:);
        % permute 2nd & 3rd (radial & circumferential) indices of CcM
    CcMx = permute(CcM,[1 3 2]);
    CcMy = CcMx(:,:,1); % 2D array for inner radial elements
        % average of CcMx around circumference
    CcMave=pagetranspose(mean(pagetranspose(CcMx)));
    
    CeM = Ce(2:end-1,1);
    
    InP3M = InP3(2:end-1,:,:);
    
%% graphing
    
    % set up figure for cytosolic Ca
     if ~mod(t,tout)

%         figure(1)
%         plot(xgrid1,Cc(:,1,1),'g.-');
%         hold on
%         grid on
%         plot(xgrid1,Cc(:,2,1),'.-','color',[0,0.75,0.75]);
%         plot(xgrid1,Cc(:,3,1),'b.-');
%         plot(xgrid1,InP3(:,1,1),'r-')
%         xlabel('Distance (um)');
%         ylabel('Concentration (uM)');
%         set(gca,'fontsize',15);
%         xlim(xlims);
%         xticks(xtic);
%         ylim([0 35]);
%         hold off
        
        %Store data
        it = it+1;
        data_time = [data_time;t];
        Ccpt{it,1} = Cc;

        if doPlot
            figure(1)
            plot(xgrid2,CcMave(:,1,1),'g.-');
            hold on
            grid on
            plot(xgrid2,CcMave(:,1,2),'.-','color',[0,0.75,0.75]);
            plot(xgrid2,CcMave(:,1,3),'b.-');
            plot(xgrid1,InP3(:,1,1),'r-')
            xlabel('Distance (um)');
            ylabel('Concentration (uM)');
            set(gca,'fontsize',15);
            xlim(xlims);
            xticks(xtic);
            ylim([0 20]);
            hold off
        
        
%         figure(1)
%         plot(xgrid1,Cc(:,1,1),'b.-');
%         hold on
%         grid on
%         plot(xgrid1,Cc(:,1,2),'.-','color',[0,0.6,0.75]);
%         plot(xgrid1,Cc(:,1,3),'g.-');
%         plot(xgrid1,InP3(:,1,1),'r-')
%         xlabel('Distance (um)');
%         ylabel('Concentration (uM)');
%         set(gca,'fontsize',15);
%         xlim(xlims);
%         xticks(xtic);
%         ylim([0 30]);
            hold off
    
            % set up figure for ER Ca
            figure(2)
            plot(xgrid1,Ce(:,1),'.-','color', [0.3,0.3,0.8]);
            hold on
            grid on
            xlabel('Distance (uM)');
            ylabel('Concentration (uM)');
            set(gca,'fontsize',15);
            xlim([0 ng*dx]);
            xticks(xtic);
            ylim([0 500]);
            hold off
        end

    end

%% Compute Ca wave propagation speed
    % find max cytosolic Ca at PM side (over all cricumf. elements)
    [M,Ic] = max(CcMave(:,1,3));
    [M,I]=max(M);
    I=Ic(I);
    if ( I>=m1 ) % wait until propagating wave has time to develop
        if M > Ccmax
            Ccmax = M; % keep running tab of maximum Ca
        end
        if ~mod(t,tout)
            if doPlot
                figure(1);
                hold on
                plot(xlims, [Ccmax,Ccmax]); % plot a line at peak value of wave
                hold off
            end
        end
    end
    if ( I==m2 && IP==m1 )
        t0 = t; % mark when wave peak first occurs at axial locatiom m2
        IP = m2;
    end
    if ( I==m3 && IP==m2 && M>0.1*Ccmax)
        vel = (m3-m2)*dx/(t-t0); % compute velocity when peak arrives @ m3
        IP = 0;
    end
    
%% Compute increments for all states %%

    %% Ca leakage from ER

    DCle = kle.*(CeM-CcM(:,1,:));

    %% Ca buffering in cytosol
    
    DCfb = kbCf.*CcM.*(CalB0-CaCalB);
    DCbf = kbCb.*CaCalB;
 
    %% Calcium pumps
    
    % NCX pumps
    DCN = kpN.*CcM(:,3,:)./(kpNC+CcM(:,3,:));
    
    % PMCA pumps
    Cc2pw = CcM(:,3,:).*CcM(:,3,:);
    DCP = kpP.*Cc2pw./(kpPC2pw+Cc2pw);

    % SERCA pumps
    DCS = kpS.*CcM(:,1,:)./((kpSC+CcM(:,1,:)).*CeM);

    %% InP3R and related states
        
    % Gating function dynamics
    Cc6pw = CcM(:,1,:).^6;
    mI =  Cc6pw./(Cc6pw+km6pw);
    hI =  kh6pw./(Cc6pw+kh6pw);
    
    In3pw = InP3M(:,1,:).^3;
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
    DCI = krI*OpenI.*(CeM-CcM(:,1,:));
    
    % InP3 increment from generation at plasma membrane
    DInP3Mi = kmIf*CcM(:,3,:)./(1+kiPI*PKCv);
    % InP3 decrement, bulk degradation
    DInP3Md = kmIb*InP3M;
    
    % PKCv update
    DPKCv = kiPf*CcM(:,3,:).*(PKC0-PKCv) - kiPb*PKCv;
    
    %% RyR states
    
    rR34 = kR34.*CcMy(KRyR).^3;
    rR13 = kR13.*CcMy(KRyR).^4;
    phR13 = rR13.*rsat13./(rR13+rsat13);
    phR34 = rR34.*rsat34./(rR34+rsat34);
    
    % RyR states increments
    DR13 = phR13.*CR1 - rR31.*OR3;
    DR32 = rR32.*OR3 - rR23.*CR2;
    DR34 = phR34.*OR3 - rR43.*OR4;

    % probability of open state for RyR
    OpenR = OR3 + OR4;
    % Ca increment through open RyR channels
    DCR = krR*OpenR.*(CeM(LRyR)-CcMy(KRyR));
        
%% Apply updates to all states %%

    %% Source/sink contributions to diffusive states (at grid centers)

    % sum of ER calcium increments
    DCeM = zeros(ng,1);
        % from RyRs
    DCeM(LRyR) = DCeM(LRyR) - DCR;
        % from other sources (summed over all circumf. elements)
    for k=1:nc
        DCeM = DCeM + ( DCS(:,1,k) - DCI(:,1,k) - DCle(:,1,k) );
    end
    DCeM = volRce * DCeM;
    CeM = CeM + DCeM;
    
    % sum of cytosolic calcium increments
        % from RyRs
    CcMy(KRyR) = CcMy(KRyR) + DCR;
    CcMx(:,:,1) = CcMy;
    CcM = permute(CcMx,[1 3 2]); % permute back to standard order
        % from other ER sources
    DCcMe = -DCS + DCI + DCle;
    CcM(:,1,:) = CcM(:,1,:) + DCcMe;
        % from PM sources
    DCcMx = -DCN-DCP + DClx;
    CcM(:,3,:) = CcM(:,3,:) + DCcMx;
        % from bulk reactions
    DCcMb = DCbf-DCfb;
    CcM = CcM+DCcMb;
    
    % InP3 increments/decrements
    InP3M(:,3,:) = InP3M(:,3,:) + DInP3Mi;
    InP3M = InP3M - DInP3Md;
    
    %% diffusion processes
    
    % ER calcium diffusion (1D only)
    
    %   update Ce array
    Ce = [Ce(1,1); CeM; Ce(end,1)];
    
    % axial diffusion
    DCeA = DCae*LP*Ce;
    
    % update array
    Ce = Ce + DCeA;
    % apply boundary conditions:
    Ce(1,1) = BL*Ce; % grad(Ce) (i.e. flux) = 0, distal end
    Ce(end,1) = BR*Ce; % grad(Ce) = 0, proximal end

    % Cytoplasmic calcium & InP3 diffusion (3D)
    
    %   radial diffusion
    %       temporary permuted arrays: put radial index first
    CcMx = permute(CcM,[2 1 3]);
    InP3Mx = permute(InP3M,[2 1 3]);
    %       perform multiplications
    DCcR = DCac * pagemtimes(RD,CcMx); % cyto Ca
    DInP3R = DInP3 * pagemtimes(RD,InP3Mx); % InP3
    %       permute indices back to correct order
    DCcR = permute(DCcR,[2 1 3]); % cyto Ca
    DInP3R = permute(DInP3R,[2 1 3]); % InP3

    % circumferential diffusion
    %       temporary permuted arrays: put circ. index first, radial last
    CcMx = permute(CcM,[3 1 2]);
    InP3Mx = permute(InP3M,[3 1 2]);
    %      perform multiplications
    for k=1:3
        DCcCx(:,:,k) = DCac * LC{k}*CcMx(:,:,k);
        DInP3Cx(:,:,k) = DInP3 * LC{k}*InP3Mx(:,:,k);
    end
    %   permute indices back to correct order
    DCcC = permute(DCcCx,[2 3 1]); % cyto Ca
    DInP3C = permute(DInP3Cx,[2 3 1]); % InP3

    % Update arrays
    CcM = CcM + DCcR + DCcC;
    InP3M = InP3M + DInP3R + DInP3C;
    % Expand arrays to include boundaries    
    for k=1:nc
        Cc(:,:,k) = [ Cc(1,:,k); CcM(:,:,k); Cc(end,:,k) ];
        InP3(:,:,k) = [ InP3(1,:,k); InP3M(:,:,k); InP3(end,:,k) ];
    end
    
    % axial diffusion
    %       perform multiplications
    for k=1:nc
        DCcA(:,:,k) = DCac * LP*Cc(:,:,k); % cyto Ca
        DInP3A(:,:,k) = DInP3 * LP*InP3(:,:,k); % InP3
    end
    
    % Update Cc and InP3 arrays
    Cc = Cc + DCcA;
    InP3 = InP3 + DInP3A;
    % apply boundary conditions:
    for k=1:nc
        Cc(1,:,k) = BL*Cc(:,:,k); % grad(Cc) (i.e. flux) = 0, distal end
        Cc(end,:,k) = BR*Cc(:,:,k); % grad(Cc) = 0, prox. end
        InP3(1,:,k) = BL*InP3(:,:,k); % grad(InP3) = 0, distal end
        InP3(end,:,k) = BR*InP3(:,:,k); % grad(InP3) = 0, prox. end
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
    
end

% filename for outputs
%outfile = ['OUT',mpfile(2:end)];
outfile = ['OUT','_',mpfile];
% save data for plotting
save(outfile, 'Cc','Ce','Cct','Cet','CcMave','InP3','xgrid1','xgrid2','Ccpt','data_time');
%plot(xgrid2,CcMave(:,1,3),'b.-');
%plot(xgrid1,InP3(:,1,1),'r-')

