function [MModeBFRF, MModeBFEnv, MModeBFdB] = PostAvgBMFMZeta(PESingle,bScanParm,XPt,YPt,GraphIt,DataToLoad)
%% PostAvgBMFM.m
%% Need to change the plot ranges (x,y and color)
% Need to change the plotting with higher channels
% Need to add angular plotting option
% Need to check delay parameter/distance in Z
% Change Noise Selection ROI

% Arguments -
%   PESingle (4-6D Array of Doubles) -  (FastTime x Slow Time (i.e., Pulse Number) x Element x Average (if DataToLoad is UnAvg or UnShift))
%   bScanParm (Struc) - Structure of all scan parameters
%   XPt (double) - X point in Scan
%   YPt (double) - Y point in Scan
%   GraphIt (double) - 0, 1, or 2; 2 is graph B-Mode; 1 is graph M-Mode; 0 is do not graph
%   DataToLoad (String) - 'UnShift','UnAvg',or 'Avg' indicating which type of data is being loaded
% Outputs -
%   MModeBFRF (3-4D Array of Doubles) - [Fast Time, PE Slow Time (i.e., Pulse No), Average (if UnAvg or UnShift)]   
%   MModeBFEnv (3-4D Array of Doubles) - [Fast Time, PE Slow Time (i.e., Pulse No), Average (if UnAvg or UnShift)]   
%   MModeBFdB (3-4D Array of Doubles) - [Fast Time, PE Slow Time (i.e., Pulse No), Average (if UnAvg or UnShift)]   

% Written by Alex Alvarez and modified from MModeAnal2D %190709

BMFMTime = tic;
cd ..\PEData % change to the appropriate path location for PE data
DisStruc = load([bScanParm.Format.filename2 '_PEParm.mat']); % Load the path
bScanParm = DisStruc.bScanParm; % [struc] Defines the parameters used for the scans
Trans = DisStruc.Trans; % [struc] Defines the transducer structure used in scan
Receive = DisStruc.Rcv; % [struc] Defines the receive structure used in scan
TX = DisStruc.TX; % [struc] Defines the transmit structure used in scan
TW = DisStruc.TW; % [struc] Defines the transmit waveform structure used in scan
twoD = bScanParm.twoD;
MModeType = 'Angular';
lat_mm = linspace(-bScanParm.Scan.XDist/2+bScanParm.LatFocus,bScanParm.Scan.XDist/2+bScanParm.LatFocus,bScanParm.Scan.Xpt);
ele_mm = linspace(-bScanParm.Scan.YDist/2+bScanParm.ElFocus,bScanParm.Scan.YDist/2+bScanParm.ElFocus,bScanParm.Scan.Ypt);
% Theta = atan(lat_mm/bScanParm.AxialFocus);
% Theta = atan(lat_mm(XPt)/bScanParm.AxialFocus);
% if Theta <= 0
%     Theta = Theta + pi/2;
% else
%     Theta = pi/2-Theta;
% end
% Phi = atan(ele_mm(YPt)/bScanParm.AxialFocus);
% if Phi <= 0
%     Phi = Phi + pi/2;
% else
%     Phi = pi/2-Phi;
% end
[Theta,Phi,R] = cart2sph(lat_mm(XPt),bScanParm.AxialFocus,ele_mm(YPt));
if Trans.name == 'H235'
    if bScanParm.twoAcq
        NumAz = 18;
        NumEl = 7;

    else
        NumAz = 12;
        NumEl = 5;
    end
    Trans.spacing_el = (Trans.elevationApertureMm/bScanParm.wavlen)/NumEl;
elseif Trans.name == 'H247'
    if bScanParm.twoAcq
        NumAz = 7;
        NumEl = 18;
    else
        NumAz = 5;
        NumEl = 12;
    end
end
if twoD
    Chantaje_X = floor(NumAz/2);
    Chantaje_Y = floor(NumEl/2);
else
    Chantaje_X = floor(Trans.numelements/2);
end

c              = 1.54;   % speed of sound mm/us SAMPLE (user defined)

t_offset       = 0;%100; % offset time to start loading data in sec
InterpValue     = 41; %numel(BF.z);

% Plotting parameters
dBRng          = [-10,0];
useREF         = 0; 
refVal         = 70;%59.2175; % empirical value in dB to define absolute max for RF data plots (only applies if useRef=1)
linRngFac      = 0.9;
linedBRng      = [-10 0];

% EUNIL Beamforming parameters
delays    = 2*max(TX.Delay(:));%15.5;%10.5;%12:0.25:19; %in usec ....% optimal delay is ~16 usec for this data

% Define ROI for beamforming (otherwise you entire field)
doROI = 0;         % 0 beamforms to all possible depth and lateral points 
x0xf = [-bScanParm.Scan.XDist/2 bScanParm.Scan.XDist/2];   % mm
if twoD
    y0yf = [-bScanParm.Scan.YDist/2 bScanParm.Scan.YDist/2]; % mm
end
z0zf = [bScanParm.AxialFocus-bScanParm.PlusMinusDist bScanParm.AxialFocus+bScanParm.PlusMinusDist];    % mm

% Define parameters for beamforming
BF.dx   = 0.2; % in mm [X pixel size for beamforming]
BF.dz   = 0.2; % in mm [Z pixel size for beamforming]

if twoD
    BF.dx   = 0.2; % in mm 
    BF.dy   = 0.2; % in mm [Y pixel size for beamforming]
    BF.dz   = 0.2; % in mm
end
BF.fnum = 1.5; % f_num of synethic aperture (typically 1.5 for an imaging system)
BF.fs   = 100;  % MHz for beamforming (upsample fo beamforming--typically 10+ X Nyquist sampling)

BF.depthFiltOn   = 0;          % Depth Filter Flag
BF.latFiltSmthOn = 0;          % Lateral Filter Flag

BF.filtCutMHZ    = [0.2 0.3 1 2];  % MHZ
BF.gLat          = 7;      % lateral smoothing filter width after beamforming in pixels.

fsV     = Receive(1).samplesPerWave*Receive(1).demodFrequency; % Sample rate in MHz (Verasonics)
cV      = c; %Resource.Parameters.speedOfSound/1000;               % speed of sound mm/us VERASONICS (fixed)

if twoD
    PlotXPt = 0;
    PlotYPt = 0;
    PlotZPt = bScanParm.AxialFocus;
end
if strcmpi(DataToLoad,'UnShift') || strcmpi(DataToLoad,'UnAvg')  %%size(PESingle,4) ~= bScanParm.Scan.Xpt && size(PESingle,4) ~= bScanParm.Scan.YPt
    LoopTo = size(PESingle,4);
else
    LoopTo = 1;
end
for u = 1:LoopTo % u is the number of averages
    for x = 1:size(PESingle,2) % x is the number of pulses
        % Get Verasonics RF data parameters
        if strcmpi(DataToLoad,'UnShift') || strcmpi(DataToLoad,'UnAvg')
            RFdata   = squeeze(PESingle(:,x,:,u)); RFMAX = max(RFdata(:));
        else
            RFdata   = squeeze(PESingle(:,x,:)); RFMAX = max(RFdata(:));
        end
        % Apply depth filter if flag is on
        if BF.depthFiltOn==1
            RFdata = 2*real(fft_filt2D_new(RFdata,[],fsV,BF.filtCutMHZ,0));
        end
        if ~twoD
            RcvChNx   = size(RFdata,2);               % Total # of receive channels [X]
        else
            RcvChNx   = NumAz;
            RcvChNy   = NumEl;   
            RFdata = reshape(RFdata,size(RFdata,1),RcvChNx,RcvChNy);
        end
        % Get Verasonics RF data time coordinates
        RFt0     = bScanParm.wavlen*Receive(1).startDepth/cV;%+Receive(1).startSample/fsV; % [usec] RF data start time
        RFtf     = bScanParm.wavlen*Receive(1).endDepth/cV;%+Receive(1).endSample/fsV;   % [usec] RF data end time
        RFtL     = RFtf-RFt0;                  % [usec] RF data total time
        tUS      = linspace(RFt0,RFtf,Receive(1).endSample - Receive(1).startSample+1); % [usec] RF data all time pts

        SlowTimet0 = 0;
        SlowTimetf = bScanParm.Scan.Duration_ms;
        SlowTime = linspace(SlowTimet0,SlowTimetf,size(PESingle,2));
        
        upSampleRatio = BF.fs/fsV;
        Ntime_BF      = round(length(tUS).*upSampleRatio);          % [pts] Total time points after upsampling
        t_BF          = linspace(tUS(1),tUS(end),Ntime_BF)';    % [usec] All time points after upsampling

        % RF data depth coordinates
        RFz0MM   = RFt0*c;                     % [mm] RF data start depth Z
        RFzfMM   = RFtf*c;                     % [mm] RF data end depth Z
        RFzL     = RFzfMM-RFz0MM;              % [mm] RF data total depth Z
        RFzMM    = tUS*c;                      % [mm] RF data all depth Z

        % RF data lateral coordinates
        RFx0   = 0;                                      % [mm] RF data start lateral X
        RFxf   = bScanParm.wavlen*Trans.spacing*(RcvChNx);    % [mm] RF data end lateral X
        midPt  = (RFxf - RFx0)/2;                        % [mm] RF data mid point lateral X
        RFxL   = RFxf;                                   % [mm] RF data mid point lateral X
        RFxMM  = linspace(RFx0-midPt,RFxf-midPt,RcvChNx); % [mm] RF data all pts lateral X (assumes center spacing around array)
        
        
        % RF data elevational coordinates
        if twoD
            RFy0   = 0;                                      % [mm] RF data start lateral X
            RFyf   = bScanParm.wavlen*Trans.spacing_el*(RcvChNy);    % [mm] RF data end lateral X
            midPt  = (RFyf - RFy0)/2;                        % [mm] RF data mid point lateral X
            RFyL   = RFyf;                                   % [mm] RF data mid point lateral X
            RFyMM  = linspace(RFy0-midPt,RFyf-midPt,RcvChNy); % [mm] RF data all pts lateral X (assumes center spacing around array)
            bScanParm.p.yROI = RFyMM;
        end
        bScanParm.p.xROI = RFxMM;
        bScanParm.p.zROI = RFzMM;
        bScanParm.p.tROI = tUS;

        % Do EUNIL beamforming for each delay 'delays' specified above
        delayIdx       = 1;
        for delayLoop=delays
            BF.offsetDelay = delayLoop;             % offset term in usec [not clear why we have so much offset--need to investigate]
            minOffset      = max(1,ceil(BF.offsetDelay*BF.fs/upSampleRatio+1)); % in case there is a delay from zero--needed for positive delays
            BF.x = bScanParm.p.xROI;                 % [mm] lateral ROI to beamform to
            BF.z = bScanParm.p.zROI(minOffset:end);  % [mm] Depth ROI to beamform to----Because a offset shift required, don't try to beamform shallow (otherwise would need zero padding)
%             BF.x = BF.x(1):BF.dx:BF.x(end);
%             BF.z = BF.z(1):BF.dz:BF.z(end);
            BF.x = linspace(BF.x(1),BF.x(end),InterpValue);
            BF.dx = BF.x(end)-BF.x(1)/InterpValue;
            BF.z = linspace(BF.z(1),BF.z(end),InterpValue);
            BF.dz = BF.z(end)-BF.z(1)/InterpValue;

            if twoD
                BF.y = bScanParm.p.yROI;
%                 BF.y = BF.y(1):BF.dy:BF.y(end);
                BF.y = linspace(BF.y(1),BF.y(end),InterpValue);
                BF.dy = BF.y(end)-BF.y(1)/InterpValue;
            end

            if doROI==1
                BF.x = BF.x(find(BF.x>=x0xf(1) & BF.x<=x0xf(2)));
                if twoD
                    BF.y = BF.y(find(BF.y>=y0yf(1) & BF.y<=y0yf(2)));
                end
                BF.z = BF.z(find(BF.z>=z0zf(1) & BF.z<=z0zf(2)));
            else
                x0xf = [BF.x(1) BF.x(end)];
                if twoD
                    y0yf = [BF.y(1) BF.y(end)];
                end
                z0zf = [BF.z(1) BF.z(end)];                     
            end
            if ~twoD
                BF.X = repmat(BF.x,[length(BF.z) 1]);
                BF.Z = repmat(BF.z',[1 length(BF.x)]);
            else
                BF.X = repmat(BF.x,[length(BF.z) 1 length(BF.y)]);
                BF.Y = repmat(BF.y,[length(BF.z) 1 length(BF.x)]);
                BF.Y = permute(BF.Y, [1 3 2]);
                BF.Z = repmat(BF.z',[1 length(BF.x) length(BF.y)]);                        
            end

            BF.os  = BF.offsetDelay.*BF.fs;	% [index] Offset into wavefield
            BF.osZ = BF.os/BF.fs*c;%BF.os/BF.fs*c/2;  % [mm] Depth offset into wavefield

            BF.aperEx = 1:RcvChNx;            % All elements for beamforming (normally full aperture)
            BF.Nax    = length(BF.aperEx);    % Total # of elements for beamforming
            BF.xa    = RFxMM(BF.aperEx);	    % [mm] total width of aperture for beamforming
            if twoD
                BF.aperEy = 1:RcvChNy;
                BF.Nay = length(BF.aperEy);
                BF.ya = RFyMM(BF.aperEy);
            end
            % Precalculate relevant beamforming parameters
            if ~twoD
                TI0    = zeros(length(BF.z),length(BF.x),BF.Nax); % Indexes into the wavefield matrix for each beamformed location and for each element
            else
                TI0    = zeros(length(BF.z),length(BF.x),length(BF.y),BF.Nax,BF.Nay);
            end
            DT     = TI0; % Interpixel differences            for each beamformed location and for each element
            wtemp  = TI0; % Element Weighting due to F-number for each beamformed location and for each element

            for idx = BF.aperEx
                if ~twoD
                    TI             = (sqrt((BF.X - BF.xa(idx)).^2 + (BF.Z).^2))/c*BF.fs - BF.os;%(sqrt((BF.X - BF.xa(idx)).^2 + (BF.Z).^2))/c*BF.fs - BF.os;%  % calculate transmission distance
                    TI0(:,:,idx)   = floor(TI);                                           % integer distance floored
                    DT(:,:,idx)    = TI - TI0(:,:,idx);                                   % remainder
                    wtemp(:,:,idx) = ((abs(BF.X - BF.xa(idx))./(BF.Z)) < (1/BF.fnum/2)); % make weights for element  choice
                else
                    for idy = BF.aperEy
                        TI         = (sqrt((BF.X - BF.xa(idx)).^2 + (BF.Y-BF.ya(idy)).^2 + (BF.Z).^2))/c*BF.fs-BF.os;
                        TI0(:,:,:,idx,idy)   = floor(TI);                                           % integer distance floored
                        DT(:,:,:,idx,idy)    = TI - TI0(:,:,:,idx,idy);                                   % remainder
                        wtempx = ((abs(BF.X - BF.xa(idx))./(BF.Z)) < (1/BF.fnum/2)); % make weights for element  choice
                        wtempy = ((abs(BF.Y - BF.ya(idy))./(BF.Z)) < (1/BF.fnum/2)); % make weights for element  choice
                        wtemp(:,:,:,idx,idy)  = wtempx.*wtempy;                               
                    end
                end
            end
            rhs = (1 + DT).*wtemp; %right hand side used lower down for beamforming

            badIDX=find(TI0>length(t_BF));
            if isempty(badIDX)~=1
                TI0(badIDX)=length(t_BF);   % dummy values--weights will kill it.
                wtemp(badIDX)=0;            % set weight to zero for points out of bounds
            end
            % determine weights
            if ~twoD
                weightTotal    = sum(wtemp,3); % total weigth calculated to divide its contribution from the final beamformed image
            else
                weightTotal    = sum(wtemp,4);
                weightTotal    = sum(squeeze(weightTotal),4);
            end
            weightTotal(weightTotal<1)=1; %avoid devide buy zero later in case no observations at a particular pixel [value becomes 0/1 = 0]

                % uncomment to keep absolute [original] coordinates; otherwise re-zero to start of pulse
                % % % %         allparam.p.xROI=BF.xROI;
                % % % %         if twoD
                % % % %             allparam.p.yROI=BF.yROI; 
                % % % %         end
                % % % %         allparam.p.zROI=BF.zROI;

            subsetIndex_fs_BF = floor(min(TI0(:))):ceil(max(TI0(:)));
            lenIDX=length(subsetIndex_fs_BF);
            %!!---> STILL NEED TO MAKE EVEN  FOR ALL CONDITIONS
            if lenIDX/2~=floor(lenIDX/2) && subsetIndex_fs_BF(end)<length(t_BF)%check if odd and data is big enough
                subsetIndex_fs_BF=subsetIndex_fs_BF(1):subsetIndex_fs_BF(lenIDX)+1; %make even
            end

            TI0      = TI0 - subsetIndex_fs_BF(1) + 1; % subtract the lowest index into the wavefield


            sig = RFdata;                                           % reassign for beamforming and keep original
            if ~twoD
                sig = interp2(BF.aperEx,tUS,sig,BF.aperEx,t_BF,'linear'); % Upsample RF for beamforming
                sig = sig(subsetIndex_fs_BF,:);                         % Select valid time pts for beamforming based on xROI
            else
                sig = interp3(BF.aperEx, tUS, BF.aperEy,sig,BF.aperEx,t_BF,BF.aperEy,'linear',5);
%                         sig = permute(sig,[3 2 1]);
                sig = sig(subsetIndex_fs_BF,:,:);
            end
            % perform beamforming
            if ~twoD
                BFdata = zeros(length(BF.z),length(BF.x));   
            else
                BFdata = zeros(length(BF.z),length(BF.x),length(BF.y));
            end
            for xi = BF.aperEx 	         % use these channels for beamforming based on aquisition
                if ~twoD
                    tmp=sig(:,xi);
                    BFdata = BFdata + tmp(TI0(:,:,xi)).*rhs(:,:,xi);%(1 + DT(:,:,xi)).*wtemp(:,:,xi);
% % %                                         tmp=20*log10(abs(BFdata+1e-15));%-max(abs(im(:)));      % uncomment to show BF steps
% % %                                         imagesc(BF.x,BF.z,tmp-max(tmp(:)),[-25,0]); title(num2str(xi));drawnow;  % uncomment to show BF steps
% % %                                         pause(0.2)
                else
                    for yi = BF.aperEy
                        tmp=sig(:,xi,yi);
                        BFdata = BFdata + tmp(TI0(:,:,:,xi,yi)).*rhs(:,:,:,xi,yi);
                    end
                end
            end
            clear tmp
            BFdata = BFdata./weightTotal;

            if BF.latFiltSmthOn==1
               gaussFilt = fspecial('gaussian',[1 BF.gLat], BF.gLat/4);
               BFdata=imfilter(BFdata,gaussFilt,'replicate','same');
            end 

            % convert to dB
            BFenv = abs(hilbert(BFdata)); BFenvMIN = min(BFenv(:));BFenvMAX = max(BFenv(:));
            BFdB  = 20*log10(BFenv+1e-9); mxVal = max(BFdB(:)); 
            if useREF==1    % choose absolute or relative amplitudes for dB scale
                BFdB = BFdB - refVal;
            else
                BFdB = BFdB - mxVal;
            end
            if twoD
                [~,RFXPlot] = min(abs(RFxMM - PlotXPt));
                [~,RFYPlot] = min(abs(RFyMM - PlotYPt));
                [~,RFZPlot] = min(abs(tUS   - PlotZPt));
                [~,BFXPlot] = min(abs(BF.x - PlotXPt));
                [~,BFYPlot] = min(abs(BF.y - PlotYPt));
                [~,BFZPlot] = min(abs(BF.z - PlotZPt));
                RFMaxXZ = max(max(squeeze(RFdata(:,:,RFYPlot))));
                RFMaxYZ = max(max(squeeze(RFdata(:,RFXPlot,:))));
                RFMaxXY = max(max(squeeze(RFdata(RFZPlot,:,:))));                        
                BFenvMinXZ = min(min(squeeze(BFenv(:,:,BFYPlot))));
                BFenvMinYZ = min(min(squeeze(BFenv(:,BFXPlot,:))));
                BFenvMinXY = min(min(squeeze(BFenv(BFZPlot,:,:))));  
                BFenvMaxXZ = max(max(squeeze(BFenv(:,:,BFYPlot))));
                BFenvMaxYZ = max(max(squeeze(BFenv(:,BFXPlot,:))));
                BFenvMaxXY = max(max(squeeze(BFenv(BFZPlot,:,:))));
                if RFMaxXZ == 0; RFMaxXZ = 1; end
                if RFMaxYZ == 0; RFMaxYZ = 1; end
                if RFMaxXY == 0; RFMaxXY = 1; end
%                         if BFenvMinXZ == 0; BFenvMinXZ = 1; end
%                         if BFenvMinYZ == 0; BFenvMinYZ = 1; end
%                         if BFenvMinXY == 0; BFenvMinXY = 1; end
                if BFenvMaxXZ == BFenvMinXZ; BFenvMaxXZ = BFenvMinXZ+1; end
                if BFenvMaxYZ == BFenvMinYZ; BFenvMaxYZ = BFenvMinYZ+1; end
                if BFenvMaxXY == BFenvMinXY; BFenvMaxXY = BFenvMinXY+1; end
            end
            if GraphIt == 2
                if ~twoD % Start plotting
                    figure(4); title(['Pulse', num2str(x), ' Avg', num2str(u),...
                        ' XPt',num2str(XPt),' YPt',num2str(YPt)])
                    sp1=subplot(2,2,1);
                    han1.P = imagesc(RFxMM,tUS,RFdata,linRngFac*[-RFMAX RFMAX]);
                    xlabel('X [mm]');ylabel('Time [us]'); title(['RF',' Pulse', num2str(x), ' Avg', num2str(u),...
                        ' XPt',num2str(XPt),' YPt',num2str(YPt)],'FontSize',8);
                    colormap(sp1,gray);ylim([tUS(1) tUS(end)]);

                    sp2=subplot(2,2,2);
                    han2.P=imagesc(BF.x,BF.z,BFdata);
                    xlabel('Lateral X [mm]');ylabel('Depth Z [mm]'); title(['BF RF',' Pulse', num2str(x), ' Avg', num2str(u),...
                        ' XPt',num2str(XPt),' YPt',num2str(YPt)],'FontSize',8);
                    colormap(sp2, gray);set(gca,'Tag','lin'); xlim([BF.x(1) BF.x(end)]);ylim([BF.z(1) BF.z(end)]);

                    sp3=subplot(2,2,3);
                    han3.P = imagesc(BF.x,BF.z,BFenv,linRngFac*[BFenvMIN+(1-linRngFac)*BFenvMIN,BFenvMAX-(1-linRngFac)*BFenvMAX]);
                    xlabel('Lateral X [mm]');ylabel('Depth Z [mm]'); title(['BF Image Linear',' Pulse', num2str(x), ' Avg', num2str(u),...
                        ' XPt',num2str(XPt),' YPt',num2str(YPt)],'FontSize',8);
                    colormap(sp3,gray); xlim([BF.x(1) BF.x(end)]);ylim([BF.z(1) BF.z(end)]);

                    sp4=subplot(2,2,4);
                    han4.P = imagesc(BF.x,BF.z,BFdB,dBRng);
                    xlabel('Lateral X [mm]');ylabel('Depth Z [mm]'); title(['BF Image ' num2str(dBRng(1),2) ' to ' num2str(dBRng(2),2) ' dB',' Pulse', num2str(x), ' Avg', num2str(u),...
                        ' XPt',num2str(XPt),' YPt',num2str(YPt)], 'FontSize',8);
                    colormap(sp4,hot); xlim([BF.x(1) BF.x(end)]);ylim([BF.z(1) BF.z(end)]);
                else 
                    % XZ Plane
                    XZ = figure(1); XZ.Name = ['XZ Pulse @ Y=' num2str(PlotYPt) 'mm']; 
                    sp1=subplot(2,2,1);
                    han1.P = imagesc(RFxMM,tUS,squeeze(RFdata(:,:,RFYPlot)),linRngFac*[-RFMaxXZ RFMaxXZ]);
                    xlabel('X [mm]'); ylabel('Time [us]'); title(['RF',' Pulse', num2str(x), ' Avg', num2str(u),...
                        ' XPt',num2str(XPt),' YPt',num2str(YPt)],'FontSize',8);
                    colormap(sp1,gray);ylim([tUS(1) tUS(end)]);

                    sp2=subplot(2,2,2);
                    han2.P=imagesc(BF.x,BF.z,squeeze(BFdata(:,:,BFYPlot)));
                    xlabel('Lateral X [mm]');ylabel('Depth Z [mm]'); title(['BF RF',' Pulse', num2str(x), ' Avg', num2str(u),...
                        ' XPt',num2str(XPt),' YPt',num2str(YPt)],'FontSize',8);
                    colormap(sp2, gray);set(gca,'Tag','lin'); xlim([BF.x(1) BF.x(end)]);ylim([BF.z(1) BF.z(end)]);

                    sp3=subplot(2,2,3);
                    han3.P = imagesc(BF.x,BF.z,squeeze(BFenv(:,:,BFYPlot)),linRngFac*[BFenvMinXZ+(1-linRngFac)*BFenvMinXZ,BFenvMaxXZ-(1-linRngFac)*BFenvMaxXZ]);
                    xlabel('Lateral X [mm]');ylabel('Depth Z [mm]'); title(['BF Image Linear',' Pulse', num2str(x), ' Avg', num2str(u),...
                        ' XPt',num2str(XPt),' YPt',num2str(YPt)],'FontSize',8);
                    colormap(sp3,gray); xlim([BF.x(1) BF.x(end)]);ylim([BF.z(1) BF.z(end)]);

                    sp4=subplot(2,2,4);
                    han4.P = imagesc(BF.x,BF.z,squeeze(BFdB(:,:,BFYPlot)),dBRng);
                    xlabel('Lateral X [mm]');ylabel('Depth Z [mm]'); title(['BF Image ' num2str(dBRng(1),2) ' to ' num2str(dBRng(2),2) ' dB',' Pulse', num2str(x), ' Avg', num2str(u),...
                        ' XPt',num2str(XPt),' YPt',num2str(YPt)], 'FontSize',8);
                    colormap(sp4,hot); xlim([BF.x(1) BF.x(end)]);ylim([BF.z(1) BF.z(end)]);
                    % YZ Plane
                    YZ = figure(2); YZ.Name = ['YZ Pulse @ X=' num2str(PlotXPt) 'mm']; 
                    sp5=subplot(2,2,1);
                    han5.P = imagesc(RFyMM,tUS,squeeze(RFdata(:,RFXPlot,:)),linRngFac*[-RFMaxYZ RFMaxYZ]);
                    xlabel('Y [mm]');ylabel('Time [us]'); title(['RF',' Pulse', num2str(x), ' Avg', num2str(u),...
                        ' XPt',num2str(XPt),' YPt',num2str(YPt)],'FontSize',8);
                    colormap(sp5,gray);ylim([tUS(1) tUS(end)]);

                    sp6=subplot(2,2,2);
                    han6.P=imagesc(BF.y,BF.z,squeeze(BFdata(:,BFXPlot,:)));
                    xlabel('Elevational Y [mm]');ylabel('Depth Z [mm]'); title(['BF RF',' Pulse', num2str(x), ' Avg', num2str(u),...
                        ' XPt',num2str(XPt),' YPt',num2str(YPt)],'FontSize',8);
                    colormap(sp6, gray);set(gca,'Tag','lin'); xlim([BF.y(1) BF.y(end)]);ylim([BF.z(1) BF.z(end)]);

                    sp7=subplot(2,2,3);
                    han7.P = imagesc(BF.y,BF.z,squeeze(BFenv(:,BFXPlot,:)),linRngFac*[BFenvMinYZ+(1-linRngFac)*BFenvMinYZ,BFenvMaxYZ-(1-linRngFac)*BFenvMaxYZ]);
                    xlabel('Elevational Y [mm]');ylabel('Depth Z [mm]'); title(['BF Image Linear',' Pulse', num2str(x), ' Avg', num2str(u),...
                        ' XPt',num2str(XPt),' YPt',num2str(YPt)],'FontSize',8);
                    colormap(sp7,gray); xlim([BF.y(1) BF.y(end)]);ylim([BF.z(1) BF.z(end)]);

                    sp8=subplot(2,2,4);
                    han8.P = imagesc(BF.y,BF.z,squeeze(BFdB(:,BFXPlot,:)),dBRng);
                    xlabel('Elevational Y [mm]');ylabel('Depth Z [mm]'); title(['BF Image ' num2str(dBRng(1),2) ' to ' num2str(dBRng(2),2) ' dB',' Pulse', num2str(x), ' Avg', num2str(u),...
                        ' XPt',num2str(XPt),' YPt',num2str(YPt)], 'FontSize',8);
                    colormap(sp8,hot); xlim([BF.y(1) BF.y(end)]);ylim([BF.z(1) BF.z(end)]);                        
                    % XY Plane
                    XY = figure(3);XY.Name = ['XY Pulse @ Z=' num2str(PlotZPt) 'mm']; 
                    sp9=subplot(2,2,1);
                    han9.P = imagesc(RFxMM,RFyMM,squeeze(RFdata(RFZPlot,:,:)),linRngFac*[-RFMaxXY RFMaxXY]);
                    xlabel('X [mm]');ylabel('Y [mm]'); title(['RF',' Pulse', num2str(x), ' Avg', num2str(u),...
                        ' XPt',num2str(XPt),' YPt',num2str(YPt)],'FontSize',8);
                    colormap(sp9,gray);

                    sp10=subplot(2,2,2);
                    han10.P=imagesc(BF.x,BF.y,squeeze(BFdata(BFZPlot,:,:)));
                    xlabel('Lateral X [mm]');ylabel('Elevational Y [mm]'); title(['BF RF',' Pulse', num2str(x), ' Avg', num2str(u),...
                        ' XPt',num2str(XPt),' YPt',num2str(YPt)],'FontSize',8);
                    colormap(sp10, gray);set(gca,'Tag','lin'); xlim([BF.x(1) BF.x(end)]);ylim([BF.y(1) BF.y(end)]);

                    sp11=subplot(2,2,3);
                    han11.P = imagesc(BF.x,BF.y,squeeze(BFenv(BFZPlot,:,:)),linRngFac*[BFenvMinXY+(1-linRngFac)*BFenvMinXY,BFenvMaxXY-(1-linRngFac)*BFenvMaxXY]);
                    xlabel('Lateral X [mm]');ylabel('Elevational Y [mm]'); title(['BF Image Linear',' Pulse', num2str(x), ' Avg', num2str(u),...
                        ' XPt',num2str(XPt),' YPt',num2str(YPt)],'FontSize',8);
                    colormap(sp11,gray); xlim([BF.x(1) BF.x(end)]);ylim([BF.y(1) BF.y(end)]);

                    sp12=subplot(2,2,4);
                    han12.P = imagesc(BF.x,BF.y,squeeze(BFdB(BFZPlot,:,:)),dBRng);
                    xlabel('Lateral X [mm]');ylabel('Elevational Y [mm]'); title(['BF Image ' num2str(dBRng(1),2) ' to ' num2str(dBRng(2),2) ' dB',' Pulse', num2str(x), ' Avg', num2str(u),...
                        ' XPt',num2str(XPt),' YPt',num2str(YPt)], 'FontSize',8);
                    colormap(sp12,hot); xlim([BF.x(1) BF.x(end)]);ylim([BF.y(1) BF.y(end)]);                     

                end
            elseif GraphIt == 1
                if strcmpi(MModeType,'Vertical')
                    [~, BF_XPt] = min(abs(BF.x-lat_mm(XPt)));
                    if twoD
                        [~, BF_YPt]         = min(abs(BF.y-ele_mm(YPt)));
                        MModeRF(:,x,u)      = RFdata(:,Chantaje_X,Chantaje_Y);
                        MModeBFRF(:,x,u)    = BFdata(:,BF_XPt,BF_YPt); % [Fast Time, PE Slow Time (i.e., Pulse No), Average]
                        MModeBFEnv(:,x,u)   = BFenv(:,BF_XPt,BF_YPt);
                        MModeBFdB(:,x,u)    = BFdB(:,BF_XPt,BF_YPt);
                    else
                        MModeRF(:,x,u)      = RFdata(:,Chantaje_X);
                        MModeBFRF(:,x,u)    = BFdata(:,BF_XPt); % [Fast Time, PE Slow Time (i.e., Pulse No), Average]
                        MModeBFEnv(:,x,u)   = BFenv(:,BF_XPt);
                        MModeBFdB(:,x,u)    = BFdB(:,BF_XPt);                                                            
                    end  
                elseif strcmpi(MModeType,'Angular')
                    if twoD
                        [ThetaBF,PhiBF,RBF] = cart2sph(BF.x,BF.z,BF.y);
                        [~, PuntoTheta] = min(abs(ThetaBF-Theta));
                        [~, PuntoPhi]   = min(abs(PhiBF-Phi));
                        [BFPol.X, BFPol.Z, BFPol.Y] = sph2cart(ThetaBF(PuntoTheta),PhiBF(PuntoPhi),RBF);
                        for i = 1:numel(BFPol.X)
                            [~, MinX] = min(abs(BF.x-BFPol.X(i)));
                            [~, MinY] = min(abs(BF.y-BFPol.Y(i)));
                            [~, MinZ] = min(abs(BF.z-BFPol.Z(i)));
                            MModeBFRF(i,x,u) = BFdata(MinZ,MinX,MinY);
                            MModeBFEnv(i,x,u) = BFenv(MinZ,MinX,MinY);
                            MModeBFdB(i,x,u) = BFdB(MinZ,MinX,MinY);
                        end                        
                        MModeRF(:,x,u)    = RFdata(:,Chantaje_X,Chantaje_Y);
% %                         MModeBFRF(:,x,u) = BFdata(:,PuntoTheta,PuntoPhi); % [Fast Time, PE Slow Time (i.e., Pulse No), Average]
% %                         MModeBFEnv(:,x,u) = BFenv(:,PuntoTheta,PuntoPhi);
% %                         MModeBFdB(:,x,u)  = BFdB(:,PuntoTheta,PuntoPhi);
                    else
                        [ThetaBF,RBF]     = cart2pol(BF.x,BF.z);
                        [~, PuntoTheta]   = min(abs(ThetaBF-Theta));
                        [BFPol.X, BFPol.Z] = sph2cart(ThetaBF(PuntoTheta),RBF);
                        for i = 1:numel(BFPol.X)
                            [~, MinX] = min(abs(BF.x-BFPol.X(i)));
                            [~, MinZ] = min(abs(BF.z-BFPol.Z(i)));
                            MModeBFRF(i,x,u) = BFdata(MinZ,MinX);
                            MModeBFEnv(i,x,u) = BFenv(MinZ,MinX);
                            MModeBFdB(i,x,u) = BFdB(MinZ,MinX);
                        end
                        MModeRF(:,x,u)    = RFdata(:,Chantaje_X);
% %                         MModeBFRF(:,x,u)  = BFdata(:,PuntoTheta); % [Fast Time, PE Slow Time (i.e., Pulse No), Average]
% %                         MModeBFEnv(:,x,u) = BFenv(:,PuntoTheta);
% %                         MModeBFdB(:,x,u)  = BFdB(:,PuntoTheta);
                    end  
                end
            end
        end
    end
    if GraphIt == 1
        MaxMModeRF = max(MModeRF(:));
        MinMModeBFEnv = min(MModeBFEnv(:));               
        MaxMModeBFEnv = max(MModeBFEnv(:));
        figure(4);
        if ~twoD
            %XZ Plane Plot
            sp1=subplot(2,3,1);        
            han1.P = imagesc(BF.x,BF.z,BFdB,dBRng);
            xlabel('Lateral X [mm]');ylabel('Depth Z [mm]'); title(['BF Image ' num2str(dBRng(1),2) ' to ' num2str(dBRng(2),2) ' dB',' Pulse', num2str(x), ' Avg', num2str(u),...
                ' XPt',num2str(XPt),' YPt',num2str(YPt)], 'FontSize',8);
            colormap(sp1,gray); xlim([BF.x(1) BF.x(end)]);ylim([BF.z(1) BF.z(end)]); hold on;
        else
            %XZ Plane Plot
            sp1=subplot(2,3,1);        
            han1.P = imagesc(BF.x,BF.z,squeeze(BFdB(:,:,BFYPlot)),dBRng);
            xlabel('Lateral X [mm]');ylabel('Depth Z [mm]'); title(['BF Image ' num2str(dBRng(1),2) ' to ' num2str(dBRng(2),2) ' dB',' Pulse', num2str(x), ' Avg', num2str(u),...
                ' XPt',num2str(XPt),' YPt',num2str(YPt)], 'FontSize',8);
            colormap(sp1,gray); xlim([BF.x(1) BF.x(end)]);ylim([BF.z(1) BF.z(end)]); hold on;
            %YZ Plane Plot
            sp2=subplot(2,3,2);
            han2.P = imagesc(BF.y,BF.z,squeeze(BFdB(:,BFXPlot,:)),dBRng);
            xlabel('Elevational Y [mm]');ylabel('Depth Z [mm]'); title(['BF Image ' num2str(dBRng(1),2) ' to ' num2str(dBRng(2),2) ' dB',' Pulse', num2str(x), ' Avg', num2str(u),...
            ' XPt',num2str(XPt),' YPt',num2str(YPt)], 'FontSize',8);
            colormap(sp2,gray); xlim([BF.y(1) BF.y(end)]);ylim([BF.z(1) BF.z(end)]);  hold on;
            %XY Plane PLot
            sp3=subplot(2,3,3);
            han3.P = imagesc(BF.x,BF.y,squeeze(BFdB(BFZPlot,:,:)),dBRng);
            xlabel('Lateral X [mm]');ylabel('Elevational Y [mm]'); title(['BF Image ' num2str(dBRng(1),2) ' to ' num2str(dBRng(2),2) ' dB',' Pulse', num2str(x), ' Avg', num2str(u),...
                ' XPt',num2str(XPt),' YPt',num2str(YPt)], 'FontSize',8);
            colormap(sp3,gray); xlim([BF.x(1) BF.x(end)]);ylim([BF.y(1) BF.y(end)]); hold on;                   
        end
        if strcmpi(MModeType,'Angular')
%             [PlotLineX, PlotLineY, PlotLineZ] = sph2cart(ThetaBF(PuntoTheta)*ones(numel(RBF),1),PhiBF(PuntoPhi)*ones(numel(RBF),1),RBF');
            plot(sp1,BFPol.X,BFPol.Z,'g--'); hold off;
            plot(sp2,BFPol.Y,BFPol.Z,'g--'); hold off;
            text(sp3,PuntoTheta,PuntoPhi,'go');
            sp4=subplot(2,3,4);
            han4.P=imagesc(SlowTime,RBF,MModeBFRF);
            xlabel('Slow Time [ms]');ylabel('Axial [mm]'); title(['BF RF M-Mode',' Avg', num2str(u),...
                ' XPt=',num2str(XPt),' YPt=',num2str(YPt)],'FontSize',8);
            colormap(sp4, gray);set(gca,'Tag','lin'); xlim([SlowTime(1) SlowTime(end)]);ylim([RBF(1) RBF(end)]);

            sp5=subplot(2,3,5);
            han5.P = imagesc(SlowTime,RBF,MModeBFEnv,linRngFac*[MinMModeBFEnv+(1-linRngFac)*MinMModeBFEnv,MaxMModeBFEnv-(1-linRngFac)*MaxMModeBFEnv]);
            xlabel('Slow Time [ms]');ylabel('Axial [mm]'); title(['BF Linear M-Mode', ' Avg', num2str(u),...
                ' XPt=',num2str(XPt),' YPt=',num2str(YPt)],'FontSize',8);
            colormap(sp5,gray); xlim([SlowTime(1) SlowTime(end)]);ylim([RBF(1) BF.z(end)]);

            sp6=subplot(2,3,6);
            han6.P = imagesc(SlowTime,RBF,MModeBFdB,dBRng);
            xlabel('Slow Time [ms]');ylabel('Axial [mm]'); title(['BF dB M-Mode ' num2str(dBRng(1),2) ' to ' num2str(dBRng(2),2) ' dB',...
                ' Avg', num2str(u),' XPt=',num2str(XPt),' YPt=',num2str(YPt)], 'FontSize',8);
            colormap(sp6,gray); xlim([SlowTime(1) SlowTime(end)]);ylim([BF.z(1) BF.z(end)]);       
        elseif strcmpi(MModeType,'Vertical')
            plot(sp1,XPt,'g--'); hold off;
            plot(sp2,YPt,'g--'); hold off;
            text(sp3,XPt,YPt,'go'); hold off;
            sp4=subplot(2,3,4);
            han4.P=imagesc(SlowTime,RBF,MModeBFRF);
            xlabel('Slow Time [ms]');ylabel('Depth Z [mm]'); title(['BF RF M-Mode',' Avg', num2str(u),...
                ' XPt=',num2str(XPt),' YPt=',num2str(YPt)],'FontSize',8);
            colormap(sp4, gray);set(gca,'Tag','lin'); xlim([SlowTime(1) SlowTime(end)]);ylim([RBF(1) RBF(end)]);

            sp5=subplot(2,3,5);
            han5.P = imagesc(SlowTime,RBF,MModeBFEnv,linRngFac*[MinMModeBFEnv+(1-linRngFac)*MinMModeBFEnv,MaxMModeBFEnv-(1-linRngFac)*MaxMModeBFEnv]);
            xlabel('Slow Time [ms]');ylabel('Depth Z [mm]'); title(['BF Linear M-Mode', ' Avg', num2str(u),...
                ' XPt=',num2str(XPt),' YPt=',num2str(YPt)],'FontSize',8);
            colormap(sp5,gray); xlim([SlowTime(1) SlowTime(end)]);ylim([RBF(1) BF.z(end)]);

            sp6=subplot(2,3,6);
            han6.P = imagesc(SlowTime,RBF,MModeBFdB,dBRng);
            xlabel('Slow Time [ms]');ylabel('Depth Z [mm]'); title(['BF dB M-Mode ' num2str(dBRng(1),2) ' to ' num2str(dBRng(2),2) ' dB',...
                ' Avg', num2str(u),' XPt=',num2str(XPt),' YPt=',num2str(YPt)], 'FontSize',8);
            colormap(sp6,gray); xlim([SlowTime(1) SlowTime(end)]);ylim([BF.z(1) BF.z(end)]);       
        end
    end
end
disp([num2str(toc(BMFMTime)), ' Seconds to Beamform @ XPt=', num2str(XPt), ' YPt=' num2str(YPt)]);
cd ..\ExpData
end