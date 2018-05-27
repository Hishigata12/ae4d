function param = read_ucsdi_info(InfoFile)
% function param = read_ucsdi_info(InfoFile)
% load configuration data of ucsdi scan recorded using the LabVIEW DAQ standalone program
% Inputs:
%   InfoFile: scan & data acquisition configuration file (*_info.dat)
% Outputs:
%   param:  the structure of scan & daq parameters.

param = struct;

narginchk(1,1);

try
    fid = fopen(InfoFile,'rb');
    n = fread(fid,1,'int32');
catch ME
    uiwait(errordlg(ME.message,mfilename('fullpath')));
    return;
end

if(n~=10)
    fclose(fid);
    uiwait(errordlg([fname,' is not a valid UCSDI info file.'],mfilename,'modal'));
    return;
end

file_type = fgets(fid,n);
if(strcmpi(file_type,'ucsdi_info') ~= 1)
    fclose(fid);
    uiwait(errordlg([fname,' is not a valid UCSDI info file.'],mfilename,'modal'));
    return;
end

nver = fread(fid,1,'int32');
param.version = nver;

switch(nver)
    case {140908,160714}
        param = read_info_v140908(fid,param);
    case 161108
        [~,param] = read_info_v161108(fid,param);
    case 170103
        [~,param] = read_info_v170103(fid,param);
    case 170427
        [~,param] = read_info_v170427(fid,param);
    otherwise
        param = read_info_default(fid,param);
end
fclose(fid);


[f_path, f_name] = fileparts(InfoFile);
param.f_path = f_path;
param.f_root = regexprep(f_name,'_info[.dat]*$','');
end

%%
%%
function param = read_info_v140908(fid,param)
daq.SingleDev = true;
n = fread(fid,1,'int32');
daq.HFdaq.DeviceName = fgets(fid,n);
n = fread(fid,1,'int32');
daq.HFdaq.HFchannelStr = fgets(fid,n);
daq.HFdaq.fs_MHz = fread(fid,1,'float64');
daq.HFdaq.pts = fread(fid,1,'int32');
tmpdata = fread(fid,4,'float64');
daq.HFdaq.pulseRepRate_Hz = tmpdata(1);
daq.HFdaq.Gain = tmpdata(2);
daq.HFdaq.vertRange = tmpdata(3);
daq.HFdaq.vertRangeAdj = tmpdata(3);
daq.HFdaq.TrigDelay_us = tmpdata(4);

daq.LFdaq.LFchannelStr = fgets(fid,fread(fid,1,'int32'));
tmpdata = fread(fid,3,'float64');
daq.LFdaq.fs_KHz = tmpdata(1);
daq.LFdaq.fs_Hz = tmpdata(1)*1000;
daq.LFdaq.vertRangeAdj = tmpdata(2);
daq.LFdaq.duration_s = tmpdata(3)/1000;
daq.duration_ms = tmpdata(3);
daq.HFdaq.duration_ms = tmpdata(3);
daq.NBurstAvg = fread(fid,1,'int32');
daq.burstRepRate_Hz = fread(fid,1,'float64');
daq.LFdaq.channels = str2num(daq.LFdaq.LFchannelStr);
daq.LFdaq.pts = round(daq.LFdaq.fs_KHz * daq.duration_ms);
daq.HFdaq.Coupling = fread(fid,1,'int32');
daq.LFdaq.NBurstAvg = daq.NBurstAvg;

daq.HFdaq.channels = str2num(daq.HFdaq.HFchannelStr);
daq.HFdaq.Nchannels = length(daq.HFdaq.channels);
daq.HFdaq.NoBurstTriggers = round(daq.duration_ms * daq.HFdaq.pulseRepRate_Hz/1000);
daq.PreTriggerDelay_ms = 0;

param.daq = daq;
sAxes = {'X','Y','Z'};
lengthUnits = {'mm','in'};

velmex.scanDis1 = fread(fid,1,'float64');
velmex.nSteps1 = fread(fid,1,'int32');
velmex.FastAxis = sAxes{fread(fid,1,'uint16')+1};

velmex.scanDis2 = fread(fid,1,'float64');
velmex.nSteps2 = fread(fid,1,'int32');
velmex.SlowAxis = sAxes{fread(fid,1,'uint16')+1};

if(strcmpi(velmex.FastAxis,'X'))
    velmex.XDist = velmex.scanDis1;
    velmex.XNStep = velmex.nSteps1;
    velmex.YDist = velmex.scanDis2;
    velmex.YNStep = velmex.nSteps2;
    velmex.YScan = false;
elseif(strcmpi(velmex.FastAxis,'Y'))
    velmex.XDist = velmex.scanDis2;
    velmex.XNStep = velmex.nSteps2;
    velmex.YDist = velmex.scanDis1;
    velmex.YNStep = velmex.nSteps1;
    velmex.YScan = true;
end


velmex.ZNStep = 1;
velmex.ZDist = 0;
velmex.unit = lengthUnits{fread(fid,1,'uint16')+1};
param.velmex = velmex;

hp.TCPIPaddress = fgets(fid,fread(fid,1,'int32'));
CurShape = {'Sin','Square','Pulse','ECG','ARB'};
hp.shape = CurShape{fread(fid,1,'uint16')+1};
tmpdata = fread(fid,3,'float64');
hp.freq = tmpdata(1);
hp.volt = tmpdata(2);
hp.imp = tmpdata(3);

val = fread(fid,1,'int8');
hp.burstmode = logical(val);

TrigSrc = {'IMM','EXT'};
hp.trig	 = TrigSrc{fread(fid,1,'uint16')+1};
hp.cycles = fread(fid,1,'int32');
hp.posphase = fread(fid,1,'float64');
hp.period = fread(fid,1,'float64');
param.hp = hp;

clrmaps = {'gray','jet','hot','hotcold'};
sp.ColorMap = clrmaps{fread(fid,1,'uint16')+1};
sp.AEChanSel = fread(fid,1,'uint16');
sp.DemodFreq = fread(fid,1,'float64');
tmpdata = logical(fread(fid,4,'int8'));
sp.HFMatchFiltOn = tmpdata(1);
sp.HFHannWinOn = tmpdata(2);
sp.LFHannWinOn = tmpdata(3);
sp.ShowEnv = tmpdata(4);
sp.HFFiltWin = str2num(fgets(fid,fread(fid,1,'int32')));
sp.LFFiltWin = str2num(fgets(fid,fread(fid,1,'int32')));

filterType = {'None','HannWin'};
sp.FilterType = filterType{fread(fid,1,'uint16')+1};
% n = fread(fid,1,'int32');
sp.AE_tmRng = fread(fid,[1,2],'float64');
sp.LF_tmRng = round(fread(fid,[1,2],'float64'));
sp.SndSpeed = 1.48;
sp.zRange = round(sort(sp.AE_tmRng)*sp.SndSpeed);
sp.ImageAutoScale = false;
sp.AEOpt = 1;
sp.PEdBRng = [-30,0];
sp.AEdBRng = [-10,0];
sp.SaveFigs = false;

param.sp = sp;
% fclose(fid);
end

%%
function param = read_info_default(fid,param)
nver = param.version;
daq.SingleDev = false;
n = fread(fid,1,'int32');

daq.HFdaq.HFchannelStr = fgets(fid,n);
daq.HFdaq.fs_MHz = fread(fid,1,'float64');
daq.HFdaq.pts = fread(fid,1,'int32');
tmpdata = fread(fid,4,'float64');
daq.HFdaq.TrigDelay_us = tmpdata(1);
daq.HFdaq.Gain = tmpdata(2);
daq.HFdaq.vertRange = tmpdata(3);
daq.HFdaq.vertRangeAdj = tmpdata(3);
daq.HFdaq.pulseRepRate_Hz = tmpdata(4);

daq.LFdaq.LFchannelStr = fgets(fid,fread(fid,1,'int32'));
tmpdata = fread(fid,4,'float64');
daq.LFdaq.fs_KHz = tmpdata(1);
daq.LFdaq.fs_Hz = tmpdata(1)*1000;
daq.LFdaq.vertRangeAdj = tmpdata(2);
daq.LFdaq.duration_s = tmpdata(3)/1000;
daq.duration_ms = tmpdata(3);
daq.HFdaq.duration_ms = tmpdata(3);
daq.burstRepRate_Hz = tmpdata(4);
daq.LFdaq.channels = str2num(daq.LFdaq.LFchannelStr);
daq.NBurstAvg = fread(fid,1,'int32');
daq.LFdaq.pts = daq.LFdaq.fs_KHz * daq.duration_ms;
daq.LFdaq.NBurstAvg = daq.NBurstAvg;

daq.HFdaq.channels = str2num(daq.HFdaq.HFchannelStr);
daq.HFdaq.Nchannels = length(daq.HFdaq.channels)*2;
daq.HFdaq.NoBurstTriggers = round(daq.duration_ms * daq.HFdaq.pulseRepRate_Hz/1000);
daq.PreTriggerDelay_ms = 0; 
param.daq = daq;
sAxes = {'X','Y','Z'};
lengthUnits = {'mm','in'};

velmex.scanDis1 = fread(fid,1,'single');
velmex.nSteps1 = fread(fid,1,'int32');
velmex.FastAxis = sAxes{fread(fid,1,'uint16')+1};

velmex.scanDis2 = fread(fid,1,'single');
velmex.nSteps2 = fread(fid,1,'int32');
velmex.SlowAxis = sAxes{fread(fid,1,'uint16')+1};

if(strcmpi(velmex.FastAxis,'X'))
    velmex.XDist = velmex.scanDis1;
    velmex.XNStep = velmex.nSteps1;
    velmex.YDist = velmex.scanDis2;
    velmex.YNStep = velmex.nSteps2;
    velmex.YScan = false;
elseif(strcmpi(velmex.FastAxis,'Y'))
    velmex.XDist = velmex.scanDis2;
    velmex.XNStep = velmex.nSteps2;
    velmex.YDist = velmex.scanDis1;
    velmex.YNStep = velmex.nSteps1;
    velmex.YScan = true;
end

velmex.ZNStep = 1;
velmex.ZDist = 0;
velmex.unit = lengthUnits{fread(fid,1,'uint16')+1};
param.velmex = velmex;

hp.TCPIPaddress = fgets(fid,fread(fid,1,'int32'));
CurShape = {'Sin','Square','Pulse','ECG'};
hp.shape = CurShape{fread(fid,1,'uint16')+1};
tmpdata = fread(fid,3,'float64');
hp.freq = tmpdata(1);
hp.volt = tmpdata(2);
hp.imp = tmpdata(3);

val = fread(fid,1,'int8');
hp.burstmode = logical(val);

TrigSrc = {'IMM','EXT'};
hp.trig	 = TrigSrc{fread(fid,1,'uint16')+1};
hp.cycles = fread(fid,1,'int32');
hp.posphase = fread(fid,1,'float64');
hp.period = fread(fid,1,'float64');
param.hp = hp;

clrmaps = {'gray','jet','hot','hotcold'};
sp.ColorMap = clrmaps{fread(fid,1,'uint16')+1};
sp.AEChanSel = fread(fid,1,'uint16');
sp.DemodFreq = fread(fid,1,'float64');
sp.HFHannWinOn = logical(fread(fid,1,'int8'));
if(nver == 1140123)
    StrTypeFiltParam = false;
    nLen = fread(fid,1,'int32');
    if(nLen > 0 && nLen < 100)  % length of the string of HF filter parameter
        HFFiltWin = str2num(fgets(fid,nLen));
        StrTypeFiltParam = (min(HFFiltWin >= 0)) && (max(HFFiltWin) <= daq.HFdaq.fs_MHz) && ...
            min(diff(HFFiltWin) >= 0);
        if(StrTypeFiltParam)
            sp.HFFiltWin = HFFiltWin;
        else
            fseek(fid,-nLen,'cof');
        end
    end
    
    if(~StrTypeFiltParam)
        tmp = fread(fid,[1,2],'float64');
        sp.HFFiltWin = [tmp(1) - tmp(2)/2, tmp(1), tmp(1), tmp(1) + tmp(2)/2];
    end
else
    sp.HFFiltWin = str2num(fgets(fid,fread(fid,1,'int32')));
end
sp.LFHannWinOn = logical(fread(fid,1,'int8'));
if(nver == 1140123)
    if(StrTypeFiltParam)
        sp.LFFiltWin = str2num(fgets(fid,fread(fid,1,'int32')));
    else
        tmp = fread(fid,[1,2],'float64');
        sp.LFFiltWin = [tmp(1) - tmp(2)/2, tmp(1), tmp(1), tmp(1) + tmp(2)/2];
    end
else
    sp.LFFiltWin = str2num(fgets(fid,fread(fid,1,'int32')));
end

sp.HFMatchFiltOn = logical(fread(fid,1,'int8'));

filterType = {'None','HannWin'};
sp.FilterType = filterType{fread(fid,1,'uint16')+1};
n = fread(fid,1,'int32');
sp.AE_tmRng = fread(fid,[1,n],'float64');
sp.SndSpeed = 1.48;
if(isempty(sp.AE_tmRng))
    sp.AE_tmRng = [40,90]/sp.SndSpeed;
end


sp.zRange = round(sort(sp.AE_tmRng)*sp.SndSpeed);
sp.ImageAutoScale = false;
sp.AEOpt = 1;
sp.PEdBRng = [-30,0];
sp.AEdBRng = [-10,0];
sp.SaveFigs = false;
sp.LF_tmRng = [0,daq.duration_ms];

param.sp = sp;
% fclose(fid);
end


function [data, parm]  = read_info_v161108(fid,parm)

data = [];
% parm = struct;

parm.daq.SingleDev = true;
n = fread(fid,1,'int32');
parm.daq.HFdaq.DeviceName = fgets(fid,n);
n = fread(fid,1,'int32');
parm.daq.HFdaq.HFchannelStr = fgets(fid,n);
parm.daq.HFdaq.channels = str2num(parm.daq.HFdaq.HFchannelStr);
parm.daq.HFdaq.Nchannels = length(parm.daq.HFdaq.channels);
parm.daq.HFdaq.fs_MHz = fread(fid,1,'float64');       % HFFS<DBL>:  30
parm.daq.HFdaq.pts = fread(fid,1,'int32');           % HFSamples<I32>:  2048

% PulseRate<DBL>:  2000
% HFDigitGain<DBL>:  1
% HFVertRng<DBL>:  1
% HFTrigDelay<DBL>:  0

vals = fread(fid,[1,4],'float64');
parm.daq.HFdaq.pulseRepRate_Hz = vals(1);
parm.daq.HFdaq.Gain = vals(2);
parm.daq.HFdaq.vertRange = vals(3);
parm.daq.HFdaq.vertRangeAdj = parm.daq.HFdaq.vertRange;
parm.daq.HFdaq.TrigDelay_us = vals(4);

% LFChan<String>:  17
parm.daq.LFdaq.LFchannelStr = fgets(fid,fread(fid,1,'int32'));
parm.daq.LFdaq.channels = str2num(parm.daq.LFdaq.LFchannelStr);

% LFFS<DBL>:  10
% LFVertRng<DBL>:  5
% Duration<DBL>:  30
vals = fread(fid,[1,3],'float64');
parm.daq.LFdaq.fs_KHz = vals(1);
parm.daq.LFdaq.fs_Hz = parm.daq.LFdaq.fs_KHz * 1000;
parm.daq.LFdaq.LFVertRng = vals(2);
parm.daq.HFdaq.duration_ms = vals(3);
parm.daq.LFdaq.duration_s = vals(3)/1000;

parm.daq.LFdaq.pts = parm.daq.LFdaq.fs_KHz * parm.daq.LFdaq.duration_s*1000;

parm.daq.NBurstAvg = fread(fid,1,'int32'); % NumOfAvg<I32>:  50
parm.daq.LFdaq.NBurstAvg = parm.daq.NBurstAvg;
parm.daq.burstRepRate_Hz = parm.daq.HFdaq.pulseRepRate_Hz;
parm.daq.PreTriggerDelay_ms = 0;
parm.daq.duration_ms = parm.daq.HFdaq.duration_ms;

parm.daq.BurstRate = fread(fid,1,'float64');% BurstRate<DBL>:  3

parm.daq.HFdaq.Coupling = fread(fid,1,'int32'); % HFCoupling<I32>:  0
parm.daq.NumAcq = fread(fid,1,'int32'); % NumAcq<I32>:  21
parm.daq.HFdaq.NoBurstTriggers = round(parm.daq.duration_ms * parm.daq.burstRepRate_Hz/1000);

% #####PacingSignal#####
% IPAddress<String>:  192.168.1.2
% Shape<EW>: Sin
% Freq (Hz)<DBL>:  200
% Ampl (V)<DBL>:  5
% Load Imp<DBL>:  50
% BurstMode<Boolean>:  1
% TrigSrc<EW>: IMM
% Cycles<I32>:  3
% StartPha<DBL>:  0
% Period (ms)<DBL>:  100
% Enable<Boolean>:  1
% #####END OF PacingSignal#####

parm.hp.TCPIPaddress = fgets(fid,fread(fid,1,'int32'));
shape_idx = fread(fid,1,'uint16');
PulseShape = {'Sin','Square','Pulse','ECG','ARB'};
parm.hp.shape = PulseShape{shape_idx+1};

vals = fread(fid,[1,3],'float64');
parm.hp.freq = vals(1);
parm.hp.volt = vals(2);
parm.hp.imp = vals(3);

parm.hp.burstmode = logical(fread(fid,1,'char'));

trig_src = {'IMM','EXT'};
trig_idx = fread(fid,1,'uint16');

parm.hp.trig = trig_src{trig_idx+1};
parm.hp.cycles = fread(fid,1,'int32');

parm.hp.posphase = fread(fid,1,'double');
parm.hp.BurstPeriod = fread(fid,1,'double');
parm.hp.PacingEnable = logical(fread(fid,1,'char'));
parm.hp.period = [];

%% velmex
parm.velmex.XDist = fread(fid,1,'double');

if isempty(parm.velmex.XDist)
    parm.velmex.XDist = 1;
    parm.velmex.XNStep = parm.daq.NumAcq;
    parm.velmex.FastAxis = 'X';
    parm.velmex.YDist = 0;
    parm.velmex.YNStep = 1;
    parm.velmex.SlowAxis = 'Y';
else
    
    parm.velmex.XNStep = fread(fid,1,'int32');
    if parm.velmex.XNStep < 1 || parm.velmex.XNStep > 1e4
        fclose(fid);
        error('parm.velmex.XNStep is invalid.');
    end
    
    
    val = fread(fid,1,'int16');
    
    parm.velmex.FastAxis = char('X'+val);
    
    parm.velmex.YDist = fread(fid,1,'double');
    parm.velmex.YNStep = fread(fid,1,'int32');
    parm.velmex.SlowAxis = char('X'+fread(fid,1,'int16'));
end

parm.velmex.ZDist = 0;
parm.velmex.ZNStep = 1;
parm.YScan = false;


% parm.velmex = struct('FastAxis','X','SlowAxis','Y',...
%     'XDist',1,'XNStep',parm.daq.NumAcq,...
%     'YDist',0,'YNStep',1,...
%     'ZDist',0,'ZNStep',1,...
%     'YScan',false);

% parm.sp = struct('AEChanSel', 1,...
%     'AEGain', 180,...
%     'AEOpt', 1,...
%     'AE_tmRng', [17 37],...
%     'AEdBRng', [-15 0],...
%     'AEopt', 1,...
%     'Colormap', 'hot',...
%     'DemodFreq', 2.2000,...
%     'FilterType', 'HannWin',...
%     'FrameRate', 15,...
%     'FrameSelect', 1,...
%     'HFFiltWin', [1.2000 2 2.5000 3.3000],...
%     'HFHannWinOn', 1,...
%     'HFMatchFiltOn', 0,...
%     'ImageAutoScale', 0,...
%     'ImageScale', 1,...
%     'LFFiltWin', [0.1000 0.1500 0.2500 0.3000],...
%     'LFGain', 2,...
%     'LFHannWinOn', 1,...
%     'LFImp', 1000,...
%     'LFMatchFiltOn', 0,...
%     'LF_tmRng', [0 44],...
%     'MedFiltOn', 0,...
%     'MedFiltSize', [],...
%     'PE_tmRng', [34 74],...
%     'PEdBRng', [-30 0],...
%     'PlaneChoice', 1,...
%     'PulsComp', 0,...
%     'SaveFigs', 0,...
%     'ScanPt', 14,...
%     'ShowEnv', 0,...
%     'SlicePos', 1,...
%     'SndSpeed', 1.4800,...
%     'XPos', 14,...
%     'YPos', 1,...
%     'bSlowTimeMatch', 0,...
%     'bmFrames', 1,...
%     'nsRng', [50 70],...
%     'szAverageRange', [],...
%     'xdcr', 'P4-2V',...
%     'yAxRev', 0,...
%     'zRange', [25 55]);

end

function [data, parm]  = read_info_v170103(fid,parm)

data = [];
% parm = struct;

parm.daq.SingleDev = true;
n = fread(fid,1,'int32');
parm.daq.HFdaq.DeviceName = fgets(fid,n);
n = fread(fid,1,'int32');
parm.daq.HFdaq.HFchannelStr = fgets(fid,n);
parm.daq.HFdaq.channels = str2num(parm.daq.HFdaq.HFchannelStr);
parm.daq.HFdaq.Nchannels = length(parm.daq.HFdaq.channels);
parm.daq.HFdaq.fs_MHz = fread(fid,1,'float64');       % HFFS<DBL>:  30
parm.daq.HFdaq.pts = fread(fid,1,'int32');           % HFSamples<I32>:  2048

% PulseRate<DBL>:  2000
% HFDigitGain<DBL>:  1
% HFVertRng<DBL>:  1
% HFTrigDelay<DBL>:  0

vals = fread(fid,[1,4],'float64');
parm.daq.HFdaq.pulseRepRate_Hz = vals(1);
parm.daq.HFdaq.Gain = vals(2);
parm.daq.HFdaq.vertRange = vals(3);
parm.daq.HFdaq.vertRangeAdj = parm.daq.HFdaq.vertRange;
parm.daq.HFdaq.TrigDelay_us = vals(4);

% LFChan<String>:  17
parm.daq.LFdaq.LFchannelStr = fgets(fid,fread(fid,1,'int32'));
parm.daq.LFdaq.channels = str2num(parm.daq.LFdaq.LFchannelStr);

% LFFS<DBL>:  10
% LFVertRng<DBL>:  5
% Duration<DBL>:  30
vals = fread(fid,[1,3],'float64');
parm.daq.LFdaq.fs_KHz = vals(1);
parm.daq.LFdaq.fs_Hz = parm.daq.LFdaq.fs_KHz * 1000;
parm.daq.LFdaq.LFVertRng = vals(2);
parm.daq.HFdaq.duration_ms = vals(3);
parm.daq.LFdaq.duration_s = vals(3)/1000;

parm.daq.LFdaq.pts = parm.daq.LFdaq.fs_KHz * parm.daq.LFdaq.duration_s*1000;

parm.daq.NBurstAvg = fread(fid,1,'int32'); % NumOfAvg<I32>:  50
parm.daq.LFdaq.NBurstAvg = parm.daq.NBurstAvg;
parm.daq.burstRepRate_Hz = parm.daq.HFdaq.pulseRepRate_Hz;
parm.daq.PreTriggerDelay_ms = 0;
parm.daq.duration_ms = parm.daq.HFdaq.duration_ms;

parm.daq.BurstRate = fread(fid,1,'float64');% BurstRate<DBL>:  3

parm.daq.HFdaq.Coupling = fread(fid,1,'int32'); % HFCoupling<I32>:  0
parm.daq.NumAcq = fread(fid,1,'int32'); % NumAcq<I32>:  21
parm.daq.HFdaq.NoBurstTriggers = round(parm.daq.duration_ms * parm.daq.burstRepRate_Hz/1000);

% #####PacingSignal#####
% IPAddress<String>:  192.168.1.2
% Shape<EW>: Sin
% Freq (Hz)<DBL>:  200
% Ampl (V)<DBL>:  5
% Load Imp<DBL>:  50
% BurstMode<Boolean>:  1
% TrigSrc<EW>: IMM
% Cycles<I32>:  3
% StartPha<DBL>:  0
% Period (ms)<DBL>:  100
% Enable<Boolean>:  1
% #####END OF PacingSignal#####

nval = fread(fid,1,'uint16');
cur_src = {'Agilent','NIDAQmx'};

parm.hp.curSource = cur_src{nval+1};
shape_idx = fread(fid,1,'uint16');
PulseShape = {'Sin','Square','Pulse','ECG','ARB'};
parm.hp.shape = PulseShape{shape_idx+1};

vals = fread(fid,[1,3],'float64');
parm.hp.freq = vals(1);
parm.hp.volt = vals(2);
parm.hp.imp = vals(3);

parm.hp.burstmode = logical(fread(fid,1,'char'));

trig_src = {'IMM','EXT'};
trig_idx = fread(fid,1,'uint16');

parm.hp.trig = trig_src{trig_idx+1};
parm.hp.cycles = fread(fid,1,'int32');

parm.hp.posphase = fread(fid,1,'double');
parm.hp.BurstPeriod = fread(fid,1,'double');
parm.hp.PacingEnable = logical(fread(fid,1,'char'));
parm.hp.period = [];

%% velmex
parm.velmex.XDist = fread(fid,1,'double');

if isempty(parm.velmex.XDist)
    parm.velmex.XDist = 1;
    parm.velmex.XNStep = parm.daq.NumAcq;
    parm.velmex.FastAxis = 'X';
    parm.velmex.YDist = 0;
    parm.velmex.YNStep = 1;
    parm.velmex.SlowAxis = 'Y';
else
    
    parm.velmex.XNStep = fread(fid,1,'int32');
    if parm.velmex.XNStep < 1 || parm.velmex.XNStep > 1e4
        fclose(fid);
        error('parm.velmex.XNStep is invalid.');
    end
    
    
    val = fread(fid,1,'int16');
    
    parm.velmex.FastAxis = char('X'+val);
    
    parm.velmex.YDist = fread(fid,1,'double');
    parm.velmex.YNStep = fread(fid,1,'int32');
    parm.velmex.SlowAxis = char('X'+fread(fid,1,'int16'));
end

parm.velmex.ZDist = 0;
parm.velmex.ZNStep = 1;
parm.YScan = false;


% parm.velmex = struct('FastAxis','X','SlowAxis','Y',...
%     'XDist',1,'XNStep',parm.daq.NumAcq,...
%     'YDist',0,'YNStep',1,...
%     'ZDist',0,'ZNStep',1,...
%     'YScan',false);

% parm.sp = struct('AEChanSel', 1,...
%     'AEGain', 180,...
%     'AEOpt', 1,...
%     'AE_tmRng', [17 37],...
%     'AEdBRng', [-15 0],...
%     'AEopt', 1,...
%     'Colormap', 'hot',...
%     'DemodFreq', 2.2000,...
%     'FilterType', 'HannWin',...
%     'FrameRate', 15,...
%     'FrameSelect', 1,...
%     'HFFiltWin', [1.2000 2 2.5000 3.3000],...
%     'HFHannWinOn', 1,...
%     'HFMatchFiltOn', 0,...
%     'ImageAutoScale', 0,...
%     'ImageScale', 1,...
%     'LFFiltWin', [0.1000 0.1500 0.2500 0.3000],...
%     'LFGain', 2,...
%     'LFHannWinOn', 1,...
%     'LFImp', 1000,...
%     'LFMatchFiltOn', 0,...
%     'LF_tmRng', [0 44],...
%     'MedFiltOn', 0,...
%     'MedFiltSize', [],...
%     'PE_tmRng', [34 74],...
%     'PEdBRng', [-30 0],...
%     'PlaneChoice', 1,...
%     'PulsComp', 0,...
%     'SaveFigs', 0,...
%     'ScanPt', 14,...
%     'ShowEnv', 0,...
%     'SlicePos', 1,...
%     'SndSpeed', 1.4800,...
%     'XPos', 14,...
%     'YPos', 1,...
%     'bSlowTimeMatch', 0,...
%     'bmFrames', 1,...
%     'nsRng', [50 70],...
%     'szAverageRange', [],...
%     'xdcr', 'P4-2V',...
%     'yAxRev', 0,...
%     'zRange', [25 55]);

end

function [data, parm]  = read_info_v170427(fid,parm)

data = [];
% parm = struct;

parm.daq.SingleDev = true;
n = fread(fid,1,'int32');
parm.daq.HFdaq.DeviceName = fgets(fid,n);
n = fread(fid,1,'int32');
parm.daq.HFdaq.HFchannelStr = fgets(fid,n);
parm.daq.HFdaq.channels = str2num(parm.daq.HFdaq.HFchannelStr);
parm.daq.HFdaq.Nchannels = length(parm.daq.HFdaq.channels);
parm.daq.HFdaq.fs_MHz = fread(fid,1,'float64');       % HFFS<DBL>:  30
parm.daq.HFdaq.pts = fread(fid,1,'int32');           % HFSamples<I32>:  2048

% PulseRate<DBL>:  2000
% HFDigitGain<DBL>:  1
% HFVertRng<DBL>:  1
% HFTrigDelay<DBL>:  0

vals = fread(fid,[1,4],'float64');
parm.daq.HFdaq.pulseRepRate_Hz = vals(1);
parm.daq.HFdaq.Gain = vals(2);
parm.daq.HFdaq.vertRange = vals(3);
parm.daq.HFdaq.vertRangeAdj = parm.daq.HFdaq.vertRange;
parm.daq.HFdaq.TrigDelay_us = vals(4);

% LFChan<String>:  17
parm.daq.LFdaq.LFchannelStr = fgets(fid,fread(fid,1,'int32'));
parm.daq.LFdaq.channels = str2num(parm.daq.LFdaq.LFchannelStr);

% LFFS<DBL>:  10
% LFVertRng<DBL>:  5
% Duration<DBL>:  30
vals = fread(fid,[1,3],'float64');
parm.daq.LFdaq.fs_KHz = vals(1);
parm.daq.LFdaq.fs_Hz = parm.daq.LFdaq.fs_KHz * 1000;
parm.daq.LFdaq.LFVertRng = vals(2);
parm.daq.HFdaq.duration_ms = vals(3);
parm.daq.LFdaq.duration_s = vals(3)/1000;

parm.daq.LFdaq.pts = parm.daq.LFdaq.fs_KHz * parm.daq.LFdaq.duration_s*1000;

parm.daq.NBurstAvg = fread(fid,1,'int32'); % NumOfAvg<I32>:  50
parm.daq.LFdaq.NBurstAvg = parm.daq.NBurstAvg;
parm.daq.burstRepRate_Hz = parm.daq.HFdaq.pulseRepRate_Hz;
parm.daq.PreTriggerDelay_ms = 0;
parm.daq.duration_ms = parm.daq.HFdaq.duration_ms;

% parm.daq.BurstRate = fread(fid,1,'float64');% BurstRate<DBL>:  3

parm.daq.HFdaq.Coupling = fread(fid,1,'int32'); % HFCoupling<I32>:  0
% parm.daq.NumAcq = fread(fid,1,'int32'); % NumAcq<I32>:  21
parm.daq.HFdaq.NoBurstTriggers = round(parm.daq.duration_ms * parm.daq.burstRepRate_Hz/1000);
parm.daq.burstMode = logical(fread(fid,1,'char'));
% 
% #####ScanParam#####
% FastScanSize<DBL>:  0
% FastScanPt<I32>:  1
% Fast Axis<EW>: X- Axis
% SlowScanSize<DBL>:  0
% SlowScanPt<I32>:  1
% Slow Axis<EW>: Y- Axis
% Unit<EW>: Millimeter
% #####END OF ScanParam#####

parm.velmex.XDist = fread(fid,1,'double');
parm.velmex.XNStep = fread(fid,1,'int32');
parm.velmex.FastAxis = char('X' + int8(fread(fid,1,'uint16')));
parm.velmex.YDist = fread(fid,1,'double');
parm.velmex.YNStep = fread(fid,1,'int32');
parm.velmex.SlowAxis = char('X' + int8(fread(fid,1,'uint16')));

velmex_units = {'mm','in'};
parm.velmex.unit = velmex_units{fread(fid,1,'uint16')+1};

parm.velmex.ZDist = 0;
parm.velmex.ZNStep = 1;

% #####PacingSignal#####
% IPAddress<String>:  192.168.1.2
% Shape<EW>: Sin
% Freq (Hz)<DBL>:  200
% Ampl (V)<DBL>:  5
% Load Imp<DBL>:  50
% BurstMode<Boolean>:  1
% TrigSrc<EW>: IMM
% Cycles<I32>:  3
% StartPha<DBL>:  0
% Period (ms)<DBL>:  100
% Enable<Boolean>:  1
% #####END OF PacingSignal#####

% nval = fread(fid,1,'uint16');
% cur_src = {'Agilent','NIDAQmx'};
% parm.hp.curSource = cur_src{nval+1};

param.IPAddr = fgets(fid,fread(fid,1,'int32'));


shape_idx = fread(fid,1,'uint16');
PulseShape = {'Sin','Square','Pulse','ECG','ARB'};
parm.hp.shape = PulseShape{shape_idx+1};

vals = fread(fid,[1,3],'float64');
parm.hp.freq = vals(1);
parm.hp.volt = vals(2);
parm.hp.imp = vals(3);

parm.hp.burstmode = logical(fread(fid,1,'char'));

trig_src = {'IMM','EXT'};
trig_idx = fread(fid,1,'uint16');

parm.hp.trig = trig_src{trig_idx+1};
parm.hp.cycles = fread(fid,1,'int32');

parm.hp.posphase = fread(fid,1,'double');
parm.hp.BurstPeriod = fread(fid,1,'double');
% parm.hp.PacingEnable = logical(fread(fid,1,'char'));
% parm.hp.period = [];
parm.YScan = false;


% parm.velmex = struct('FastAxis','X','SlowAxis','Y',...
%     'XDist',1,'XNStep',parm.daq.NumAcq,...
%     'YDist',0,'YNStep',1,...
%     'ZDist',0,'ZNStep',1,...
%     'YScan',false);

% parm.sp = struct('AEChanSel', 1,...
%     'AEGain', 180,...
%     'AEOpt', 1,...
%     'AE_tmRng', [17 37],...
%     'AEdBRng', [-15 0],...
%     'AEopt', 1,...
%     'Colormap', 'hot',...
%     'DemodFreq', 2.2000,...
%     'FilterType', 'HannWin',...
%     'FrameRate', 15,...
%     'FrameSelect', 1,...
%     'HFFiltWin', [1.2000 2 2.5000 3.3000],...
%     'HFHannWinOn', 1,...
%     'HFMatchFiltOn', 0,...
%     'ImageAutoScale', 0,...
%     'ImageScale', 1,...
%     'LFFiltWin', [0.1000 0.1500 0.2500 0.3000],...
%     'LFGain', 2,...
%     'LFHannWinOn', 1,...
%     'LFImp', 1000,...
%     'LFMatchFiltOn', 0,...
%     'LF_tmRng', [0 44],...
%     'MedFiltOn', 0,...
%     'MedFiltSize', [],...
%     'PE_tmRng', [34 74],...
%     'PEdBRng', [-30 0],...
%     'PlaneChoice', 1,...
%     'PulsComp', 0,...
%     'SaveFigs', 0,...
%     'ScanPt', 14,...
%     'ShowEnv', 0,...
%     'SlicePos', 1,...
%     'SndSpeed', 1.4800,...
%     'XPos', 14,...
%     'YPos', 1,...
%     'bSlowTimeMatch', 0,...
%     'bmFrames', 1,...
%     'nsRng', [50 70],...
%     'szAverageRange', [],...
%     'xdcr', 'P4-2V',...
%     'yAxRev', 0,...
%     'zRange', [25 55]);

end


