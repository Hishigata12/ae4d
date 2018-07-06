function [param, HFData, LFData] = read_ucsdi_data(InfoFile, ScanPt)
% function [param, AEData, PEData, LFData] = read_ucsdi_data_v2(InfoFile, ScanPt)
% load ucsdi data recorded using the LabVIEW DAQ standalone program
% Inputs:
%   InfoFile: scan & data acquisition configuration file (*_info.dat)
%   ScanPt:   For single point M-mode, ScanPt = 1; for B-mode or C-mode
%           scan, ScanPt is the position index from 1 to the number of
%           steps.
%           If this argument is omitted, this program will read the data of
%           the first scan point
% Outputs:
%   param:  the structure of scan & daq parameters.
%   AEData: A multi-channel M-Mode UCSDI data 
%           [3-D, fast-time (column) x slow-time(row) x channels (page)]
%   PEData: Pulse echo data, [2-D, fast-time (column) x slow-time(row)]
%   LFData: Low frequency data

HFData = [];
LFData = [];
param = struct;

narginchk(1,2);

if(~exist(InfoFile,'file'))
    uiwait(errordlg([InfoFile,' does not exist'],mfilename,'modal'));
    return;
end

if(exist('ScanPt','var')==0 || isempty(ScanPt))
    ScanPt = 1;
end

[f_path, f_name,f_ext] = fileparts(InfoFile);
if(isempty(regexp([f_name,f_ext],'_info.dat$', 'once')))
    uiwait(errordlg([InfoFile,' is not a valid UCSDI info file.'],mfilename,'modal'));
    return;    
end

try
    param = read_ucsdi_info(InfoFile);
catch ME
    stk_info = dbstack('-completenames');
    uiwait(errordlg([ME.message, ' in ', stk_info(1).file,...
        ' at line: ', num2str(stk_info(1).line-2)], mfilename,'modal'));
    return;
end

assignin('base','parascan',param);

if(nargout==1)
    return;
end

% param.nScanPt = ScanPt;
f_root = regexprep([f_name,f_ext],'_info.dat$','');

nver = param.version;
switch(nver)
    case {140908,160714}
%     case 160714
        hffile_d1 = [f_root,num2str(ScanPt,'_P%04d'),'_HF_Avg.dat'];
        hffile_d2 = [];
%         lffile = [f_root,num2str(ScanPt,'_P%04d'),'_LF_Avg.dat'];
    case {161108,170103}
        hffile_d1 = [f_root,num2str(ScanPt,'_P%04d'),'_HF_Avg.dat'];
        hffile_d2 = [];
%         lffile = [f_root,num2str(ScanPt,'_P%04d'),'_LF_Avg.dat'];
    case 170427
        hffile_d1 = [f_root,'_HF_Avg.dat'];
        hffile_d2 = [];
    otherwise
        hffile_d1 = [f_root,num2str(ScanPt,'_P%04d'),'_D1_Avg.bin'];
        hffile_d2 = regexprep(hffile_d1,'_D1_Avg.bin$','_D2_Avg.bin');
%         lffile = regexprep(hffile_d1,'_D1_Avg.bin$','_LF_Avg.bin');
end


%%

[HFData] = read_hfdata(fullfile(f_path,hffile_d1),ScanPt);

% param.daq = daqparam;

if(~isempty(hffile_d2))
    hfdata2 = read_hfdata(fullfile(f_path,hffile_d2),ScanPt);
    
    dsize1 = size(HFData);
    chList1 = (1:dsize1(end))*2-1;

    dsize2 = size(hfdata2);
    chList2 = (1:dsize2(end))*2;
    chList = [chList1,chList2];
    
    HFData(:,:,dsize1(end)+(1:dsize2(end))) = hfdata2;
    HFData = HFData(:,:,chList);
    
    clear hfdata2;
end

HFData = HFData/param.daq.HFdaq.Gain;


%% read low frequency data
if nargout >2
    [LFData] = read_lfdata(f_root,ScanPt);
end

end


%% read high frequency data
function [hfdata,param] = read_hfdata(hffile,ScanPt)

blk_idx = ScanPt;

if exist(hffile,'file')
    fid = fopen(hffile,'rb');
    blk_idx = 1;
else
    hf_avg_file = regexprep(hffile,'_P[0-9]{4,4}_','_');
    fid = fopen(hf_avg_file,'rb');
end

if fid < 0
    error('Cannot open HF data file');
end

sztype = fgets(fid,fread(fid,1,'int32'));
if(strcmpi(sztype(1:5),'ucsdi') == 0)
    fclose(fid);
    return;
end

nver = fread(fid,1,'int32');
switch(nver)
    case {140908,160714}
        [hfdata,param] = read_hfdata_v140908(fid);
    case {161108,170103}
        [hfdata,param] = read_hfdata_v161108(fid,blk_idx);
    case 170427
        [hfdata,param] = read_hfdata_v170427(fid,ScanPt);
    otherwise
        [hfdata,param] = read_hfdata_default(fid);
end
param.Version = nver;
fclose(fid);
end
%%

function [hfdata,param] = read_hfdata_v140908(fid)

n = fread(fid,1,'int32');
HFdaq.DeviceName = fgets(fid,n);


HFdaq.ChanList = fgets(fid,fread(fid,1,'int32'));
HFdaq.fs_MHz = fread(fid,1,'float64');
HFdaq.nPts = fread(fid,1,'int32');
vdata = fread(fid,4,'float64');
% HFdaq.TrigDelay_us = vdata(1)*1e6;
% HFdaq.Gain = vdata(2);
% HFdaq.VertRange = vdata(3);
% HFdaq.PulseRepRate = vdata(4);
HFdaq.pulseRepRate_Hz = vdata(1);
HFdaq.Gain = vdata(2);
HFdaq.vertRange = vdata(3);
HFdaq.vertRangeAdj = vdata(3);
HFdaq.TrigDelay_us = vdata(4);


LFdaq.ChanList = fgets(fid,fread(fid,1,'int32'));
vdata = fread(fid,3,'float64');
LFdaq.fs_KHz = vdata(1);
LFdaq.fs_Hz = vdata(1)*1000;
LFdaq.VertRange = vdata(2);
LFdaq.vertRangeAdj = vdata(2);
LFdaq.duration_ms = vdata(3);
LFdaq.duration_s = vdata(3)/1000;

% LFdaq.BurstRate_Hz = vdata(4);

NAvgNum = fread(fid,1,'int32');
LFdaq.burstRepRate_Hz = fread(fid,1,'float64');
% LFdaq.TrigDelay_ms = fread(fid,1,'float64')*1e3;
HFdaq.Coupling = fread(fid,1,'int32');

HFdaq.NChannels = length(str2num(HFdaq.ChanList)); 
% vdata = fread(fid,5,'int32');
% HFdaq.NumOfChans = vdata(1);
% HFdaq.NumOfFrames = vdata(2);
% LFdaq.NumOfChans = vdata(3);
% LFdaq.nPts = vdata(4);
% 
% FileRoot = fgets(fid,fread(fid,1,'int32'));
% IsAvgOnly = logical(fread(fid,1,'int8'));



n = fread(fid,1,'int32');

dsize = fread(fid,[1,n],'int32');
data = single(fread(fid,[dsize(end),prod(dsize(1:end-1))],'single'));
hfdata = permute(reshape(data,fliplr(dsize)),[1,3,2]);
% figure(1);clf;
% imagesc(data);
% axis xy;

% figure(2);clf;
% plot(mean(data,1));
% 
% fclose(fid);

param.HFdaq = HFdaq;
param.LFdaq = LFdaq;
param.NAvgNum = NAvgNum;
% param.SaveAvgOnly = IsAvgOnly;
% param.FileRoot = FileRoot;
% param.Version = nver;
end

%%
function [hfdata,param] = read_hfdata_default(fid)

HFdaq.ChanList = fgets(fid,fread(fid,1,'int32'));
HFdaq.fs_MHz = fread(fid,1,'float64')/1e6;
HFdaq.nPts = fread(fid,1,'int32');
vdata = fread(fid,4,'float64');
HFdaq.TrigDelay_us = vdata(1)*1e6;
HFdaq.Gain = vdata(2);
HFdaq.VertRange = vdata(3);
HFdaq.PulseRepRate = vdata(4);

LFdaq.ChanList = fgets(fid,fread(fid,1,'int32'));
vdata = fread(fid,4,'float64');
LFdaq.fs_kHz = vdata(1)/1e3;
LFdaq.VertRange = vdata(2);
LFdaq.Durattion_ms = vdata(3)*1e3;
LFdaq.BurstRate_Hz = vdata(4);

NAvgNum = fread(fid,1,'int32');
LFdaq.TrigDelay_ms = fread(fid,1,'float64')*1e3;

vdata = fread(fid,5,'int32');
HFdaq.NumOfChans = vdata(1);
HFdaq.NumOfFrames = vdata(2);
LFdaq.NumOfChans = vdata(3);
LFdaq.nPts = vdata(4);

FileRoot = fgets(fid,fread(fid,1,'int32'));
IsAvgOnly = logical(fread(fid,1,'int8'));

n = fread(fid,1,'int32');

dsize = fread(fid,n,'int32');
data = single(fread(fid,[dsize(end),prod(dsize(1:end-1))],'single'));
hfdata = permute(reshape(data,[HFdaq.nPts,HFdaq.NumOfChans,HFdaq.NumOfFrames]),[1,3,2]);
% figure(1);clf;
% imagesc(data);
% axis xy;

% figure(2);clf;
% plot(mean(data,1));
% 
% fclose(fid);

param.HFdaq = HFdaq;
param.LFdaq = LFdaq;
param.NAvgNum = NAvgNum;
param.SaveAvgOnly = IsAvgOnly;
param.FileRoot = FileRoot;
% param.Version = nver;
end

function [data,parm] = read_hfdata_v161108(fid,ScanPt)

% nPos = ftell(fid);
n = fread(fid,1,'int32');
dsize = fread(fid,[1,n],'int32');

% blk_size = (length(dsize) + 1)*4 + prod(dsize)*4;
blk_size = prod(dsize)*4;
fseek(fid,blk_size*(ScanPt-1),'cof');
data = fread(fid,[prod(dsize),1],'single');
tsize = dsize(1)*dsize(2)*dsize(3);
if tsize ~= length(data)
    data = zeros(tsize,1);
    disp(['Missing point ' num2str(ScanPt)]);
end
data = reshape(data,fliplr(dsize));

data = permute(data,[1,3,2]);

parm.dsize = dsize;
parm.ScanPt = ScanPt;

% parm.daq.SingleDev = true;
% n = fread(fid,1,'int32');
% parm.daq.HFdaq.DeviceName = fgets(fid,n);
% n = fread(fid,1,'int32');
% parm.daq.HFdaq.HFchannelStr = fgets(fid,n);
% parm.daq.HFdaq.channels = str2num(parm.daq.HFdaq.HFchannelStr);
% parm.daq.HFdaq.Nchannels = length(parm.daq.HFdaq.channels);
% parm.daq.HFdaq.fs_MHz = fread(fid,1,'float64');       % HFFS<DBL>:  30
% parm.daq.HFdaq.pts = fread(fid,1,'int32');           % HFSamples<I32>:  2048
% 
% % PulseRate<DBL>:  2000
% % HFDigitGain<DBL>:  1
% % HFVertRng<DBL>:  1
% % HFTrigDelay<DBL>:  0
% 
% vals = fread(fid,[1,4],'float64');
% parm.daq.HFdaq.pulseRepRate_Hz = vals(1);
% parm.daq.HFdaq.Gain = vals(2);
% parm.daq.HFdaq.vertRange = vals(3);
% parm.daq.HFdaq.vertRangeAdj = parm.daq.HFdaq.vertRange;
% parm.daq.HFdaq.TrigDelay_us = vals(4);
% 
% % LFChan<String>:  17
% parm.daq.LFdaq.LFchannelStr = fgets(fid,fread(fid,1,'int32'));
% parm.daq.LFdaq.channels = str2num(parm.daq.LFdaq.LFchannelStr);
% 
% % LFFS<DBL>:  10
% % LFVertRng<DBL>:  5
% % Duration<DBL>:  30
% vals = fread(fid,[1,3],'float64');
% parm.daq.LFdaq.fs_KHz = vals(1);
% parm.daq.LFdaq.LFVertRng = vals(2);
% parm.daq.HFdaq.duration_ms = vals(3);
% parm.daq.LFdaq.duration_s = vals(3)/1000;
% 
% parm.daq.LFdaq.pts = parm.daq.LFdaq.fs_KHz * parm.daq.LFdaq.duration_s*1000;
% 
% parm.daq.NBurstAvg = fread(fid,1,'int32'); % NumOfAvg<I32>:  50
% parm.daq.LFdaq.NBurstAvg = parm.daq.NBurstAvg;
% parm.daq.burstRepRate_Hz = parm.daq.HFdaq.pulseRepRate_Hz;
% parm.daq.PreTriggerDelay_ms = 0;
% parm.daq.duration_ms = parm.daq.HFdaq.duration_ms;
% 
% parm.BurstRate = fread(fid,1,'float64');% BurstRate<DBL>:  3
% 
% parm.daq.HFdaq.Coupling = fread(fid,1,'int32'); % HFCoupling<I32>:  0
% parm.NumAcq = fread(fid,1,'int32'); % NumAcq<I32>:  21
% parm.daq.HFdaq.NoBurstTriggers = round(parm.daq.duration_ms * parm.daq.burstRepRate_Hz/1000);
% 
% % #####PacingSignal#####
% % IPAddress<String>:  192.168.1.2
% % Shape<EW>: Sin
% % Freq (Hz)<DBL>:  200
% % Ampl (V)<DBL>:  5
% % Load Imp<DBL>:  50
% % BurstMode<Boolean>:  1
% % TrigSrc<EW>: IMM
% % Cycles<I32>:  3
% % StartPha<DBL>:  0
% % Period (ms)<DBL>:  100
% % Enable<Boolean>:  1
% % #####END OF PacingSignal#####
% 
% parm.hp.TCPIPaddress = fgets(fid,fread(fid,1,'int32'));
% shape_idx = fread(fid,1,'uint16');
% PulseShape = {'Sin','Square','Pulse','ECG','ARB'};
% parm.hp.shape = PulseShape{shape_idx+1};
% 
% vals = fread(fid,[1,3],'float64');
% parm.hp.freq = vals(1);
% parm.hp.volt = vals(2);
% parm.hp.imp = vals(3);
% 
% parm.hp.burstmode = logical(fread(fid,1,'char'));
% 
% trig_src = {'IMM','EXT'};
% trig_idx = fread(fid,1,'uint16');
% 
% parm.hp.trig = trig_src{trig_idx+1};
% parm.hp.cycles = fread(fid,1,'int32');
% 
% parm.hp.posphase = fread(fid,1,'double');
% parm.BurstPeriod = fread(fid,1,'double');
% parm.PacingEnable = logical(fread(fid,1,'char'));
% parm.hp.period = [];
% 
% parm.headerlen = ftell(fid);
% 
% n = fread(fid,1,'int32');
% dsize = fread(fid,[1,n],'int32');
% 
% blk_size = (length(dsize)+1)*4 + prod(dsize)*4;
% 
% cur_pos = ftell(fid);
% fseek(fid,0,'eof');
% flen = ftell(fid);
% fseek(fid,cur_pos,'bof');
% 
% NumAcq = floor((flen - parm.headerlen)/blk_size);
% 
% data = zeros(prod(dsize),NumAcq);
% for ni = 1:NumAcq
%     if ni > 1
%         fseek(fid,(numel(dsize)+1)*4,'cof');
%     end
%     data(:,ni) = fread(fid,[prod(dsize),1],'single');
% end
% 
% data = reshape(data,[fliplr(dsize),NumAcq]);
% if ndims(data) == 3
%     data = permute(data,[1,3,2]);
% elseif ndims(data) == 4
%     data = permute(data,[1,3,2,4]);
% end
% parm.fast_t = (1:size(data,1))/parm.daq.HFdaq.fs_MHz + parm.daq.HFdaq.TrigDelay_us;
% daqParm = parm.daq;
end

function [data,parm] = read_hfdata_v170427(fid,ScanPt)

% nPos = ftell(fid);
n = fread(fid,1,'int32');
dsize = fread(fid,[1,n],'int32');

% blk_size = (length(dsize) + 1)*4 + prod(dsize)*4;
blk_size = prod(dsize)*4;
fseek(fid,blk_size*(ScanPt-1),'cof');

data = fread(fid,[prod(dsize),1],'single');
data = reshape(data,fliplr(dsize));

data = permute(data,[1,3,2]);

parm.dsize = dsize;
parm.ScanPt = ScanPt;

% parm.daq.SingleDev = true;
% n = fread(fid,1,'int32');
% parm.daq.HFdaq.DeviceName = fgets(fid,n);
% n = fread(fid,1,'int32');
% parm.daq.HFdaq.HFchannelStr = fgets(fid,n);
% parm.daq.HFdaq.channels = str2num(parm.daq.HFdaq.HFchannelStr);
% parm.daq.HFdaq.Nchannels = length(parm.daq.HFdaq.channels);
% parm.daq.HFdaq.fs_MHz = fread(fid,1,'float64');       % HFFS<DBL>:  30
% parm.daq.HFdaq.pts = fread(fid,1,'int32');           % HFSamples<I32>:  2048
% 
% % PulseRate<DBL>:  2000
% % HFDigitGain<DBL>:  1
% % HFVertRng<DBL>:  1
% % HFTrigDelay<DBL>:  0
% 
% vals = fread(fid,[1,4],'float64');
% parm.daq.HFdaq.pulseRepRate_Hz = vals(1);
% parm.daq.HFdaq.Gain = vals(2);
% parm.daq.HFdaq.vertRange = vals(3);
% parm.daq.HFdaq.vertRangeAdj = parm.daq.HFdaq.vertRange;
% parm.daq.HFdaq.TrigDelay_us = vals(4);
% 
% % LFChan<String>:  17
% parm.daq.LFdaq.LFchannelStr = fgets(fid,fread(fid,1,'int32'));
% parm.daq.LFdaq.channels = str2num(parm.daq.LFdaq.LFchannelStr);
% 
% % LFFS<DBL>:  10
% % LFVertRng<DBL>:  5
% % Duration<DBL>:  30
% vals = fread(fid,[1,3],'float64');
% parm.daq.LFdaq.fs_KHz = vals(1);
% parm.daq.LFdaq.LFVertRng = vals(2);
% parm.daq.HFdaq.duration_ms = vals(3);
% parm.daq.LFdaq.duration_s = vals(3)/1000;
% 
% parm.daq.LFdaq.pts = parm.daq.LFdaq.fs_KHz * parm.daq.LFdaq.duration_s*1000;
% 
% parm.daq.NBurstAvg = fread(fid,1,'int32'); % NumOfAvg<I32>:  50
% parm.daq.LFdaq.NBurstAvg = parm.daq.NBurstAvg;
% parm.daq.burstRepRate_Hz = parm.daq.HFdaq.pulseRepRate_Hz;
% parm.daq.PreTriggerDelay_ms = 0;
% parm.daq.duration_ms = parm.daq.HFdaq.duration_ms;
% 
% parm.BurstRate = fread(fid,1,'float64');% BurstRate<DBL>:  3
% 
% parm.daq.HFdaq.Coupling = fread(fid,1,'int32'); % HFCoupling<I32>:  0
% parm.NumAcq = fread(fid,1,'int32'); % NumAcq<I32>:  21
% parm.daq.HFdaq.NoBurstTriggers = round(parm.daq.duration_ms * parm.daq.burstRepRate_Hz/1000);
% 
% % #####PacingSignal#####
% % IPAddress<String>:  192.168.1.2
% % Shape<EW>: Sin
% % Freq (Hz)<DBL>:  200
% % Ampl (V)<DBL>:  5
% % Load Imp<DBL>:  50
% % BurstMode<Boolean>:  1
% % TrigSrc<EW>: IMM
% % Cycles<I32>:  3
% % StartPha<DBL>:  0
% % Period (ms)<DBL>:  100
% % Enable<Boolean>:  1
% % #####END OF PacingSignal#####
% 
% parm.hp.TCPIPaddress = fgets(fid,fread(fid,1,'int32'));
% shape_idx = fread(fid,1,'uint16');
% PulseShape = {'Sin','Square','Pulse','ECG','ARB'};
% parm.hp.shape = PulseShape{shape_idx+1};
% 
% vals = fread(fid,[1,3],'float64');
% parm.hp.freq = vals(1);
% parm.hp.volt = vals(2);
% parm.hp.imp = vals(3);
% 
% parm.hp.burstmode = logical(fread(fid,1,'char'));
% 
% trig_src = {'IMM','EXT'};
% trig_idx = fread(fid,1,'uint16');
% 
% parm.hp.trig = trig_src{trig_idx+1};
% parm.hp.cycles = fread(fid,1,'int32');
% 
% parm.hp.posphase = fread(fid,1,'double');
% parm.BurstPeriod = fread(fid,1,'double');
% parm.PacingEnable = logical(fread(fid,1,'char'));
% parm.hp.period = [];
% 
% parm.headerlen = ftell(fid);
% 
% n = fread(fid,1,'int32');
% dsize = fread(fid,[1,n],'int32');
% 
% blk_size = (length(dsize)+1)*4 + prod(dsize)*4;
% 
% cur_pos = ftell(fid);
% fseek(fid,0,'eof');
% flen = ftell(fid);
% fseek(fid,cur_pos,'bof');
% 
% NumAcq = floor((flen - parm.headerlen)/blk_size);
% 
% data = zeros(prod(dsize),NumAcq);
% for ni = 1:NumAcq
%     if ni > 1
%         fseek(fid,(numel(dsize)+1)*4,'cof');
%     end
%     data(:,ni) = fread(fid,[prod(dsize),1],'single');
% end
% 
% data = reshape(data,[fliplr(dsize),NumAcq]);
% if ndims(data) == 3
%     data = permute(data,[1,3,2]);
% elseif ndims(data) == 4
%     data = permute(data,[1,3,2,4]);
% end
% parm.fast_t = (1:size(data,1))/parm.daq.HFdaq.fs_MHz + parm.daq.HFdaq.TrigDelay_us;
% daqParm = parm.daq;
end

