function varargout = AvgViewer(varargin)
% AVGVIEWER MATLAB code for AvgViewer.fig
%      AVGVIEWER, by itself, creates a new AVGVIEWER or raises the existing
%      singleton*.
%
%      H = AVGVIEWER returns the handle to a new AVGVIEWER or the handle to
%      the existing singleton*.
%
%      AVGVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AVGVIEWER.M with the given input arguments.
%
%      AVGVIEWER('Property','Value',...) creates a new AVGVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AvgViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AvgViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AvgViewer

% Last Modified by GUIDE v2.5 17-Sep-2019 16:40:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AvgViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @AvgViewer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before AvgViewer is made visible.
function AvgViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AvgViewer (see VARARGIN)

% Choose default command line output for AvgViewer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AvgViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AvgViewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in selectfile.
function selectfile_Callback(hObject, eventdata, handles)
[f, d] = uigetfile(fullfile(pwd,'*_info.mat'));
fl = fullfile(d,f);
set(handles.filedir,'String',fl);
p = open(fl);
param = p.bScanParm;
set(handles.depthval,'String',floor(param.Scan.Zpt));
set(handles.timeval,'String',param.Scan.Tpt);
set(handles.channelval,'String',numel(str2num(param.Daq.HF.Channels)));
set(handles.avgval,'String',param.Scan.Avg);
set(handles.xval,'String',param.Scan.Xpt);
set(handles.yval,'String',param.Scan.Ypt);
set(handles.Xpt,'String',1);
set(handles.Ypt,'String',1);
hfchans = str2num(param.Daq.HF.Channels);
set(handles.HFchan,'String',num2str(hfchans(1)));
set(handles.AvgSelect,'String',['1:',num2str(param.Scan.Avg)]);
set(handles.lfchans,'String',numel(str2num(param.Daq.LF.Channels)));
lfchans = str2num(param.Daq.LF.Channels);
set(handles.LFchan,'String',num2str(lfchans(1)));
assignin('base','param',param);
% hObject    handle to selectfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in openfilehf.
function openfilehf_Callback(hObject, eventdata, handles)
param = evalin('base','param');
X = str2num(get(handles.Xpt,'String'));
Y = str2num(get(handles.Ypt,'String'));
Xr = str2num(get(handles.xval,'String'));
Yr = str2num(get(handles.yval,'String'));
chan2use = str2num(get(handles.HFchan,'String'));
HFchans = str2num(param.Daq.HF.Channels);
[~,chan] = find(HFchans == chan2use);
Avg = str2num(handles.AvgSelect.String);
ScanPt = (Y-1)*Xr+X;
froot = get(handles.filedir,'String');
froot = froot(1:end-9);
fHF = [froot '_HF.dat'];

%Opens HF data file

blk_idx = ScanPt;

if exist(fHF,'file')
    fid = fopen(fHF,'rb');
    blk_idx = 1;
else
    hf_avg_file = regexprep(fHF,'_P[0-9]{4,4}_','_');
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

% part 2
nver = fread(fid,1,'int32');
n = fread(fid,1,'int32');
dsize = fread(fid,[1,n],'int32');

% blk_size = (length(dsize) + 1)*4 + prod(dsize)*4;
blk_size = prod(dsize)*4;
% fseek(fid,blk_size*(ScanPt-1),'cof');
% data = fread(fid,[prod(dsize),1],'single');

for i = 1:length(Avg)
blk_size2 = prod(dsize(1:3))*4;
fseek(fid,blk_size*((ScanPt-1)+blk_size2*(Avg(i)-1)),'cof');
data = fread(fid,[prod(dsize(1:3)),1],'single');
data = reshape(data,[dsize(3),dsize(2),dsize(1)]);
data = permute(data,[3,1,2]);
data = data(:,:,chan);
HF(:,:,i) = data;
end
assignin('base','HF',HF);

%Filters HF Data
time_lc = str2num(handles.timelow.String);
time_hc = str2num(handles.timehigh.String);
depth_lc = str2num(handles.depthlow.String);
depth_hc = str2num(handles.depthhigh.String);
time_band = get(handles.timeband,'Value');
depth_band = get(handles.depthband,'Value');
time_match = get(handles.timematch,'Value');
depth_match = get(handles.depthmatch,'Value');

if time_band == 1
    time_type = 0;
end

if time_match == 1
    time_type = 1;
end

if depth_band == 1
    depth_type = 0;
end
if depth_match == 1
    depth_type = 1;
end


HF = permute(HF,[3 4 2 1]);
if depth_match == 1 || depth_band == 1
    pathstop = find(froot == '\',2,'last');
    path = froot(1:pathstop(1));
    
    path2 = [path 'PEData\'];
    file2 = uigetfile(fullfile(path2,'*PEParm.mat')); %gets US pulse waveform
    
    PE = open([path2 file2]);
    handles.tc.Value = 0;
    HF = w_ae_filt2(param,HF,PE,depth_type,handles,[depth_lc depth_hc]);
    HF = real(HF);
end

if time_match == 1 || time_band == 1
    LF = evalin('base','LF2');
    param.full_sm = 1;
    param.medfilt = 0;
    [HF, LF] = w_slow_filt2(param,HF,LF,time_type,[time_lc time_hc]);
    HF = real(HF);
end

multiWaitbar('CLOSEALL');

%Plots HF Data
HF = permute(HF,[4 3 1 2]); 
figure(11);
total = numel(str2num(handles.AvgSelect.String));
spn = floor(sqrt(total));

if spn ~= sqrt(total)
    spn2 = spn+1;
else 
    spn2 = spn;
end
for i = 1:total
subplot(spn,spn2,i)
clims = [min(min(HF(:,200:end))) max(max(HF(:,200:end)))];
imagesc((HF(:,200:end,i)'),clims./5)

end



% tsize = dsize(1)*dsize(2)*dsize(3)*dsize(4); %Slow Time, Channels, Fast Time, Averages
% 
% if tsize ~= length(data)
%     data = zeros(tsize,1);
%     disp(['Missing point ' num2str(ScanPt)]);
% end

% data = reshape(data,[dsize(3),dsize(2),dsize(1),dsize(4)]); %FT,Chan,ST,Avg
% data = permute(data,[3,1,4,2]); %ST,FT,AVG,Chan
% data = data(:,:,:,p);

fclose(fid);


% hObject    handle to openfilehf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function Xpt_Callback(hObject, eventdata, handles)
% hObject    handle to Xpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xpt as text
%        str2double(get(hObject,'String')) returns contents of Xpt as a double


% --- Executes during object creation, after setting all properties.
function Xpt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xpt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Ypt_Callback(hObject, eventdata, handles)
% hObject    handle to Ypt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Ypt as text
%        str2double(get(hObject,'String')) returns contents of Ypt as a double


% --- Executes during object creation, after setting all properties.
function Ypt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Ypt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function HFchan_Callback(hObject, eventdata, handles)
% hObject    handle to HFchan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HFchan as text
%        str2double(get(hObject,'String')) returns contents of HFchan as a double


% --- Executes during object creation, after setting all properties.
function HFchan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HFchan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function AvgSelect_Callback(hObject, eventdata, handles)
% hObject    handle to AvgSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AvgSelect as text
%        str2double(get(hObject,'String')) returns contents of AvgSelect as a double


% --- Executes during object creation, after setting all properties.
function AvgSelect_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AvgSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LFchan_Callback(hObject, eventdata, handles)
% hObject    handle to LFchan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LFchan as text
%        str2double(get(hObject,'String')) returns contents of LFchan as a double


% --- Executes during object creation, after setting all properties.
function LFchan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LFchan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in openfilelf.
function openfilelf_Callback(hObject, eventdata, handles)
param = evalin('base','param');
X = str2num(get(handles.Xpt,'String'));
Y = str2num(get(handles.Ypt,'String'));
Xr = str2num(get(handles.xval,'String'));
Yr = str2num(get(handles.yval,'String'));
ScanPt = (Y-1)*Xr+X;
froot = get(handles.filedir,'String');
froot = froot(1:end-9);
fLF = [froot '_LF.dat'];

% Gets LF data
fid = fopen(fLF,'rb');
nBytes = fread(fid,1,'int32');
nPos = ftell(fid);
stype = strtrim(fgets(fid,nBytes)); % LVDATA

nver = fread(fid,1,'int32');
n = fread(fid,1,'int32');
dsize = fread(fid,[1,n],'int32');

cur_pos = ftell(fid);

% first trace
param.nPts = dsize(2);
param.numChan = dsize(1);

% blk_size = (length(dsize) + 1)*4 + prod(dsize)*4;
blk_size = prod(dsize)*4;
offset_n = cur_pos + blk_size * (ScanPt-1);
fseek(fid,offset_n,'bof');

data2 = fread(fid,fliplr(dsize),'single');
dec = param.Scan.Avg;
leng = size(data2,1)/dec;
numchans = dsize(1);
chan2use = str2num(handles.LFchan.String);
LFchans = str2num(param.Daq.LF.Channels);
[~,chan] = find(LFchans == chan2use);
data2 = reshape(data2,[prod(dsize),1]);
for i = 1:dec*numchans
    SingleAvg(:,i) = data2(leng*(i-1)+1:leng*i);
end
for i = 1:numchans
    for j = 1:dec
        data(:,j,i) = SingleAvg(:,i+(j-1)*(numchans));
    end
end

fclose(fid);

%Plots LF Data
figure(10);
Avg = str2num(handles.AvgSelect.String);
total = numel(str2num(handles.AvgSelect.String));
spn = floor(sqrt(total));

if spn ~= sqrt(total)
    spn2 = spn+1;
else 
    spn2 = spn;
end
for i = 1:total
subplot(spn,spn2,i)
plot(data(:,Avg(i),chan))
end
assignin('base','LF2',data(:,1));
    
% hObject    handle to openfilelf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in newpath.
function newpath_Callback(hObject, eventdata, handles)
fpathfull = get(handles.filedir,'String');
stop = find(fpathfull == '\',1,'last');
froot = fpathfull(1:stop);
cd(froot);
% hObject    handle to newpath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function timelow_Callback(hObject, eventdata, handles)
% hObject    handle to timelow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timelow as text
%        str2double(get(hObject,'String')) returns contents of timelow as a double


% --- Executes during object creation, after setting all properties.
function timelow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timelow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function timehigh_Callback(hObject, eventdata, handles)
% hObject    handle to timehigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of timehigh as text
%        str2double(get(hObject,'String')) returns contents of timehigh as a double


% --- Executes during object creation, after setting all properties.
function timehigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timehigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function depthlow_Callback(hObject, eventdata, handles)
% hObject    handle to depthlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of depthlow as text
%        str2double(get(hObject,'String')) returns contents of depthlow as a double


% --- Executes during object creation, after setting all properties.
function depthlow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to depthlow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function depthhigh_Callback(hObject, eventdata, handles)
% hObject    handle to depthhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of depthhigh as text
%        str2double(get(hObject,'String')) returns contents of depthhigh as a double


% --- Executes during object creation, after setting all properties.
function depthhigh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to depthhigh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in timeband.
function timeband_Callback(hObject, eventdata, handles)
% hObject    handle to timeband (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of timeband


% --- Executes on button press in timematch.
function timematch_Callback(hObject, eventdata, handles)
% hObject    handle to timematch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of timematch


% --- Executes on button press in depthmatch.
function depthmatch_Callback(hObject, eventdata, handles)
% hObject    handle to depthmatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of depthmatch


% --- Executes on button press in depthband.
function depthband_Callback(hObject, eventdata, handles)
% hObject    handle to depthband (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of depthband
 