function varargout = ae4d(varargin)
% AE4D MATLAB code for ae4d.fig
%      AE4D, by itself, creates a new AE4D or raises the existing
%      singleton*.
%
%     open  H = AE4D returns the handle to a new AE4D or the handle to
%      the existing singleton*.
%
%      AE4D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AE4D.M with the given input arguments.
%
%      AE4D('Property','Value',...) creates a new AE4D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ae4d_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ae4d_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ae4d

% Last Modified by GUIDE v2.5 12-Oct-2018 13:25:53

% Begin initialization code - DO NOT EDIT

% Search Guide
% @001 = Opening Function
% @002 = Plot 1
% @003 = Load AE Data
% @004 = Create AE Data
% @005 = Movie
% @006 = Enhance
% @007 = Chop
% @008 = Max Ranges
% @009 = Envelope
% @010 = iRadon
% @011 = dB
% @012 = Plot LF
% @013 = Sensitivity
% @014 = FWHM
% @015 = Reset Plots
% @016 = Modify (baseband)
% @017 = 4D Overlay
% @018 = 2D Overlay
% @019 = Get Noise
% @020 = Create PE
% @021 = Load PE
% @022 = Use PE
% @023 = Save Figure
% @024 = Return to PE
% @025 = Plot data in all 4 windows
% @026 = Merge AE datasets into structure
% @027 = Use particular channel for merged data
% @028 = Create LF Movie
% @029 = Plot on Click
% @030 = T shift slider
% @031 = Z shift slider
% @032 = Y shift slider
% @033 = Stitch for fusing two images

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ae4d_OpeningFcn, ...
    'gui_OutputFcn',  @ae4d_OutputFcn, ...
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


% --- Executes just before ae4d is made visible.
function ae4d_OpeningFcn(hObject, eventdata, handles, varargin) % @001
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ae4d (see VARARGIN)

% Choose default command line output for ae4d
handles.output = hObject;
set(handles.xR,'String','0 0');
set(handles.yR,'String','0 0');
set(handles.zR,'String','0 0');
set(handles.tR,'String','0 0');
set(handles.xP,'String','0');
set(handles.yP,'String','0');
set(handles.zP,'String','0');
set(handles.tP,'String','5.5');
%set(handles.aeR, 'String','-9 0');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ae4d wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ae4d_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in plotbox1.
function plotbox1_Callback(hObject, eventdata, handles)
% hObject    handle to plotbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plotbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotbox1


% --- Executes during object creation, after setting all properties.
function plotbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plotbox2.
function plotbox2_Callback(hObject, eventdata, handles)
% hObject    handle to plotbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plotbox2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotbox2


% --- Executes during object creation, after setting all properties.
function plotbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xR_Callback(hObject, eventdata, handles)
% hObject    handle to xR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xR as text
%        str2double(get(hObject,'String')) returns contents of xR as a double


% --- Executes during object creation, after setting all properties.
function xR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yR_Callback(hObject, eventdata, handles)
% hObject    handle to yR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yR as text
%        str2double(get(hObject,'String')) returns contents of yR as a double


% --- Executes during object creation, after setting all properties.
function yR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zR_Callback(hObject, eventdata, handles)
% hObject    handle to zR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zR as text
%        str2double(get(hObject,'String')) returns contents of zR as a double


% --- Executes during object creation, after setting all properties.
function zR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tR_Callback(hObject, eventdata, handles)
% hObject    handle to tR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tR as text
%        str2double(get(hObject,'String')) returns contents of tR as a double


% --- Executes during object creation, after setting all properties.
function tR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_ae.
function plot_ae_Callback(hObject, eventdata, handles) % @002
% hObject    handle to plot_ae (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

param = evalin('base','param');
if handles.use_chop.Value == 1
    Xfilt = evalin('base','X_c');
    ax = evalin('base','ax_c');
else
    Xfilt = evalin('base','Xfilt');
    ax = evalin('base','ax');
end

Xfilt = real(Xfilt);
xP = str2double(handles.xP.String);
yP = str2double(handles.yP.String);
zP = str2double(handles.zP.String);
tP = str2double(handles.tP.String);

xR = str2num(handles.xR.String);
if length(xR) == 1
    xR = [xR xR xR];
else
    xR = [xR(1) xR(2) xP];
end
yR = str2num(handles.yR.String);
if length(yR) == 1
    yR = [yR yR yR];
else
    yR = [yR(1) yR(2) yP];
end
zR = str2num(handles.zR.String);
if length(zR) == 1
    zR = [zR zR zR];
else
    zR = [zR(1) zR(2) zP];
end
tR = str2num(handles.tR.String);
if length(tR) == 1
    tR = [tR tR tR];
else
    tR = [tR(1) tR(2) tP];
end
aeR = str2num(handles.aeR.String);
dims = size(Xfilt);
if length(dims) < 3
    dims(3) = 1;
end
%[~,ax] = make_axes(param,dims,[1 2],1);
q.x = 1:dims(1);
q.y = 1:dims(2);
q.z = 1:dims(3);

if handles.plotbox1.Value == 1
    zR(:) = zR(3);
    yR(:) = yR(3);
end

if handles.plotbox1.Value == 2
    tR(:) = tR(3);
    zR(:) = zR(3);
end

if handles.plotbox1.Value == 3
    tR(:) = tR(3);
    yR(:) = yR(3);
end

if handles.plotbox1.Value == 4
    xR(:) = xR(3);
    zR(:) = zR(3);
end

if handles.plotbox1.Value == 5
    tR(:) = tR(3);
    xR(:) = xR(3);
end

if handles.plotbox1.Value == 6
    xR(:) = xR(3);
    yR(:) = yR(3);
end


xInd = find(ax.x >= xR(1)):find(ax.x >= xR(2));
yInd = find(ax.y >= yR(1)):find(ax.y >= yR(2));
zInd = find(ax.depth >= zR(1)):find(ax.depth >= zR(2));

if length(size(Xfilt)) < 4
    Y = squeeze(Xfilt(xInd,yInd,zInd));
else
    q.t = 1:dims(4);
    tInd = q.t(find(ax.stime >= tR(1)):find(ax.stime >= tR(2)));
    Y = squeeze(Xfilt(xInd,yInd,zInd,tInd));
end

if length(size(Y)) > 2
    errordlg('Too many dimensions; check ranges')
    return
end
% if handles.med_box.Value == 1
%     Y = medfilt2(Y,[3 3]);
% end

if handles.hotcold.Value == 1
    h = hotcoldDB;
elseif handles.graybox.Value == 1
    h = 'gray';
else
    h = 'hot';
end

if handles.use_ext_fig.Value == 0
    axes(handles.axes1)
    if handles.plotbox1.Value == 1
        imagesc(ax.stime(tInd),ax.x(xInd),(Y));
        colormap(h)
        if ~isempty(aeR)
            caxis(aeR)
        end
        handles.axes1.XLabel.String = 'Time (ms)';
        handles.axes1.YLabel.String = 'Lateral (mm)';
    end
    if handles.plotbox1.Value == 2
        imagesc(ax.x(xInd),ax.y(yInd),(Y'),'ButtonDownFcn',{@Plot4OnClickXY,handles})
        colormap(h)
        if ~isempty(aeR)
            caxis(aeR)
        end
        handles.axes1.XLabel.String = 'Lateral (mm)';
        handles.axes1.YLabel.String = 'Elevational (mm)';
    end
    if handles.plotbox1.Value == 3
        imagesc(ax.x(xInd),ax.depth(zInd),(Y'),'ButtonDownFcn',{@Plot4OnClickXZ,handles})
        colormap(h)
        if ~isempty(aeR)
            caxis(aeR)
        end
        handles.axes1.XLabel.String = 'Lateral (mm)';
        handles.axes1.YLabel.String = 'Depth (mm)';
    end
    if handles.plotbox1.Value == 4
        imagesc(ax.stime(tInd),ax.y(yInd),(Y))
        colormap(h)
        if ~isempty(aeR)
            caxis(aeR)
        end
        handles.axes1.XLabel.String = 'Time (ms)';
        handles.axes1.YLabel.String = 'Elevational (mm)';
    end
    if handles.plotbox1.Value == 5
        imagesc(ax.y(yInd),ax.depth(zInd),(Y'),'ButtonDownFcn',{@Plot4OnClickYZ,handles})
        colormap(h)
        if ~isempty(aeR)
            caxis(aeR)
        end
        handles.axes1.XLabel.String = 'Elevational (mm)';
        handles.axes1.YLabel.String = 'Depth (mm)';
    end
    if handles.plotbox1.Value == 6
        imagesc(ax.stime(tInd),ax.depth(zInd),Y,'ButtonDownFcn',{@Plot4OnClickTZ,handles})
        colormap(h)
        if ~isempty(aeR)
            caxis(aeR)
        end
        handles.axes1.XLabel.String = 'Time (ms)';
        handles.axes1.YLabel.String = 'Depth (mm)';
    end
else
    figure(2)
    imshow(Y')
    colormap(gca,h)
    if ~isempty(aeR)
        caxis(aeR)
    end
end

%assignin('base','ax',ax);
% if handles.save_fig.Value == 1
%     saveas('figure.png',handles.axes1)
% end





% if handles.plotbox1.Value == 1
%     handles.plotbox1.CData = squeeze(XdB(:,
%     handles.xR
%




function fname_Callback(hObject, eventdata, handles)
% hObject    handle to fname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fname as text
%        str2double(get(hObject,'String')) returns contents of fname as a double


% --- Executes during object creation, after setting all properties.
function fname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_data.
function load_data_Callback(hObject, eventdata, handles) %@003
% hObject    handle to load_data (see GCBO)
[f,  p] = uigetfile(fullfile(pwd,'*4d_data.mat'));
cd(p)
fprintf('Loading 4D Dataset...')
load([p f]);
%[~, ax] = make_axes(param,size(Xfilt),[1 2],1);
set(handles.fname,'String',file);
fprintf('Done\n')
assignin('base','Xfilt',Xfilt);
assignin('base','fpath',[path file]);
assignin('base','param',param);
assignin('base','ax',ax);
assignin('base','LF',LF);
if handles.onemhz.Value == 0
    assignin('base','PEparam',PE);
else
    if exist('PE','var')
        errordlg('Uncheck the 1MHz box')
    end
end
set(handles.active_ae,'String',num2str(size(Xfilt)));
set(handles.active_xfilt,'String','AE');
set(handles.LF_chan,'String',num2str(size(LF,2)));
set(handles.tms,'String',num2str([ax.stime(1) ax.stime(end)]));
set(handles.tsamp,'String',num2str([1 length(ax.stime)]));
set(handles.xmm,'String',num2str([ax.x(1) ax.x(end)]));
set(handles.xsamp,'String',num2str([1 length(ax.x)]));
set(handles.ymm,'String',num2str([ax.y(1) ax.y(end)]));
set(handles.ysamp,'String',num2str([1 length(ax.y)]));
set(handles.zmm,'String',num2str([round(ax.depth(1)) round(ax.depth(end))]));
set(handles.zsamp,'String',num2str([1 length(ax.depth)]));

if handles.reset_axes.Value == 1
    set(handles.xR,'String', num2str([ax.x(1) ax.x(end)]));
    set(handles.yR,'String', num2str([ax.y(1) ax.y(end)]));
    set(handles.zR,'String', num2str([ax.depth(1) floor(ax.depth(end))]));
    set(handles.tR,'String', num2str([ax.stime(1) ax.stime(end)]));
    
end



% --- Executes during object creation, after setting all properties.
function load_data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in create_4d.
function create_4d_Callback(hObject, eventdata, handles) % @004
% hObject    handle to create_4d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.newaescan.Value
    [file, path] = uigetfile(fullfile(pwd,'*_info.mat'));
    param = load([path file]);
    param = param.bScanParm;
    [~,LF] = Read_Data([path file],1);
else
    [file, path] = uigetfile(fullfile(pwd,'*_info.dat')); %Gets file location
    param = read_ucsdi_info([path file]); %Gets scan parameters
    [~,~,LF] = read_ucsdi_data([path file],1); %Gets input current waveform
end
cd(path);
PE = [];
t_delay = 5.2;
if handles.trans.Value ~= 4
    path2 = [path(1:end-8) 'PEData\'];
    file2 = uigetfile(fullfile(path2,'*PEParm.mat')); %gets US pulse waveform
    
    PE = open([path2 file2]);
    US = PE.TW.Wvfm1Wy;
    t_delay = length(US)*(1/PE.bScanParm.vsx_fs); %adjusts depth based on length of pulse
    if handles.trans.Value == 3
        t_delay = t_delay + 12;
    end
end
ax.HFfreq = linspace(0,param.daq.HFdaq.fs_MHz,param.daq.HFdaq.pts); %Creates fast frequency axis
ax.LFfreq = linspace(0,param.daq.HFdaq.pulseRepRate_Hz,param.daq.HFdaq.NoBurstTriggers); %creates slow frequency axis


a_full = str2num(handles.hfchans.String);
hf_num = length(a_full);

%**************************************************************************
%Builds 4D Matrix***************~~~~~~~~~~~~~~~~**************

if isempty(handles.hfchans.String)
    a = 2;
    hf_num = 1;
end
for p = 1:hf_num
    if ~isempty(handles.hfchans.String)
        a = a_full(p);
    end
    
    %Gets Raw Data into cell array
    [~, HF1] = full_signal([path file],param,a,handles.onemhz.Value,handles.newaescan.Value); %Gets the raw data
    
    
    %Filters in Fast Time
    if handles.ft_on.Value == 1
        if handles.match_box.Value == 1
            X = w_ae_filt2(param,HF1,PE,1,handles); %Filters in fast time; 0 is match, 1 uses cutoffs
        end
        if handles.match_box.Value == 0
            X = w_ae_filt2(param,HF1,1,0,handles,[str2double(handles.fast_cut1.String) str2double(handles.fast_cut2.String)]); %Filters in fast time; 0 is match, 1 uses cutoffs
        end
    end
    
    LF = padarray(LF,length(LF));
      %Filters in Slow Time
      if handles.st_on.Value == 1
          if ~isempty(handles.slow_cut2.String) && ~isempty(handles.slow_cut1.String) && handles.slow_box.Value == 0
              [X, LF] = w_slow_filt2(param,X,LF,handles.slow_box.Value,[str2double(handles.slow_cut1.String) str2double(handles.slow_cut2.String)]); %Filters in slow time 0 is match, 1 uses cutoffs
          elseif handles.slow_box.Value == 1
              [X, LF] = w_slow_filt2(param, X,LF,handles.slow_box.Value);
          end
      end
      if handles.st_on.Value == 0 & handles.ft_on.Value == 0
          X = HF1;
      end
      
  
    %%%%%%%%%%%%%%%%%%%% Gets new axes for z and t %%%%%%%%%%%%%%%%%%
    % b = waitbar(0,'Matrix Conversion');
    %    waitbar(0,b,'Enveloping and converting to 4D matrix')
    HF = zeros(size(X,1),size(X,2),size(X{1},1),size(X{1},2));
    s = size(X);
    for i = 1:s(1)
        for j = 1:s(2)
            %HF(i,j,:,:) = envelope(real(X{i,j})); %Converts cell array to double
            HF(i,j,:,:) = X{i,j}; %Converts cell array to double
        end
        % waitbar(i/param.velmex.XNStep,b,'Converting to 4D matrix');
        multiWaitbar('Converting to 4D Matrix', i/s(1));
    end
    
    if handles.onemhz.Value == 1
        if param.velmex.XNStep ~= 1 && param.velmex.YNStep ~= 1
            if param.velmex.SlowAxis == 'X'
                for i = 1:size(HF,1)
                    if mod(i,2) == 0
                        HF(i,:,:,:) = fliplr(HF(i,:,:,:));
                    end
                end
            else
                for i = 1:size(HF,2)
                    if mod(i,2) == 0
                        HF(:,i,:,:) = fliplr(HF(:,i,:,:));
                    end
                end
            end
        end
        %     else
        %         if param.velmex.XNStep ~= 1 && param.velmex.YNStep ~= 1
        %              for i = 1:size(HF,2)
        %                     if mod(i,2) == 0
        %                         HF(:,i,:,:) = fliplr(HF(:,i,:,:));
        %                     end
        %              end
        %         end
    end
    
    
    %%%%%%%%
    if isempty(num2str(handles.depR.String))
        qq = [10 round(1.48*param.daq.HFdaq.pts/param.daq.HFdaq.fs_MHz-10)];
    else
        qq = str2num(handles.depR.String);
    end
    
    dims = size(HF);
    [~, ax] = make_axes(param,dims,t_delay,qq, 12.3); %selects range for dB calculation exlcuding the 10mm around each border
    % XdB = real(20*log10(real(HF)./max(max(max(max(real(HF(:,:,M.xT,:))))))));
    % Xfilt = filts2D(XdB,[0 1 12],[0 2 2]);
    
    %%%%% Using depth calculation from SEP
    if handles.onemhz.Value == 0
        lensDelay = PE.Trans.lensCorrection.*PE.PData.Lambda/1.49; %modified to 1.49
        RefPulseOneWay = PE.TW.Wvfm1Wy;
        RefPulse       = resample(RefPulseOneWay,PE.bScanParm.HFSamplingRate,PE.bScanParm.vsx_fs);
        [~,b]          = max(abs(hilbert(RefPulse-mean(RefPulse))));
        RefPulseT      = (0:size(RefPulse,1)-1)/PE.bScanParm.HFSamplingRate;
        peakDelay       = RefPulseT(b);
        z_delay = (lensDelay + peakDelay).*1.485;
        ax.depth = ax.depth - z_delay;
    end
    %%%%%%%%%%%%%
    
    if handles.tc.Value
        tc_params = evalin('base','tc_params');
        dep_sub = tc_params.thickness/1.48 - tc_params.thickness/tc_params.speed;
        ax.depth = ax.depth-dep_sub;
    end
    %     if ~isempty('PE')
    %         for i = 1:length(PE.TXArray)
    %             fd(i) = max(PE.TXArray(i).Delay);
    %         end
    %         form_delay = max(fd);
    %         ax.depth = ax.depth - form_delay;
    %     end
    
    %%%%% Cyl to Cart %%% This is very rough draft
    ptheta = atan((ax.x-mean(ax.x))./PE.TX.focus*PE.PData.Lambda);
    px = ax.x.*cos(ptheta);
    pz = ax.depth.*cos(ptheta)';
    for i = 1:length(ax.x)
        pXind(i) = find(ax.x >= px(i),1);
        XF2(i,:,:,:) = circshift(HF(i,:,:,:),abs(round(length(ax.x)/2)-i)*-1,3);
        for j = 1:length(ax.depth)
            pZind(i,j) = find(ax.depth >= pz(i,j),1);
        end
    end
    
    if length(size(XF2)) < 4
        HF = permute(XF2,[1 4 2 3]);
    else
        HF = XF2;
    end
    % for m = 1:size(HF,4)
    %     Pfin = zeros(length(pXind),length(pZind));
    %     Pnum = Pfin;
    %    % for i = 1:length(pXind)
    %    for i = 1:length(ax.x)
    %         for j = 1:length(pZind)
    %             Pfin(i,pZind(i,j)) = Pfin(i,pZind(i,j)) + squeeze(HF(i,1,j,m));
    %             Pnum(i,pZind(i,j)) = Pnum(i,pZind(i,j)) + 1;
    %         end
    %     end
    %     Ptot = Pfin./Pnum;
    %     Ptot(isnan(Ptot)) = 0;
    %     XF(:,:,m) = Ptot;
    %
    % end
    
    %%%%%
    Xfilt = HF;
    multiWaitbar('CLOSEALL');
    % delete(b)
    
    if handles.save_4d.Value == 1
        multiWaitbar('Saving','busy');
        f = file(1:end-4);
        if ~isempty(handles.hfchans.String)
            hchan = num2str(a);
            f2 = [f '_chan_' hchan '_4d_data.mat'];
        else
            f2 = [f '_4d_data.mat'];
        end
        fprintf('Saving 4D file...')
        if handles.large_box.Value == 1
            clearvars -except Xfilt file path param ax LF PE f2
            eval([ 'save ' f2 ' -v7.3']);
        else
            %  clearvars -except Xfilt file path param ax LF PE f2
            save(f2,'Xfilt','path','file','param','ax','LF','PE','f2');
        end
        clearvars -except a_full hf_num file path file2 path2 PE US LF param handles
        fprintf('Done\n')
        multiWaitbar('Saving','Close');
    end
end




% --- Executes on button press in bad.
function bad_Callback(hObject, eventdata, handles)
% hObject    handle to bad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bad


% --- Executes on button press in save_4d.
function save_4d_Callback(hObject, eventdata, handles)
% hObject    handle to save_4d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_4d


% --- Executes on button press in reset_axes.
function reset_axes_Callback(hObject, eventdata, handles)
% hObject    handle to reset_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of reset_axes


% --- Executes on button press in movie_button.
function movie_button_Callback(hObject, eventdata, handles) %@005
% hObject    handle to movie_button (see GCBO)
if handles.use_chop.Value == 0
    Xfilt = evalin('base','Xfilt');
    param = evalin('base','param');
    ax = evalin('base','ax');
else
    Xfilt = evalin('base','X_c');
    param = evalin('base','param');
    ax = evalin('base','ax_c');
end

Xfilt = real(Xfilt);
if handles.showlf.Value == 1
    LF = evalin('base','LF');
    LF = LF(:,str2double(handles.LF_chan.String));
    lf_ax = linspace(0,30,length(LF));
end

xR = str2num(handles.xR.String);
if length(xR) == 1
    xR = [xR xR];
end
yR = str2num(handles.yR.String);
if length(yR) == 1
    yR = [yR yR];
end
zR = str2num(handles.zR.String);
if length(zR) == 1
    zR = [zR zR];
end
tR = str2num(handles.tR.String);
if length(tR) == 1
    tR = [tR tR];
end
aeR = str2num(handles.aeR.String);
dims = size(Xfilt);
%[~,ax] = make_axes(param,dims,[1 2],1);
q.x = 1:dims(1);
q.y = 1:dims(2);
q.z = 1:dims(3);
if length(dims) == 4
    q.t = 1:dims(4);
else
    q.t = 1;
end
xInd = q.x(find(ax.x >= xR(1)):find(ax.x >= xR(2)));
yInd = q.y(find(ax.y >= yR(1)):find(ax.y >= yR(2)));
zInd = q.z(find(ax.depth >= zR(1)):find(ax.depth >= zR(2)));
tInd = q.t(find(ax.stime >= tR(1)):find(ax.stime >= tR(2)));

xP = str2double(handles.xP.String);
yP = str2double(handles.yP.String);
zP = str2double(handles.zP.String);
tP = str2double(handles.tP.String);
xpoint = find(ax.x >= xP,1);
ypoint = find(ax.y >= yP,1);
zpoint = find(ax.depth >= zP,1);
tpoint = find(ax.stime >= tP,1);



if handles.showlf.Value == 1
    lfInd = find(lf_ax >= tR(1)):find(lf_ax >= tR(2));
    lfdif = length(lfInd)/length(tInd);
end


if handles.use_ext_fig.Value == 0
    axes(handles.axes2)
else
    figure(1);
end
if handles.hotcold.Value == 1
    h = hotcoldDB;
elseif handles.graybox.Value == 1
    h = 'gray';
else
    h = 'hot';
end
if handles.plotbox2.Value == 1 && handles.all_movie.Value == 0
    if handles.save_fig.Value == 0
        for k = tInd %Mod loop
            %             if handles.med_box.Value == 1
            %                 J = medfilt2((squeeze(Xfilt(xInd,ypoint,zInd,k)))',[5 5]);
            %             else
            J = squeeze(Xfilt(xInd,ypoint,zInd,k))';
            %             end
            
            if handles.use_ext_fig.Value == 0
                imagesc(ax.x(xInd),ax.depth(zInd),J) % mod plots
                colormap(h)
                if handles.showlf.Value == 1
                    
                    plot(handles.axes4,lf_ax(lfInd),LF(lfInd),'k')
                    hold(handles.axes4,'on')
                    plot(handles.axes4,lf_ax(lfInd(round(k*lfdif))),LF(lfInd(round(k*lfdif))),'ro','MarkerFaceColor','r')
                    hold(handles.axes4,'off')
                    
                end
                
            else
                imshow(J)
                colormap(h)
            end
            if ~isempty(aeR)
                caxis(aeR)
            end
            title(['t = ' num2str(ax.stime(k))]);
            % text(15,15,['t = ' num2str(ax.stime(k))]);
            handles.axes2.XLabel.String = 'Lateral (mm)'; %Mod Axes
            handles.axes2.YLabel.String = 'Depth (mm)';
            drawnow
        end
    elseif handles.save_fig.Value == 1
        vidwrite(param,ax,Xfilt,handles)
    end
    
elseif handles.plotbox2.Value == 2 && handles.all_movie.Value == 0
    if handles.save_fig.Value == 0
        for k = tInd %Mod loop
            %             if handles.med_box.Value == 1
            %                 J = medfilt2((squeeze(Xfilt(xpoint,yInd,zInd,k)))',[5 5]);
            %             else
            J = squeeze(Xfilt(xpoint,yInd,zInd,k))';
            %             end
            
            if handles.use_ext_fig.Value == 0
                imagesc(ax.y(yInd),ax.depth(zInd),J) % mod plots
                %   h = hotcoldDB;
                colormap(h)
                if handles.showlf.Value == 1
                    
                    plot(handles.axes4,lf_ax(lfInd),LF(lfInd),'k')
                    hold(handles.axes4,'on')
                    plot(handles.axes4,lf_ax(lfInd(round(k*lfdif))),LF(lfInd(round(k*lfdif))),'ro','MarkerFaceColor','r')
                    hold(handles.axes4,'off')
                    
                end
                
            else
                imshow(J)
                colormap(h)
            end
            if ~isempty(aeR)
                caxis(aeR)
            end
            title(['t = ' num2str(ax.stime(k))]);
            % text(15,15,['t = ' num2str(ax.stime(k))]);
            handles.axes2.XLabel.String = 'Elevational (mm)'; %Mod Axes
            handles.axes2.YLabel.String = 'Depth (mm)';
            drawnow
        end
        
        
    elseif handles.save_fig.Value == 1
        vidwrite(param,ax,Xfilt,handles)
    end
    
elseif handles.plotbox2.Value == 4 && handles.all_movie.Value == 0
    if handles.save_fig.Value == 0
        for k = tInd %Mod loop
            %             if handles.med_box.Value == 1
            %                 J = medfilt2((squeeze(Xfilt(xInd,yInd,zpoint,k)))',[5 5]);
            %             else
            J = squeeze(Xfilt(xInd,yInd,zpoint,k))';
            %             end
            
            if handles.use_ext_fig.Value == 0
                imagesc(ax.x(xInd),ax.y(yInd),J) % mod plots
                %   h = hotcoldDB;
                colormap(h)
                if handles.showlf.Value == 1
                    
                    plot(handles.axes4,lf_ax(lfInd),LF(lfInd),'k')
                    hold(handles.axes4,'on')
                    plot(handles.axes4,lf_ax(lfInd(round(k*lfdif))),LF(lfInd(round(k*lfdif))),'ro','MarkerFaceColor','r')
                    hold(handles.axes4,'off')
                    
                end
                
            else
                imshow(J)
                colormap(h)
            end
            if ~isempty(aeR)
                caxis(aeR)
            end
            title(['t = ' num2str(ax.stime(k))]);
            % text(15,15,['t = ' num2str(ax.stime(k))]);
            handles.axes2.XLabel.String = 'Lateral (mm)'; %Mod Axes
            handles.axes2.YLabel.String = 'Elevational (mm)';
            drawnow
        end
        
        
    elseif handles.save_fig.Value == 1
        vidwrite(param,ax,Xfilt,handles)
    end
    
    
elseif handles.plotbox2.Value == 3 && handles.all_movie.Value == 0
    if handles.save_fig.Value == 0
        for k = zInd %Mod loop
            %             if handles.med_box.Value == 1
            %                 J = medfilt2((squeeze(Xfilt(xInd,yInd,k,tpoint)))',[5 5]);
            %             else
            J = squeeze(Xfilt(xInd,yInd,k,tpoint))';
            %             end
            
            if handles.use_ext_fig.Value == 0
                imagesc(ax.x(xInd),ax.y(yInd),J) % mod plots
                %     h = hotcoldDB;
                colormap(h)
            else
                imshow(J)
                colormap(h)
            end
            if ~isempty(aeR)
                caxis(aeR)
            end
            title(['z = ' num2str(ax.depth(k))]);
            % text(15,15,['t = ' num2str(ax.stime(k))]);
            handles.axes2.XLabel.String = 'Lateral (mm)'; %Mod Axes
            handles.axes2.YLabel.String = 'Elevational (mm)';
            drawnow
        end
        
        
    elseif handles.save_fig.Value == 1
        vidwrite(param,ax,Xfilt,handles)
    end
end
if handles.all_movie.Value == 1
    
    for k = tInd
        J1 = squeeze(Xfilt(xInd,ypoint,zInd,k))';
        J2 = squeeze(Xfilt(xpoint,yInd,zInd,k))';
        J3 = squeeze(Xfilt(xInd,yInd,zpoint,k))';
        
        axes(handles.axes1)
        imagesc(ax.x(xInd),ax.depth(zInd),J1) % mod plots
        colormap(h)
        if ~isempty(aeR)
            caxis(aeR)
        end
        
        axes(handles.axes3)
        imagesc(ax.y(yInd),ax.depth(zInd),J2)
        colormap(h)
        if ~isempty(aeR)
            caxis(aeR)
        end
        
        axes(handles.axes2)
        imagesc(ax.x(xInd),ax.y(yInd),J3)
        colormap(h)
        if ~isempty(aeR)
            caxis(aeR)
        end
        drawnow
        %                     handles.axes1.XLabel.String = 'Lateral (mm)'; %Mod Axes
        %             handles.axes1.YLabel.String = 'Depth (mm)';
        %                         handles.axes3.XLabel.String = 'Elevational (mm)'; %Mod Axes
        %             handles.axes3.YLabel.String = 'Depth (mm)';
        %                         handles.axes2.XLabel.String = 'Lateral (mm)'; %Mod Axes
        %             handles.axes2.YLabel.String = 'Elevational (mm)';
        
        
        if handles.showlf.Value == 1
            
            plot(handles.axes4,lf_ax(lfInd),LF(lfInd),'k')
            hold(handles.axes4,'on')
            %  plot(handles.axes4,lf_ax(lfInd(round(k*lfdif))),LF(lfInd(round(k*lfdif))),'ro','MarkerFaceColor','r')
            plot(handles.axes4,lf_ax(round(k*lfdif)),LF(round(k*lfdif)),'ro','MarkerFaceColor','r')
            hold(handles.axes4,'off')
            
        end
    end
end



% --- Executes on button press in save_fig.
function save_fig_Callback(hObject, eventdata, handles)
% hObject    handle to save_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_fig



function aeR_Callback(hObject, eventdata, handles)
% hObject    handle to aeR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of aeR as text
%        str2double(get(hObject,'String')) returns contents of aeR as a double


% --- Executes during object creation, after setting all properties.
function aeR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to aeR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in mean_box.
function mean_box_Callback(hObject, eventdata, handles)
% hObject    handle to mean_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of mean_box



function mean_x_Callback(hObject, eventdata, handles)
% hObject    handle to mean_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mean_x as text
%        str2double(get(hObject,'String')) returns contents of mean_x as a double


% --- Executes during object creation, after setting all properties.
function mean_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mean_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mean_y_Callback(hObject, eventdata, handles)
% hObject    handle to mean_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mean_y as text
%        str2double(get(hObject,'String')) returns contents of mean_y as a double


% --- Executes during object creation, after setting all properties.
function mean_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mean_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function mean_z_Callback(hObject, eventdata, handles)
% hObject    handle to mean_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mean_z as text
%        str2double(get(hObject,'String')) returns contents of mean_z as a double


% --- Executes during object creation, after setting all properties.
function mean_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mean_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function int_y_Callback(hObject, eventdata, handles)
% hObject    handle to int_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of int_y as text
%        str2double(get(hObject,'String')) returns contents of int_y as a double


% --- Executes during object creation, after setting all properties.
function int_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to int_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function int_z_Callback(hObject, eventdata, handles)
% hObject    handle to int_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of int_z as text
%        str2double(get(hObject,'String')) returns contents of int_z as a double


% --- Executes during object creation, after setting all properties.
function int_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to int_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in int_box.
function int_box_Callback(hObject, eventdata, handles)
% hObject    handle to int_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of int_box



function int_x_Callback(hObject, eventdata, handles)
% hObject    handle to int_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of int_x as text
%        str2double(get(hObject,'String')) returns contents of int_x as a double


% --- Executes during object creation, after setting all properties.
function int_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to int_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Enhance_Sig.
function Enhance_Sig_Callback(hObject, eventdata, handles) % @006
% hObject    handle to Enhance_Sig (see GCBO)
if handles.use_chop.Value == 0
    param = evalin('base','param');
    Xfilt = evalin('base','Xfilt');
    m = [handles.mean_box.Value str2double(handles.mean_x.String) str2double(handles.mean_y.String) str2double(handles.mean_z.String)];
    n = [handles.int_box.Value str2double(handles.int_x.String) str2double(handles.int_y.String) str2double(handles.int_z.String) get(handles.squarify_box,'Value')];
    o = [handles.med_box.Value str2double(handles.med_x.String) str2double(handles.med_y.String) str2double(handles.med_z.String)];
    Xfilt = filts3D(Xfilt,m,n,o,param);
    [~,ax] = make_axes(param,size(Xfilt));
    assignin('base','ax',ax);
    assignin('base','Xfilt',Xfilt)
else
    param = evalin('base','param');
    Xfilt = evalin('base','X_c');
    ax = evalin('base','ax_c');
    m = [handles.mean_box.Value str2double(handles.mean_x.String) str2double(handles.mean_y.String) str2double(handles.mean_z.String)];
    n = [handles.int_box.Value str2double(handles.int_x.String) str2double(handles.int_y.String) str2double(handles.int_z.String) get(handles.squarify_box,'Value')];
    o = [handles.med_box.Value str2double(handles.med_x.String) str2double(handles.med_y.String) str2double(handles.med_z.String)];
    Xfilt = filts3D(Xfilt,m,n,o,param);
    
    dims = size(Xfilt);
    xR = [ax.x(1) ax.x(end)];
    yR = [ax.y(1) ax.y(end)];
    zR = [ax.depth(1) ax.depth(end)];
    tR = [ax.stime(1) ax.stime(end)];
    ax.depth = linspace(zR(1),zR(2),dims(3));
    if length(dims) > 3
        ax.stime = linspace(tR(1),tR(2),dims(4));
    else
        ax.stime = 1;
    end
    if mean(abs(xR)) > 1
        ax.x = linspace(xR(1),xR(2),dims(1));
    else
        ax.x = 1;
    end
    if mean(abs(yR)) > 1
        ax.y = linspace(yR(1),yR(2),dims(2));
    else
        ax.y = 1;
    end
    
    assignin('base','ax_c',ax);
    assignin('base','X_c',Xfilt)
end



% --- Executes on button press in chop.
function chop_Callback(hObject, eventdata, handles) %@007
% hObject    handle to chop (see GCBO)
Xfilt = evalin('base','Xfilt');
param = evalin('base','param');
ax = evalin('base','ax');
clear X_c
xR = str2num(handles.xR.String);
if length(xR) == 1
    xR(2) = xR(1);
end
yR = str2num(handles.yR.String);
zR = str2num(handles.zR.String);
tR = str2num(handles.tR.String);
if length(yR) == 1
    yR(2) = yR(1);
end
if length(tR) == 1
    tR(2) = tR(1);
end
if length(zR) == 1
    zR(2) = zR(1);
end
aeR = str2num(handles.aeR.String);
dims = size(Xfilt);
%[~,ax] = make_axes(param,dims,[1 2],1);
q.x = 1:dims(1);
q.y = 1:dims(2);
q.z = 1:dims(3);
if length(size(Xfilt)) > 3
    q.t = 1:dims(4);
    tInd = q.t(find(ax.stime >= tR(1)):find(ax.stime >= tR(2)));
else
    tInd = 1;
end
xInd = q.x(find(ax.x >= xR(1)):find(ax.x >= xR(2)));
yInd = q.y(find(ax.y >= yR(1)):find(ax.y >= yR(2)));
zInd = q.z(find(ax.depth >= zR(1),1):find(ax.depth >= zR(2),1));

X = Xfilt(xInd,yInd,zInd,tInd);
if length(tInd) == 1
    X = permute(X,[1 2 3 4]);
end
%[~,ax] = make_axes(param,size(X));
dims = size(X);
if length(dims) < 4
    if length(tInd) == 1 && length(zInd) == 1
        dims = [dims(1) dims(2) 1 1];
    else
        dims = [dims(1) dims(2) dims(3) 1];
    end
end


ax.depth = linspace(zR(1),zR(2),dims(3));
ax.stime = linspace(tR(1),tR(2),dims(4));
if mean(abs(xR)) > 1
    ax.x = linspace(xR(1),xR(2),dims(1));
else
    ax.x = 1;
end
if mean(abs(yR)) > 1
    ax.y = linspace(yR(1),yR(2),dims(2));
else
    ax.y = 1;
end
if handles.active_xfilt.String == 'PE'
    set(handles.active_pe,'String',num2str(size(X)))
else
    set(handles.active_ae,'String',num2str(size(X)))
end
assignin('base','ax_c',ax);
assignin('base','X_c',X);


% --- Executes on button press in max_box.
function max_box_Callback(hObject, eventdata, handles) %@008
% hObject    handle to max_box (see GCBO)
if get(hObject,'Value') == 1
    if handles.use_chop.Value == 0
        ax = evalin('base','ax');
        x1 = num2str(round(ax.x(1),1)); x2 = num2str(round(ax.x(end),1));
        z2 = floor(ax.depth(end));
        set(handles.xR,'String', [x1 ' ' x2]);
        set(handles.yR,'String', num2str([ax.y(1) ax.y(end)]));
        set(handles.zR,'String', num2str([0 z2]));
        set(handles.tR,'String', num2str([ax.stime(1) ax.stime(end)]));
    else
        ax_c = evalin('base','ax_c');
        set(handles.xR,'String', num2str([ax_c.x(1) ax_c.x(end)]));
        set(handles.yR,'String', num2str([ax_c.y(1) ax_c.y(end)]));
        set(handles.zR,'String', num2str([ax_c.depth(1) floor(ax_c.depth(end))]));
        set(handles.tR,'String', num2str([ax_c.stime(1) ax_c.stime(end)]));
    end
end

% Hint: get(hObject,'Value') returns toggle state of max_box


% --- Executes on button press in use_chop.
function use_chop_Callback(hObject, eventdata, handles)
% hObject    handle to use_chop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_chop


% --- Executes on button press in match_box.
function match_box_Callback(hObject, eventdata, handles)
% hObject    handle to match_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of match_box



function fast_cut1_Callback(hObject, eventdata, handles)
% hObject    handle to fast_cut1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fast_cut1 as text
%        str2double(get(hObject,'String')) returns contents of fast_cut1 as a double


% --- Executes during object creation, after setting all properties.
function fast_cut1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fast_cut1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fast_cut2_Callback(hObject, eventdata, handles)
% hObject    handle to fast_cut2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fast_cut2 as text
%        str2double(get(hObject,'String')) returns contents of fast_cut2 as a double


% --- Executes during object creation, after setting all properties.
function fast_cut2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fast_cut2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in slices.
function slices_Callback(hObject, eventdata, handles)
% hObject    handle to slices (see GCBO)
if handles.use_chop.Value == 0
    X = evalin('base','Xfilt');
    ax = evalin('base','ax');
else
    X = evalin('base','X_c');
    ax = evalin('base','ax_c');
end
tR = str2num(handles.tR.String);
if length(tR) > 1;
    tR = tR(1);
end
aeR = str2num(handles.aeR.String);
% zR = str2num(handles.zR.String);
% xR = str2num(handles.xR.String);
% yR = str2num(handles.yR.String);
% dims = size(X);
% %[~,ax] = make_axes(param,dims);
%
% ax.depth = linspace(zR(1),zR(2),dims(3));
% ax.stime = linspace(tR(1),tR(2),dims(4));
% if mean(abs(xR)) > 1
%     ax.x = linspace(xR(1),xR(2),dims(1));
% else
%     ax.x = 1;
% end
% if mean(abs(yR)) > 1
%     ax.y = linspace(yR(1),yR(2),dims(2));
% else
%     ax.y = 1;
% end
%
% q.x = 1:dims(1);
% q.y = 1:dims(2);
% q.z = 1:dims(3);
% q.t = 1:dims(4);
t = find(ax.stime >= tR,1);
% xInd = q.x(find(ax.x >= xR(1)):find(ax.x >= xR(2)));
% yInd = q.y(find(ax.y >= yR(1)):find(ax.y >= yR(2)));
% zInd = q.z(find(ax.depth >= zR(1)):find(ax.depth >= zR(2)));
% tInd = q.t(find(ax.stime >= tR(1)):find(ax.stime >= tR(2)));

X = squeeze(X(:,:,:,t));
[x,y,z] = meshgrid(linspace(ax.y(1),ax.y(end),length(ax.y)),linspace(ax.x(1),ax.x(end),length(ax.x)),...
    linspace(ax.depth(1),ax.depth(end),length(ax.depth)));
yslice = ax.y(end);
%zslice = [ax.depth(1),median(ax.depth),ax.depth(end)];
zslice = ax.depth(1);
axes(handles.axes2)
for i = ax.x
    slice(x,y,z,X,i,yslice,zslice);
    if ~isempty(aeR)
        caxis(aeR)
    end
    drawnow
end



function slow_cut1_Callback(hObject, eventdata, handles)
% hObject    handle to slow_cut1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slow_cut1 as text
%        str2double(get(hObject,'String')) returns contents of slow_cut1 as a double


% --- Executes during object creation, after setting all properties.
function slow_cut1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slow_cut1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function slow_cut2_Callback(hObject, eventdata, handles)
% hObject    handle to slow_cut2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of slow_cut2 as text
%        str2double(get(hObject,'String')) returns contents of slow_cut2 as a double


% --- Executes during object creation, after setting all properties.
function slow_cut2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slow_cut2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in env_button.
function env_button_Callback(hObject, eventdata, handles) %@009
% hObject    handle to env_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.use_chop.Value == 0
    HF = evalin('base','Xfilt');
else
    HF = evalin('base','X_c');
end
dims = (size(HF));
%%%TESTING THIS OUT
b = waitbar(0);
for i = 1:dims(1)
    for j = 1:dims(2)
        % for k = 1:dims(4)
        HF(i,j,:,:) = envelope(squeeze(real(HF(i,j,:,:))));
        %end
    end
    waitbar(i/dims(1),b,'2D Enveloping');
end
delete(b)
if handles.use_chop.Value == 0
    assignin('base','Xfilt',HF);
else
    assignin('base','X_c',HF);
end




% --- Executes on button press in slow_box.
function slow_box_Callback(hObject, eventdata, handles)
% hObject    handle to slow_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of slow_box



function depR_Callback(hObject, eventdata, handles)
% hObject    handle to depR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of depR as text
%        str2double(get(hObject,'String')) returns contents of depR as a double


% --- Executes during object creation, after setting all properties.
function depR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to depR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in irad.
function irad_Callback(hObject, eventdata, handles) %@010
param = evalin('base','param');
if handles.use_chop.Value == 1
    HF = evalin('base','X_c');
    ax = evalin('base','ax_c');
else
    HF = evalin('base','Xfilt');
    ax = evalin('base','ax');
end

HF = squeeze(HF); %Converts down to lateral+depth+time
t = str2num(handles.tR.String); %gets range of time points
y = str2num(handles.yR.String);

if length(size(HF)) == 2
    HF = permute(HF,[2 1]);
    R = iradon(HF,str2double(handles.dtheta.String),'None');
    R = permute(R,[2 1]);
    dims = size(R);
    ax.depth = linspace(0,ax.depth(end),dims(2));
    ax.x = linspace(ax.x(1),ax.x(end),dims(1));
    low_x = find(ax.x >= min(ax.x),1);
    high_x = find(ax.x >= max(ax.x),1);
    R = R(low_x:high_x,:,:,:);
    ax.x = linspace(ax.x(low_x),ax.x(high_x),size(R,1));
    R = permute(R,[1 3 2]);
    
elseif length(size(HF)) == 3
    q.t = 1:size(HF,3);
    tInd =  q.t(find(ax.stime >= t(1)):find(ax.stime >= t(2)));
    HF = permute(HF,[2 1 3]);
    b = waitbar(0);
    
    for i = tInd
        R(:,:,i) = iradon(HF(:,:,i),str2double(handles.dtheta.String));
        waitbar(i/tInd(end),b,'Computing inverse radon transform');
    end
    
    R = permute(R,[1 4 2 3]);
    delete(b);
    dims = size(R); %gets new dimensions of reconstructed data
    ax.depth = linspace(0,ax.depth(end),dims(3));
    ax.x = linspace(ax.x(1),ax.x(end),dims(1));
    low_x = find(ax.x >= min(ax.x),1);
    high_x = find(ax.x >= max(ax.x),1);
    R = R(low_x:high_x,:,:,:);
    ax.x = linspace(ax.x(low_x),ax.x(high_x),size(R,1));
    
    
else
    q.t = 1:size(HF,4);
    q.y = 1:size(HF,2);
    tInd =  q.t(find(ax.stime >= t(1)):find(ax.stime >= t(2)));
    yInd = q.y(find(ax.y >= y(1)):find(ax.y >= y(2)));
    HF = permute(HF,[3 1 2 4]);
    b = waitbar(0);
    
    for j = yInd
        for i = tInd
            R(:,:,j,i) = iradon(HF(:,:,j,i),str2double(handles.dtheta.String),'None');
            waitbar((j)/yInd(end),b,'Computing inverse radon transform');
        end
    end
    
    
    R = permute(R,[2 3 1 4]);
    delete(b);
    dims = size(R); %gets no dimensions of reconstructed data
    ax.depth = linspace(0,ax.depth(end),dims(3));
    ax.x = linspace(ax.x(1),ax.x(end),dims(1));
    ax.y = linspace(ax.y(1),ax.y(end),dims(2));
end


if handles.use_chop.Value == 1
    assignin('base','ax_c',ax);
    assignin('base','X_c',R);
else
    assignin('base','ax',ax);
    assignin('base','Xfilt',R);
end







function dtheta_Callback(hObject, eventdata, handles)
% hObject    handle to dtheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dtheta as text
%        str2double(get(hObject,'String')) returns contents of dtheta as a double


% --- Executes during object creation, after setting all properties.
function dtheta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dtheta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dB_Button.
function dB_Button_Callback(hObject, eventdata, handles) %@011
% hObject    handle to dB_Button (see GCBO)
if handles.use_chop.Value == 0
    HF = evalin('base','Xfilt');
    ax = evalin('base','ax');
else
    HF = evalin('base','X_c');
    ax = evalin('base','ax_c');
end
U = length(ax.depth);
if isempty(num2str(handles.depR.String))
    qq = [ax.depth(floor(U./5)) ax.depth(floor(U.*0.9))];
else
    qq = str2num(handles.depR.String);
end
dims = (size(HF));
q.z = 1:size(HF,3);
zInd =  q.z(find(ax.depth >= qq(1)):find(ax.depth >= qq(2)));
HFdB = 20*log10(HF./(max(max(max(max(HF(:,:,zInd,:)))))));
fprintf('Finished converting to dB\n');
if handles.use_chop.Value == 0
    assignin('base','Xfilt',HFdB);
else
    assignin('base','X_c',HFdB);
end


% --- Executes on button press in use_ext_fig.
function use_ext_fig_Callback(hObject, eventdata, handles)
% hObject    handle to use_ext_fig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_ext_fig


% --- Executes on selection change in LF_Box.
function LF_Box_Callback(hObject, eventdata, handles)
% hObject    handle to LF_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns LF_Box contents as cell array
%        contents{get(hObject,'Value')} returns selected item from LF_Box


% --- Executes during object creation, after setting all properties.
function LF_Box_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LF_Box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LF_chan_Callback(hObject, eventdata, handles)
% hObject    handle to LF_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LF_chan as text
%        str2double(get(hObject,'String')) returns contents of LF_chan as a double


% --- Executes during object creation, after setting all properties.
function LF_chan_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LF_chan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in LF_FFT.
function LF_FFT_Callback(hObject, eventdata, handles)
% hObject    handle to LF_FFT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LF_FFT


% --- Executes on button press in LF_butt.
function LF_butt_Callback(hObject, eventdata, handles) %@012
param = evalin('base','param');
LF = evalin('base','LF');
LF = LF(:,str2double(handles.LF_chan.String));
LF = LF/str2double(handles.lfgain.String)*1000;
if handles.LF_FFT.Value == 1
    lf = fft(LF);
    x = linspace(0,param.daq.LFdaq.fs_Hz,length(lf));
    if handles.use_ext_fig.Value == 1
        figure(33)
        plot(x,abs(lf))
    else
        axes(handles.axes3)
        plot(x,abs(lf))
    end
    if ~isempty(handles.xlims3.String)
        xlim(str2num(handles.xlims3.String));
    end
    if ~isempty(handles.ylims3.String)
        ylim(str2num(handles.ylims3.String));
    end
    
end
x = linspace(0,param.daq.HFdaq.duration_ms,length(LF));
if handles.use_ext_fig.Value == 1
    figure(3)
    plot(x,LF)
else
    axes(handles.axes2)
    plot(x,LF)
end
ylabel('mA');
xlabel('ms');
if ~isempty(handles.xlims.String)
    xlim(str2num(handles.xlims.String));
end
if ~isempty(handles.ylims.String)
    ylim(str2num(handles.ylims.String));
end




function xlims_Callback(hObject, eventdata, handles)
% hObject    handle to xlims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xlims as text
%        str2double(get(hObject,'String')) returns contents of xlims as a double


% --- Executes during object creation, after setting all properties.
function xlims_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xlims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function ylims_Callback(hObject, eventdata, handles)
% hObject    handle to ylims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ylims as text
%        str2double(get(hObject,'String')) returns contents of ylims as a double


% --- Executes during object creation, after setting all properties.
function ylims_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ylims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%% Add to here at end %%%%%%%%%%% 05/22/18



% --- Executes on button press in sense_button.
function sense_button_Callback(hObject, eventdata, handles) %@013
% hObject    handle to sense_button (see GCBO)
%datacursormode on
%dcm = datacursormode(gcf);
if handles.keep_box.Value == 0
    [x,y] = ginput(1);
else
    x = str2double(get(handles.tP,'String'));
    y = str2double(get(handles.zP,'String'));
end
    
if handles.use_chop.Value == 0
    Xfilt = evalin('base','Xfilt');
    param = evalin('base','param');
    ax = evalin('base','ax');
else
    Xfilt = evalin('base','X_c');
    param = evalin('base','param');
    ax = evalin('base','ax_c');
end
dims = size(Xfilt);
LF = evalin('base','LF');
lchan = str2double(handles.LF_chan.String);
LF = LF(:,lchan);
LF_axis = linspace(0,param.daq.HFdaq.duration_ms,param.daq.LFdaq.pts);
q.lf = 1:length(LF);
C = handles.axes1.Children.CData; %Gets image data, Y is dim 1, X is dim 2
xax = handles.axes1.Children.XData;
yax = handles.axes1.Children.YData;
%Get whether from from space varying or time varying plot
if handles.plotbox1.Value == 1 || handles.plotbox1.Value == 4 || handles.plotbox1.Value == 6
    yloc = find(yax >= y,1);
    S = C(yloc,:);
    tR = str2num(handles.tR.String);
    if length(tR) <2
        errordlg('Time input must be a range');
        return
    end
    q.t = 1:dims(4);
    tInd = q.t(find(ax.stime >= tR(1),1):find(ax.stime >= tR(2),1));
    lfInd = q.lf(find(LF_axis >= tR(1),1):find(LF_axis >= tR(2),1));
    Lae = LF(lfInd);
    Sae2 = resample(S,length(Lae),length(S));
end

if handles.plotbox1.Value == 2 || handles.plotbox1.Value == 3 || handles.plotbox1.Value == 5
    tR = str2num(handles.tR.String);
    if length(tR) <2
        errordlg('Time input must be a range');
        return
    end
    q.t = 1:dims(4);
    tInd = q.t(find(ax.stime >= tR(1),1):find(ax.stime >= tR(2),1));
    lfInd = q.lf(find(LF_axis >= tR(1),1):find(LF_axis >= tR(2),1));
    Lae = LF(lfInd);
    if handles.plotbox1.Value == 3
        %         yR = str2num(handles.yR.String);
        %         if length(yR) > 1
        %             yR = yR(1);
        %         end
        yP = str2double(handles.yP.String);
        sy = find(ax.y >=yP,1);
        sx = find(ax.x>=x,1);
        sz = find(ax.depth>=y,1);
        Sae = squeeze(Xfilt(sx,sy,sz,tInd));
        Sae2 = resample(Sae,length(Lae),length(Sae));  %This might need to be interp1
        %Sae2 = interp1(linspace(0,dims(4),dims(4)),Sae,linspace(0,dims(4),length(Lae)));
        
    end
    if handles.plotbox1.Value == 2
        %         zR = str2num(handles.zR.String);
        %         if length(zR) > 1
        %             zR = zR(1);
        %         end
        zP = str2double(handles.zP.String);
        sZ = find(ax.y >=zP,1);
        sx = find(ax.x>=x,1);
        sy = find(ax.y>=y,1);
        Sae = squeeze(Xfilt(sx,sy,sZ,tInd));
        Sae2 = resample(Sae,length(Lae),length(S));
        
    end
    if handles.plotbox1.Value == 5
        xR = str2num(handles.xR.String);
        if length(xR) > 1
            xR = xR(1);
        end
        xP = str2double(handles.xP.String);
        sX = find(ax.y >=xP,1);
        sx = find(ax.y>=x,1);
        sz = find(ax.depth>=y,1);
        Sae = squeeze(Xfilt(sX,sx,sz,tInd));
        Sae2 = resample(Sae,length(Lae),length(S));
    end
end

Lnorm = Lae-min(Lae);
Lnorm = Lnorm./max(Lnorm);
Snorm = Sae2 - min(Sae2);
Snorm = Snorm./max(Snorm);
T_axis = linspace(tR(1),tR(2),length(Lae));

Sae2 = Sae2/str2double(handles.hfgain.String)*1000000;
Lae = Lae/str2double(handles.lfgain.String)*1000;

%R = corrcoef(abs(Sae2),abs(Lae));
R = corrcoef(Sae2,Lae);
R = R(2);
if size(Sae2,1) ~= size(Lae,1)
    Sae2 = Sae2';
end
fit  = polyfit(Lae,Sae2,1);
taxis = linspace(min(Lae),max(Lae),length(Lae));
yfit = fit(1)*taxis+fit(2);

if handles.use_ext_fig.Value == 1
    figure(6);
    hold off;
    plot(0)
    % scatter(abs(Lae),abs(Sae2));
    scatter(Lae,Sae2,13,'r','filled')
    hold on
    plot(taxis,yfit,'Color','k','LineWidth',2.5)
    hold off
    ylabel('\muV')
    xlabel('mA')
    figure(66)
    hold on
    plot(T_axis,Lnorm,'k')
    plot(T_axis,Snorm,'r')
    title(['R^2 = ' num2str(R)]);
    hold off
    xlabel('ms')
else
    axes(handles.axes2)
    hold off
    plot(0)
    % scatter(abs(Lae),abs(Sae2));
    scatter(Lae,Sae2,13,'r','filled')
    hold on
    plot(taxis,yfit,'Color','k','LineWidth',2.5)
    hold off
    ylabel('\muV')
    xlabel('mA')
    axes(handles.axes3);
    hold off
    plot(0)
    hold on
    plot(T_axis,Lnorm,'k')
    plot(T_axis,Snorm,'r')
    title(['R^2 = ' num2str(R)]);
    xlabel('ms')
    xlim(tR)
    hold off
end

%m = round(mean(abs(Sae2)./abs(Lae)),2);
m = fit(1);
m2 = (max(abs(Sae2))-min(abs(Sae2)))/(max(abs(Lae))-min(abs(Lae)));
ave = round(mean(abs(Sae2)),2);
dev = round(std((Sae2)),2);
if ~isempty(handles.output5.String)
    pres = str2num(handles.output5.String);
else
    pres = 1;
end
%m3 = m/pres;
m3 = m/pres;

set(handles.param1,'String','slope')
set(handles.param2,'String','mean')
set(handles.param3,'String','std')
set(handles.output1,'String',num2str(m));
set(handles.output2,'String',num2str(ave));
set(handles.output3,'String',num2str(dev));
set(handles.output6,'String',num2str(m2));
set(handles.param4,'String','truslp')
set(handles.output4,'String',num2str(m3));
set(handles.param6,'String','slope2')
set(handles.param5,'String','prsr');


p = 7;



% --- Executes on button press in fwhm_button.
function fwhm_button_Callback(hObject, eventdata, handles) %@014
if handles.keep_box.Value == 0
[x,y] = ginput(1);
else
    if handles.plotbox1.Value == 3
        x = str2double(get(handles.xP,'String'));
        y = str2double(get(handles.zP,'String'));
    elseif handles.plotbox1.Value == 5
         x = str2double(get(handles.yP,'String'));
        y = str2double(get(handles.zP,'String'));
    elseif handles.plotbox1.Value == 2
         x = str2double(get(handles.xP,'String'));
        y = str2double(get(handles.yP,'String'));
    end
end

C = handles.axes1.Children.CData; %Gets image data, Y is dim 1, X is dim 2
xax = handles.axes1.Children.XData;
yax = handles.axes1.Children.YData;

[~,yloc] = find(yax>y,1);
S = C(yloc,:);


if handles.use_ext_fig.Value == 1
    figure(5);
    if handles.hold_box.Value == 0
        hold off
        plot(0)
    end
else
    axes(handles.axes2);
    if handles.hold_box.Value == 0
        hold off
        plot(0)
    end
end

if max(S) <= 0
    ydb(1:length(S)) = max(S)-6;
    ylabel('dB')
else
    S = S/str2double(handles.hfgain.String)*1000000;
    ydb(1:length(S)) = max(S)/2;
    ylabel('\muV');
end

hold on
plot(xax,S,'k')
plot(xax,ydb,'r--')
xlabel('mm');

P = round(max(S),2);
f1 = find(S>=ydb(1),1);
f2 = find(S(f1:end) <= ydb(1),1) +f1;
cut1 = round(xax(f1),3);
cut2 = round(xax(f2),3);
% S2 = fliplr(S);
% cutt = (find(S2>=ydb,1)); %Finish this later
% cut2 = round(xax(length(xax)-cutt));


set(handles.param1,'String','peak X')
set(handles.output1,'String',num2str(P));
set(handles.param2,'String','minX')
set(handles.param3,'String','maxX')
set(handles.param7,'String','FWHM X')
set(handles.output2,'String',num2str(cut1));
set(handles.output3,'String',num2str(cut2));
set(handles.output7,'String',num2str(cut2-cut1));

if ~isempty(handles.xlims.String)
    xlim(str2num(handles.xlims.String));
end
if ~isempty(handles.ylims.String)
    ylim(str2num(handles.ylims.String));
end

if handles.use_ext_fig.Value == 1
    figure(55);
    if handles.hold_box.Value == 0
        hold off
        plot(0)
    end
else
    axes(handles.axes3);
    if handles.hold_box.Value == 0
        hold off
        plot(0)
    end
end

hold on
clear S
clear ydb
xloc = find(xax>x,1);
S = C(:,xloc);
if max(S) <= 0
    ydb(1:length(S)) = max(S)-6;
    ylabel('dB')
else
    S = S/str2double(handles.hfgain.String)*1000000;
    ydb(1:length(S)) = max(S)/2;
    ylabel('\muV')
end

plot(yax,S,'k')
plot(yax,ydb,'r--')
xlabel('mm')

P = round(max(S),2);
f1 = find(S>=ydb(1),1);
f2 = find(S(f1:end) <= ydb(1),1) +f1;
cut1 = round(yax(f1),3);
cut2 = round(yax(f2),3);

% cut1 = round(yax(find(S>=ydb,1)),3);
% S2 = fliplr(S);
% cutt = (find(S2>=ydb,1)); %Finish this later
% cut2 = round(yax(length(yax)-cutt));

set(handles.param4,'String','peak Z')
set(handles.output4,'String',num2str(P));
set(handles.param5,'String','minZ')
set(handles.param6,'String','maxZ')
set(handles.param8,'String','FWHM Z')
set(handles.output5,'String',num2str(cut1));
set(handles.output6,'String',num2str(cut2));
set(handles.output8,'String',num2str(cut2-cut1));


if ~isempty(handles.xlims3.String)
    xlim(str2num(handles.xlims3.String));
else
    xlim(str2num(handles.zR.String));
end
if ~isempty(handles.ylims3.String)
    ylim(str2num(handles.ylims3.String));
end


% --- Executes on button press in hold_box.
function hold_box_Callback(hObject, eventdata, handles)
% hObject    handle to hold_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hold_box


% --- Executes on button press in reset_button.
function reset_button_Callback(hObject, eventdata, handles) %@015
% hObject    handle to reset_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes2)
hold off
%plot(0);
cla reset;
axes(handles.axes3)
hold off
%plot(0);
cla reset
axes(handles.axes4)
hold off
%plot(0)
cla reset


% --- Executes during object creation, after setting all properties.
function text6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in modify_button.
function modify_button_Callback(jObject, eventdata, handles) %@016
if handles.use_chop.Value == 0
    Xfilt = evalin('base','Xfilt');
    param = evalin('base','param');
    ax = evalin('base','ax');
else
    Xfilt = evalin('base','X_c');
    param = evalin('base','param');
    ax = evalin('base','ax_c');
end
X = Xfilt;
X = circshift(X,str2double(handles.tshift.String),4);
X = circshift(X,str2double(handles.dshift.String),3);
wc1 = str2double(handles.fast_cut1.String);
wc2 = str2double(handles.fast_cut2.String);
bb = str2num(handles.baseb.String);

if handles.bb_win.Value == 1 && isempty(handles.baseb.String)
    errordlg('Enter any number great than 0 into BB freq (temp fix) to allow for windowed basenbanding')
end

if length(str2num(handles.baseb.String)) == 1
    if bb(1) > 0
        if handles.bb_win.Value == 0
            X = baseband2(X,str2double(handles.baseb.String),param.daq.HFdaq.fs_MHz,wc1,wc2);
        else
            win_n = str2double(handles.bb_win_num.String); %number of windows
            win = size(X,1)/win_n; %number of points per window
            freqs = str2num(handles.bbvar.String);
            if length(freqs) == 1
                freqs2 = ones(win_n);
                freqs = freqs2.*freqs;
            end
            if win_n > size(X,1)
                errordlg('Number of windows too large')
                return
            elseif length(freqs) ~= win_n
                errordlg('Number of baseband frequencies need to match number of windows')
                return
            end
            
            if mod(win,win_n) ~= 0
                win = floor(win);
                win_ex = round(win_n*mod(win,win_n));
                for i = 1:(win_n-1)
                    X((i-1)*win+1:(i*win),:,:,:) = baseband2(X((i-1)*win+1:(i*win),:,:,:),freqs(i),param.daq.HFdaq.fs_MHz,wc1,wc2);
                end
                i = i+1;
                X((i-1)*win+1:end,:,:,:) = baseband2(X((i-1)*win+1:end,:,:,:),freqs(end),param.daq.HFdaq.fs_MHz,wc1,wc2);
            else
                for i = 1:win_n
                    X = baseband2(X((i-1)*win+1:(i*win),:,:,:),freqs(i),param.daq.HFdaq.fs_MHz,wc1,wc2);
                end
            end
            
        end
        
        b = waitbar(0);
        if handles.signed_env.Value == 1
            S = sign(imag(X));
            if length(size(squeeze(S))) < 4
                for i = 1:size(S,4)
                    S = squeeze(S);
                    S2(:,1,:,i) = medfilt2(S(:,:,i),[5 5]);
                    S2(:,1,:,i) = medfilt2(S(:,:,i),[3,3]);
                end
                S = S2;
            else
                for i = 1:size(S,4)
                    S(:,:,:,i) = medfilt3(S(:,:,:,i),[5 5 5]);
                    S(:,:,:,i) = medfilt3(S(:,:,:,i),[3 3 3]);
                end
            end
            dims = size(Xfilt);
            b = waitbar(0);
            for i = 1:dims(1)
                for j = 1:dims(2)
                    Xfilt(i,j,:,:) = envelope(squeeze(real(Xfilt(i,j,:,:))));
                end
                waitbar(i/dims(1),b,'Basebanding');
            end
            if handles.bbdb.Value == 1
                Xfilt = 20*log10(Xfilt./max(max(max(max(abs(Xfilt))))));
            end
            if length(size(squeeze(Xfilt))) < 4
                Xf = real(squeeze(Xfilt(:,1,:,:)));
                for i = 1:size(Xfilt,4)
                    Xf2(:,:,i) = imboxfilt(Xf(:,:,i),[5 5]);
                end
                Xfilt = permute(Xf2,[1 4 2 3]);
            else
                for i = 1:size(Xfilt,4)
                    Xfilt(:,:,:,i) = imboxfilt3(real(Xfilt(:,:,:,i)),[5 5 5]);
                end
            end
            
            X = S.*real(Xfilt);
            
            %   X = S.*envelope(real(X));
        end
        delete(b)
    end
    
    if handles.invertbox.Value == 1
        X = X*(-1);
    end
    
    
    if handles.use_chop.Value == 0
        assignin('base','Xfilt',X)
    else
        assignin('base','X_c',X)
    end
    
elseif length(str2num(handles.baseb.String)) == 3
    axes(handles.axes4)
    h = hotcoldDB;
    cfreq = str2num(handles.baseb.String);
    rfreq = cfreq(1):cfreq(3):cfreq(2);
    
    for n  = 1:length(rfreq)
        R = Xfilt;
        X = baseband2(R,rfreq(n),param.daq.HFdaq.fs_MHz,wc1,wc2);
        
        if handles.signed_env.Value == 1
            S = sign(imag(X));
            dims = size(R);
            b = waitbar(0);
            for i = 1:dims(1)
                for j = 1:dims(2)
                    R(i,j,:,:) = envelope(squeeze(real(R(i,j,:,:))));
                end
                waitbar(i/dims(1),b,'Basebanding');
            end
            if handles.bbdb.Value == 1
                R = 20*log10(R./max(max(max(max(R)))));
            end
            X = S.*abs(R);
            %   X = S.*envelope(real(X));
        end
        
        xR = str2num(handles.xR.String);
        if length(xR) == 1
            xR = [xR xR];
        end
        yR = str2num(handles.yR.String);
        if length(yR) == 1
            yR = [yR yR];
        end
        zR = str2num(handles.zR.String);
        if length(zR) == 1
            zR = [zR zR];
        end
        tR = str2num(handles.tR.String);
        if length(tR) == 1
            tR = [tR tR];
        end
        aeR = str2num(handles.aeR.String);
        dims = size(X);
        %[~,ax] = make_axes(param,dims,[1 2],1);
        q.x = 1:dims(1);
        q.y = 1:dims(2);
        q.z = 1:dims(3);
        
        xInd = q.x(find(ax.x >= xR(1)):find(ax.x >= xR(2)));
        yInd = q.y(find(ax.y >= yR(1)):find(ax.y >= yR(2)));
        zInd = q.z(find(ax.depth >= zR(1)):find(ax.depth >= zR(2)));
        
        if length(size(X)) == 3
            Y = squeeze(X(xInd,yInd,zInd));
        else
            q.t = 1:dims(4);
            tInd = q.t(find(ax.stime >= tR(1)):find(ax.stime >= tR(2)));
            Y = squeeze(X(xInd,yInd,zInd,tInd));
        end
        
        if length(size(Y)) > 2
            errordlg('Too many dimensions; check ranges')
            return
        end
        %         if handles.med_box.Value == 1
        %             Y = medfilt2(Y,[3 3]);
        %         end
        delete(b)
        axes(handles.axes4)
        imagesc(ax.x(xInd),ax.depth(zInd),real(Y'))
        if ~isempty(aeR)
            caxis(aeR);
        end
        colormap(h)
        
        title(['wc = ' num2str(rfreq(n))])
        text(mean(ax.x(xInd)),ax.depth(zInd(12)),['wc = ' num2str(rfreq(n))],'Color','white');
    end
end





%
% else
%     Xfilt = evalin('base','X_c');
%     param = evalin('base','param');
%     X = Xfilt;
%     X = circshift(X,str2double(handles.tshift.String),4);
%     X = circshift(X,str2double(handles.dshift.String),3);
%     if str2double(handles.baseb.String) > 0
%          X = baseband2(X,str2double(handles.baseb.String),param.daq.HFdaq.fs_MHz);
%
%         %X = baseband_russ3(Xfilt,param.daq.HFdaq.fs_MHz,str2double(handles.baseb.String));
%         % p = squeeze(X(:,1,:,22));
%        %  figure; plot(real(ifft(p)))
%         if handles.signed_env.Value == 1
%             S = sign(imag(X));
%             dims = size(Xfilt);
%             b = waitbar(0);
%             for i = 1:dims(1)
%                 for j = 1:dims(2)
%                     Xfilt(i,j,:,:) = envelope(real(squeeze(Xfilt(i,j,:,:))));
%                    % X(i,j,:,:) = squeeze(S(i,j,:,:)).*envelope(squeeze(real(X(i,j,:,:))));
%                 end
%                 waitbar(i/dims(1),b,'Finalzing');
%             end
%             if handles.bbdb.Value == 1
%                 Xfilt = 20*log10(Xfilt./max(max(max(max(Xfilt)))));
%             end
%            X = S.*abs(Xfilt);
%         %  X = Xfilt;
%             delete(b)
%         end
%     end
%     if handles.invertbox.Value == 1
%         X = X*(-1);
%     end
%     assignin('base','X_c',X)
% end



function tshift_Callback(hObject, eventdata, handles)
% hObject    handle to tshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tshift as text
%        str2double(get(hObject,'String')) returns contents of tshift as a double


% --- Executes during object creation, after setting all properties.
function tshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function output1_Callback(hObject, eventdata, handles)
% hObject    handle to output1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output1 as text
%        str2double(get(hObject,'String')) returns contents of output1 as a double


% --- Executes during object creation, after setting all properties.
function output1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function output2_Callback(hObject, eventdata, handles)
% hObject    handle to output2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output2 as text
%        str2double(get(hObject,'String')) returns contents of output2 as a double


% --- Executes during object creation, after setting all properties.
function output2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function output3_Callback(hObject, eventdata, handles)
% hObject    handle to output3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output3 as text
%        str2double(get(hObject,'String')) returns contents of output3 as a double


% --- Executes during object creation, after setting all properties.
function output3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in IJ_butts.
function IJ_butts_Callback(hObject, eventdata, handles)
% hObject    handle to IJ_butts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SendToImageJ(handles,handles.overlay4d.Value);


function overlay4d_Callback(hObject, eventdata, handles) %@017
Overlay4D(handles);


% --- Executes on button press in AE_4dbox.
function AE_4dbox_Callback(hObject, eventdata, handles)
% hObject    handle to AE_4dbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of AE_4dbox


% --- Executes on button press in PE_4dbox.
function PE_4dbox_Callback(hObject, eventdata, handles)
% hObject    handle to PE_4dbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PE_4dbox

% --- Executes on button press in overlay.
function overlay_Callback(hObject, eventdata, handles) %@018
if handles.use_chop.Value == 1
    Xfilt = evalin('base','X_c');
    ax = evalin('base','ax_c');
else
    Xfilt = evalin('base','Xfilt');
    ax = evalin('base','ax');
end
pex = evalin('base','pex');
param = evalin('base','param');
PEData = evalin('base','PEdata');
a = str2double(handles.alph.String);
ax.pe = pex.depth;

xR = str2num(handles.xR.String);
if length(xR) == 1
    xR = [xR xR];
end
yR = str2num(handles.yR.String);
if length(yR) == 1
    yR = [yR yR];
end
zR = str2num(handles.zR.String);
if length(zR) == 1
    zR = [zR zR];
end
tR = str2num(handles.tR.String);
if length(tR) == 1
    tR = [tR tR];
end
aeR = str2num(handles.aeR.String);
peR = str2num(handles.peR.String);
dims = size(Xfilt);
%[~,ax] = make_axes(param,dims,[1 2],1);
q.x = 1:dims(1);
q.y = 1:dims(2);
q.z = 1:dims(3);
q.pe = 1:size(PEData,3);

xInd = q.x(find(ax.x >= xR(1)):find(ax.x >= xR(2)));
yInd = q.y(find(ax.y >= yR(1)):find(ax.y >= yR(2)));
zInd = q.z(find(ax.depth >= zR(1)):find(ax.depth >= zR(2)));
peInd = q.pe(find(ax.pe >= zR(1)):find(ax.pe >= zR(2)));
%peInd = zInd*2;

if handles.hotcold.Value == 1
    h = hotcoldDB;
elseif handles.graybox.Value == 1
    h = 'gray';
else
    h = 'hot';
end

if length(size(Xfilt)) == 3
    Y = squeeze(Xfilt(xInd,yInd,zInd));
else
    q.t = 1:dims(4);
    tInd = q.t(find(ax.stime >= tR(1)):find(ax.stime >= tR(2)));
    Y = squeeze(Xfilt(xInd,yInd,zInd,tInd));
end
Y = circshift(Y,str2double(handles.dshift.String),2);

if length(size(PEData)) == 3
    P = squeeze(PEData(1:size(PEData,1),yInd,peInd));
else
    
    P = squeeze(PEData(1:size(PEData,1),yInd,peInd,tInd));
end



% if size(Y) ~= size(P)
%     errordlg('PE and AE datasets must be same size')
%     return
% end

if length(size(Y)) > 2
    errordlg('Too many dimensions; check ranges')
    return
end
% if handles.med_box.Value == 1
%     Y = medfilt2(Y,[3 3]);
% end

[x, y] = meshgrid(1:size(Y,2),1:size(Y,1));
[xq, yq] = meshgrid(linspace(1,size(Y,2),size(P,2)),linspace(1,size(Y,1),size(P,1)));
Y2 = interp2(x,y,Y,xq,yq);
G = imfuse(Y2',P','ColorChannels',[1 0 2]);
H = im2double(G);
if handles.use_ext_fig.Value == 0
    axes(handles.axes1) % Plots AE
    if handles.plotbox1.Value == 1
        imagesc(ax.stime(tInd),ax.x(xInd),(Y));
        colormap(h)
        if ~isempty(aeR)
            caxis(aeR)
        end
        handles.axes1.XLabel.String = 'Time (ms)';
        handles.axes1.YLabel.String = 'Lateral (mm)';
    end
    if handles.plotbox1.Value == 2
        imagesc(ax.x(xInd),ax.y(yInd),(Y'))
        colormap(h)
        if ~isempty(aeR)
            caxis(aeR)
        end
        handles.axes1.XLabel.String = 'Lateral (mm)';
        handles.axes1.YLabel.String = 'Elevational (mm)';
    end
    if handles.plotbox1.Value == 3
        imagesc(ax.x(xInd),ax.depth(zInd),(Y'))
        colormap(h)
        if ~isempty(aeR)
            caxis(aeR)
        end
        handles.axes1.XLabel.String = 'Lateral (mm)';
        handles.axes1.YLabel.String = 'Depth (mm)';
    end
    if handles.plotbox1.Value == 4
        imagesc(ax.stime(tInd),ax.y(yInd),(Y))
        colormap(h)
        if ~isempty(aeR)
            caxis(aeR)
        end
        handles.axes1.XLabel.String = 'Time (ms)';
        handles.axes1.YLabel.String = 'Elevational (mm)';
    end
    if handles.plotbox1.Value == 5
        imagesc(ax.y(yInd),ax.depth(zInd),(Y'))
        colormap(h)
        if ~isempty(aeR)
            caxis(aeR)
        end
        handles.axes1.XLabel.String = 'Elevational (mm)';
        handles.axes1.YLabel.String = 'Depth (mm)';
    end
    if handles.plotbox1.Value == 6
        imagesc(ax.stime(tInd),ax.depth(zInd),Y)
        colormap(h)
        if ~isempty(aeR)
            caxis(aeR)
        end
        handles.axes1.XLabel.String = 'Time (ms)';
        handles.axes1.YLabel.String = 'Depth (mm)';
    end
    axes(handles.axes3) %Plots Pulse Echo
    if handles.plotbox1.Value == 1
        imagesc(ax.stime(tInd),ax.x(xInd),(P));
        colormap(gca,'gray')
        if ~isempty(peR)
            caxis(peR)
        end
        handles.axes3.XLabel.String = 'Time (ms)';
        handles.axes3.YLabel.String = 'Lateral (mm)';
    end
    if handles.plotbox1.Value == 2
        imagesc(ax.x(xInd),ax.y(yInd),(P'))
        colormap(gca,'gray')
        if ~isempty(peR)
            caxis(peR)
        end
        handles.axes3.XLabel.String = 'Lateral (mm)';
        handles.axes3.YLabel.String = 'Elevational (mm)';
    end
    if handles.plotbox1.Value == 3
        imagesc(ax.x(xInd),ax.pe(peInd),(P'))
        colormap(gca,'gray')
        if ~isempty(peR)
            caxis(peR)
        end
        handles.axes3.XLabel.String = 'Lateral (mm)';
        handles.axes3.YLabel.String = 'Depth (mm)';
    end
    if handles.plotbox1.Value == 4
        imagesc(ax.stime(tInd),ax.y(yInd),(P))
        colormap(gca,'gray')
        if ~isempty(peR)
            caxis(peR)
        end
        handles.axes3.XLabel.String = 'Time (ms)';
        handles.axes3.YLabel.String = 'Elevational (mm)';
    end
    if handles.plotbox1.Value == 5
        imagesc(ax.y(yInd),ax.pe(peInd),(P'))
        colormap(gca,'gray')
        if ~isempty(peR)
            caxis(peR)
        end
        handles.axes3.XLabel.String = 'Elevational (mm)';
        handles.axes3.YLabel.String = 'Depth (mm)';
    end
    if handles.plotbox1.Value == 6
        imagesc(ax.stime(tInd),ax.pe(peInd),P)
        colormap(gca,'gray')
        if ~isempty(peR)
            caxis(peR)
        end
        handles.axes3.XLabel.String = 'Time (ms)';
        handles.axes3.YLabel.String = 'Depth (mm)';
    end
    
    
    axes(handles.axes4)
    % H = 20*log10(H./max(max(max((H)))));
    gmax1 = max(max(G(:,:,1)));
    gmax2 = max(max(G(:,:,3)));
    for i = 1:size(G,1)
        for j = 1:size(G,2)
            
            if G(i,j,1) < gmax1*.5
                G(i,j,1) = 0;
            else
                G(i,j,1) = G(i,j,1);
            end
            if G(i,j,3) < gmax2*.5
                G(i,j,3) = 0;
            else
                G(i,j,3) = G(i,j,3);
            end
            
        end
    end
    
    imagesc(ax.x(xInd),ax.depth(zInd),G);
    handles.axes4.XLabel.String = 'Lateral (mm)';
    handles.axes4.YLabel.String = 'Depth (mm)';
    caxis(gca,[0.5 1])
    %    if ~isempty(aeR)
    %        caxis(aeR)
    %    end
    
    %    hold off
    %    plot(0)
    %
    %    yim = image(Y');
    %    pim = image(P');
    %
    %    imagesc(P')
    %       set(handles.axes4.Children,'AlphaData',0.5);
    %     hold on
    %    imagesc(Y');
    %    set(handles.axes4.Children,'AlphaData',0.5);
    %   % set(handles.axes4.Children,'AlphaData',0.5);
    h=5;
    
    
    
    
    
    
else
    figure(2)
    imagesc(G)
    %colormap(gca,'hot')
    if ~isempty(aeR)
        caxis(aeR)
    end
end




%
% axes(handles.axes2)
% imagesc(squeeze(X))
% axes(handles.axes3)
% imagesc(squeeze(PEData))






function ylims3_Callback(hObject, eventdata, handles)
% hObject    handle to ylims3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ylims3 as text
%        str2double(get(hObject,'String')) returns contents of ylims3 as a double


% --- Executes during object creation, after setting all properties.
function ylims3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ylims3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xlims3_Callback(hObject, eventdata, handles)
% hObject    handle to xlims3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xlims3 as text
%        str2double(get(hObject,'String')) returns contents of xlims3 as a double


% --- Executes during object creation, after setting all properties.
function xlims3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xlims3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ylims4_Callback(hObject, eventdata, handles)
% hObject    handle to ylims4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ylims4 as text
%        str2double(get(hObject,'String')) returns contents of ylims4 as a double


% --- Executes during object creation, after setting all properties.
function ylims4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ylims4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xlims4_Callback(hObject, eventdata, handles)
% hObject    handle to xlims4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xlims4 as text
%        str2double(get(hObject,'String')) returns contents of xlims4 as a double


% --- Executes during object creation, after setting all properties.
function xlims4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xlims4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function output6_Callback(hObject, eventdata, handles)
% hObject    handle to output6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output6 as text
%        str2double(get(hObject,'String')) returns contents of output6 as a double


% --- Executes during object creation, after setting all properties.
function output6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function output5_Callback(hObject, eventdata, handles)
% hObject    handle to output5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output5 as text
%        str2double(get(hObject,'String')) returns contents of output5 as a double


% --- Executes during object creation, after setting all properties.
function output5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function output4_Callback(hObject, eventdata, handles)
% hObject    handle to output4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output4 as text
%        str2double(get(hObject,'String')) returns contents of output4 as a double


% --- Executes during object creation, after setting all properties.
function output4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function baseb_Callback(hObject, eventdata, handles)
% hObject    handle to baseb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of baseb as text
%        str2double(get(hObject,'String')) returns contents of baseb as a double


% --- Executes during object creation, after setting all properties.
function baseb_CreateFcn(hObject, eventdata, handles)
% hObject    handle to baseb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function invert_Callback(hObject, eventdata, handles)
% hObject    handle to invert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of invert as text
%        str2double(get(hObject,'String')) returns contents of invert as a double


% --- Executes during object creation, after setting all properties.
function invert_CreateFcn(hObject, eventdata, handles)
% hObject    handle to invert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in noise_button.
function noise_button_Callback(hObject, eventdata, handles) %@019
param = evalin('base','param');
if handles.use_chop.Value == 1
    Xfilt = evalin('base','X_c');
    ax = evalin('base','ax_c');
else
    Xfilt = evalin('base','Xfilt');
    ax = evalin('base','ax');
end
Xfilt = Xfilt.*1000000./str2double(handles.hfgain.String);
dims = size(Xfilt);
ave = mean(mean(mean(mean(abs(Xfilt)))));
b = waitbar(0,'Calculating Noise');
for i = 1:dims(2)
    for j = 1:dims(3)
        for k = dims(1)
            dev1(i,j,k) = std(Xfilt(k,i,j,:));
        end
        waitbar((1-i/dims(2))+j/dims(3),b,'Calculating Noise')
    end
end
delete(b)
dev = mean(mean(mean(dev1)));
if isempty(handles.output5.String)
    pres = 1;
else
    pres = str2double(handles.output5.String);
end
slope = str2double(handles.output4.String);
LF = evalin('base','LF');
thresh = (2*dev+ave); %/pres
detect = thresh/slope; %maybe multi by pres
set(handles.param1,'String','mean')
set(handles.param2,'String','std')
set(handles.param3,'String','thresh')
set(handles.output1,'String',num2str(ave));
set(handles.output2,'String',num2str(dev));
set(handles.output3,'String',num2str(thresh));
set(handles.param4,'String','sense')
set(handles.param5,'String','pres')
set(handles.param6,'String','detect')
set(handles.output6,'String',num2str(detect));


% --- Executes on button press in onemhz.
function onemhz_Callback(hObject, eventdata, handles)
if handles.onemhz.Value == 1;
    set(handles.trans,'Value',4)
    trans_Callback(hObject, eventdata, handles)
end
% hObject    handle to onemhz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of onemhz



function lfgain_Callback(hObject, eventdata, handles)
% hObject    handle to lfgain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lfgain as text
%        str2double(get(hObject,'String')) returns contents of lfgain as a double


% --- Executes during object creation, after setting all properties.
function lfgain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lfgain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hfgain_Callback(hObject, eventdata, handles)
% hObject    handle to hfgain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hfgain as text
%        str2double(get(hObject,'String')) returns contents of hfgain as a double


% --- Executes during object creation, after setting all properties.
function hfgain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hfgain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in med_box.
function med_box_Callback(hObject, eventdata, handles)
% hObject    handle to med_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of med_box



function reals_Callback(hObject, eventdata, handles)
% hObject    handle to reals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reals as text
%        str2double(get(hObject,'String')) returns contents of reals as a double


% --- Executes during object creation, after setting all properties.
function reals_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in invertbox.
function invertbox_Callback(hObject, eventdata, handles)
% hObject    handle to invertbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of invertbox


% --- Executes on button press in realbox.
function realbox_Callback(hObject, eventdata, handles)
% hObject    handle to realbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of realbox


% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in create_pe.
function create_pe_Callback(hObject, eventdata, handles) %@020
% hObject    handle to create_pe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB

if handles.onemhz.Value == 1
    [file, path] = uigetfile(fullfile(pwd,'*_info.dat')); %Gets file location
    param = read_ucsdi_info([path file]); %Gets scan parameters
    cd(path);
    
    ax.HFfreq = linspace(0,param.daq.HFdaq.fs_MHz,param.daq.HFdaq.pts); %Creates fast frequency axis
    ax.LFfreq = linspace(0,param.daq.HFdaq.pulseRepRate_Hz,param.daq.HFdaq.NoBurstTriggers); %creates slow frequency axis
    
    [~, HF1] = full_signal([path file],param,1); %Gets the raw data
    
    X = w_ae_filt2(param,HF1,1,0,[str2double(handles.fast_cut1.String) str2double(handles.fast_cut2.String)]); %Filters in fast time; 0 is match, 1 uses cutoffs
    
    
    %%%%%%%%%%%%%%%%%%%% Gets new axes for z and t %%%%%%%%%%%%%%%%%%
    b = waitbar(0,'Matrix Conversion');
    waitbar(0,b,'Enveloping and converting to 4D matrix')
    HF = zeros(size(X,1),size(X,2),size(X{1},1),size(X{1},2));
    
    if param.velmex.SlowAxis == 'X'
        LatStep = param.velmex.YNStep;
        EleStep = param.velmex.XNStep;
    else
        LatStep = param.velmex.XNStep;
        EleStep = param.velmex.YNStep;
    end
    
    for i = 1:LatStep
        for j = 1:EleStep
            %HF(i,j,:,:) = envelope(real(X{i,j})); %Converts cell array to double
            HF(i,j,:,:) = X{i,j}; %Converts cell array to double
        end
        waitbar(i/LatStep,b,'Converting to 4D matrix');
    end
    delete(b)
    
    %     if param.velmex.XNStep ~= 1 && param.velmex.YNStep ~= 1
    %     %    if param.velmex.SlowAxis == 'X'
    %             for i = 1:size(HF,1)
    %                 if mod(i,2) == 0
    %                     HF(i,:,:,:) = fliplr(HF(i,:,:,:));
    %                 end
    %             end
    %      %   end
    %     end
    
    if param.velmex.XNStep ~= 1 && param.velmex.YNStep ~= 1
        if param.velmex.SlowAxis == 'X'
            for i = 1:size(HF,1)
                if mod(i,2) == 0
                    HF(i,:,:,:) = fliplr(HF(i,:,:,:));
                end
            end
        else
            for i = 1:size(HF,2)
                if mod(i,2) == 0
                    HF(:,i,:,:) = fliplr(HF(:,i,:,:));
                end
            end
        end
    end
    
    
    
    % PEdata = HF(:,:,:,1:2);
    PEdata = HF;
    [~, pex] = make_axes(param,size(HF));
    pex.depth = linspace(0,1.48*param.daq.HFdaq.pts/param.daq.HFdaq.fs_MHz/2,size(PEdata,3));
    if handles.save_4d.Value == 1
        clearvars -except file path param pex PEdata
        f = file(1:end-4);
        f2 = [f '_4d_PE.mat'];
        fprintf('Saving 4D file...')
        save(f2);
        fprintf('Done\n')
    end
    
    %FOR PROCESSING BSQDATA
elseif handles.bsq.Value == 1
    [f, p]  = uigetfile(fullfile(pwd,'*_PEParm.mat'));
    f_root = f(1:(end-11));
    load([p f]);
    cd(p)
    PEbsqFile = fullfile(p,[f_root,'_PEImgData.bsq']);
    fid = fopen(PEbsqFile,'rb');
    nScanPt = 1:bScanParm.nScanPt;
    n = fread(fid,1,'int32');
    dsize = fread(fid,[1,n],'int32');
    nOffset = (n+1)*4;
    fclose(fid);
    
    PEImage = multibandread(PEbsqFile,[dsize(1:2),prod(dsize(3:end))],...
        'single',nOffset,'bsq','ieee-le',{'Band','Direct',nScanPt});
    if bScanParm.ylen > 1
        PEImage = reshape(PEImage,[dsize(1),dsize(2),bScanParm.XSteps,bScanParm.YSteps]);
    end
    pedata = squeeze(mean(PEImage,3));
    wavlen = PData.Lambda;
    x_mm = ((1:dsize(2))*wavlen*PData(1).PDelta(1));
    x_mm = x_mm - mean(x_mm);
    z_mm = (1:dsize(1))*wavlen*PData(1).PDelta(3);
    pex.x = x_mm;
    pex.depth = z_mm;
    if bScanParm.ylen == 0
        pex.y = 1;
    else
        pex.y = linspace(-bScanParm.ylen/2,bScanParm.ylen/2,bScanParm.YSteps);
    end
    pedata = permute(pedata,[2 3 1]);
    if handles.save_4d.Value == 1
        clearvars -except f bScanParm pedata PData Trans TW TX pex TXArray
        f = f(1:end-10);
        file02 = [f '_vbf_4d_PE.mat'];
        fprintf('Saving 4D file...')
        save(file02);
        fprintf('Done\n')
    end
    
    %FOR PROCESSING CHIRP DATA VERASONICS DATA
elseif handles.onemhz.Value == 0 && handles.bsq.Value == 0
    [f, p]  = uigetfile(fullfile(pwd,'*_PEParm.mat'));
    load([p f]);
    cd(p)
    
    %  b = waitbar(0,'Filtering');
    US = TW.Wvfm2Wy;
    sz = size(PEMatrix);
    %  PEM = zeros(sz(1)
    %     PEMatrix = permute(PEMatrix,[2 3 1]);
    %     PEMatrix = reshape(PEMatrix,[bScanParm.XSteps bScanParm.YSteps sz(1) sz(2)]);
    
    bScanParm.depth = Rcv(1).endDepth*PData.Lambda;
    
    
    % Do some filterings here!!! bad in mind that PE is downsampled from
    % 250MHz to something around 10 MHz
    if handles.match_box.Value == 0
        N = size(PEMatrix,1);
        wc1 = str2double(handles.fast_cut1.String);
        wc2 = str2double(handles.fast_cut2.String);
        fast_axis = linspace(0,round(Rcv(1).decimSampleRate),size(PEMatrix,1));
        f1 = find(fast_axis >= wc1,1);
        f2 = find(fast_axis >= wc2,1);
        H(1:f1-1) = 0;
        H(f1:f2) = 1;
        H(f2+1:N-f2-1) = 0;
        H(N-f2:N-f1) = 1;
        H(N-f1+1:N) = 0;
        h = ifft(H);
        h2 = circshift(h,round(N/2)).*hamming(N)';
        H2 = fft(h2);
        
        Hn = round(length(find(H))/2);
        ham = hamming(Hn);
        fh1 = find(fast_axis >= mean([f1,f2]),1);
        fh2 = N-fh1;
        
        if length(sz) > 2
            for i = 1:sz(3)
                Q{i} = fft(squeeze(PEMatrix(:,:,i)));
            end
            
            
            
            
            %Insert FFT filter here
            for i = 1:sz(3)
                for j = 1:sz(2)
                    PE{i}(:,j) = Q{i}(:,j).*H2';
                end
                %  waitbar(i/size(PEMatrix,3),b,'Filtering')
                multiWaitbar('Filtering',i/size(PEMatrix,3));
            end
            for i = 1:sz(3)
                %  pe2{i} = real(flipud(ifft(PE{i})));
                pe{i} = circshift(real(ifft(PE{i})),round(size(PE{i},1)/2));
            end
        else
            Q = fft(PEMatrix);
            for j = 1:sz(2)
                PE(:,j) = Q(:,j).*H2';
            end
            pe = circshift(real(ifft(PE)),round(size(PE,1)/2));
        end
        
        
    else
        FsAE = round(Rcv(1).decimSampleRate);
        FsUS = bScanParm.vsx_fs;
        t_delay = length(US)*(1/FsUS);
        % US = interp1(linspace(0,100,length(US)),US,linspace(0,100,length(US)*2))';
        if FsUS~=FsAE
            RefPulse       = resample(US,FsAE,FsUS);
        end
        RefPulse = RefPulse/(sum(abs(RefPulse)));
        RefPulse = flipud(conj(RefPulse));
        hwin = hamming(length(RefPulse));
        %RefPulse = hwin.*RefPulse;
        
        %d_fix = round(sz(1)*.75); %0.75 is a scaling factor chosen to stretch the decompressed image back to regular size
        
        for i = 1:sz(3)
            for j = 1:sz(2)
                Q{i}(:,j) = conv(squeeze(PEMatrix(:,j,i)),RefPulse);
            end
            %   Q{i} = Q{i}(1:round(size(Q{i},1)*.75),:); %scaling new data range
            %  waitbar(i/sz(3),b,'Filtering');
            multiWaitbar('Filtering',i/sz(3));
        end
        
        %         for i = 1:sz(3)
        %             for j = 1:sz(2)
        %
        %                 pe{i}(:,j) = interp1(linspace(0,10,size(Q{1},1)),Q{i}(:,j),linspace(0,10,sz(1)));
        %
        %             end
        %             %   waitbar(i/sz(3),b,'Compressing Depth Axis');
        %             multiWaitbar('Compressing Depth Axis',i/sz(3));
        %         end
        pe = Q;
        
    end
    
    if length(sz)>2
        for i = 1:sz(3)
            pdata(:,:,i) = pe{i};
        end
    else
        pdata=pe;
    end
    sz = size(pdata);
    
    pdata = permute(pdata,[3 1 2]);
    pdata = reshape(pdata,[bScanParm.XSteps bScanParm.YSteps sz(1) sz(2)]);
    fstele = find(TX.VDASApod,1)-1;
    if length(TX.VDASApod) == size(pdata,4)
    for i = 1:length(TX.VDASApod)
        if TX.VDASApod(i) == 1        
            pedata(:,:,:,i-fstele) = pdata(:,:,:,i);              
        end
    end
    else
        pedata  = pdata;
    end
    %
    %     PEformed = mean(abs(pedata),4);
    %     figure;
    %     for i = 1:64
    %         imagesc(squeeze(pedata(:,1,:,i))')
    %         drawnow;
    %         pause(.1)
    %         title(num2str(i));
    %     end
    
    
    
    % Create axes data
    if bScanParm.xlen ~= 1
        pex.x = linspace(-bScanParm.xlen/2,bScanParm.xlen/2,size(pedata,1));
    else
        pex.x = 1;
    end
    if bScanParm.ylen ~= 1
        pex.y = linspace(-bScanParm.ylen/2,bScanParm.ylen/2,size(pedata,2));
    else
        pex.y = 1;
    end
    pex.element = Trans.ElementPos(:,1)*PData.Lambda;
    pex.element = pex.element.*TX.Apod';
    pex.element = pex.element(pex.element~=0);
    cent = round(length(pex.element)/2);
    pex.depth = linspace(0,bScanParm.depth,size(pedata,3));
    pex.stime = linspace(0,bScanParm.Duration,size(pedata,4));
    apods = find(TX.VDASApod);
    delay.x = Trans.ElementPos(:,1)*PData.Lambda;
    %for i = 1:
    
    
    
    
    
    for i = 1:size(pedata,4) %Element
        for j = 1:size(pedata,1) %Lateral Position
            %bf.r(:,i,j) = sqrt((pex.element(i)-pex.x(j))^2+pex.depth.^2);
            pex.theta(i,j) = abs((atan((pex.x(j)-pex.element(i))/(TX.focus*PData.Lambda))));
            %pex.theta(i,j) = abs((atan((pex.x(j)-pex.element(round(end/2)))/(TX.focus*PData.Lambda))));
            if pex.x(j)-pex.element(i) < 0
                pex.theta(i,j) = pex.theta(i,j)*-1;
            end
        end
    end
    
    
    
    % LETS DO SOME BEAM FORMING!!!
    for m = 1:size(pedata,2) %Elevational
        for a = 1:length(pex.element) %Element
            %a = 1+c(k);
            
            %Decimation and Interpolation to reshape data
            if handles.match_box.Value == 1
                D = real(squeeze(pedata(:,m,:,a)));
            else
                %  D = squeeze(PEMatrix(:,a+32,:))'; %This is not currently using filtered data
                D = real(squeeze(pedata(:,m,:,a)));
            end
            [grdx1, grdz1] = meshgrid(1:size(D,1),1:size(D,2));
            [grdx2, grdz2] = meshgrid(linspace(1,size(D,1),PData(1).Size(2)),linspace(1,size(D,2),PData(1).Size(1)));
            Dint = interp2(grdx1, grdz1,D',grdx2, grdz2);
            
            clear D;
            D = Dint';
            
            C = interp1(1:size(pedata,1),pex.theta(cent,:),linspace(1,size(pedata,1),size(D,1)));  % Adjust where Center is here
            %   C = pex.theta(cent,:);
            wavelen = PData.Lambda;
            %pex.x = (1:size(pedata,1))*wavlen*PData(1).PDelta(1).*(PData.Size(2)/size(pedata,1));
            pex.x = (1:size(D,1))*wavelen*PData(1).PDelta(1);
            pex.x = pex.x - mean(pex.x);
            %  pex.depth = (1:size(pedata,3))*wavlen*PData(1).PDelta(3).*(PData.Size(1)/size(pedata,3));
            pex.depth = (1:size(D,2))*wavelen*PData(1).PDelta(3); % + delay.x(1,1);
            
            t_delay = length(US)*(1/bScanParm.vsx_fs); %adjusts depth based on length of pulse
            depth2 = pex.depth; % - t_delay; %Correcting
            Qx = pex.depth.*sin(C');
            Qz = pex.depth.*cos(C');
            x2 = linspace(min(min(Qx)),max(max(Qx)),size(Qx,1));
            foc = find(pex.depth >= TX.focus*PData.Lambda,1);
            Qxind = zeros(size(Qz,1),size(Qz,2));
            Qzind = Qxind;
            Qweight = Qzind;
            for i = 1:size(Qz,1) %Lateral
                for j = 1:size(Qz,2) %Depth
                    Qxind(i,j) = find(pex.x >= Qx(i,j),1); %or use x2
                    %   Qzind(i,j) = find(pex.depth+(TXArray(i).Delay(a)*PData.Lambda) >= Qz(i,j),1); %Experimental
                    if max(depth2) >= Qz(i,j)
                        Qzind(i,j) = find(depth2 >= Qz(i,j),1); %Pair with pex.depth = pex.depth-max(TX.Delay)*PData.Lambda
                    else
                        Qzind(i,j) = length(depth2);
                    end
                    if Qzind(i,j) <= foc
                        Qweight(i,j) = 1+Qzind(i,j)/foc;
                    else
                        Qweight(i,j) = 1+foc/Qzind(i,j);
                    end
                end
                multiWaitbar(['Element ' num2str(a)],i/size(Qx,1));
            end
            
            
            Qfin = zeros(size(Qx));
            Qnum = Qfin;
            for i = 1:size(Qz,1) % X
                for j = 1:size(Qz,2) %Depth
                    Qfin(Qxind(i,j),Qzind(i,j),m) = Qfin(Qxind(i,j),Qzind(i,j),m)+D(i,j,m).*Qweight(i,j);
                    Qnum(Qxind(i,j),Qzind(i,j),m) = Qnum(Qxind(i,j),Qzind(i,j),m) +1;
                end
            end
            Qtot = Qfin./Qnum;
            Qtot(isnan(Qtot)) = 0;
            Qtot = abs(Qtot)';
            BF(:,:,m,a) = Qtot;
            
            if a < length(pex.element)
                multiWaitbar(['Element ' num2str(a)],'Relabel',['Element ' num2str(a+1)]);
            else
                multiWaitbar(['Element ' num2str(a)],'Close');
            end
            
            % figure; imagesc(x2,pex.depth,Qtot)
            % ylim([35 55]);
            %  waitbar(a/length(pex.element)*(m/size(pedata,2)),b,'Beamforming');
            multiWaitbar('Beamforming Elements',a/length(pex.element));
        end
        multiWaitbar('Beamforming Elevation',m/size(pedata,2));
    end
    multiWaitbar('CLOSEALL');
    bf = mean(BF,4);
    %pex.x = x2;
    
    %figure; imagesc(x2,pex.depth,bf2)
    
    if length(size(bf)) < 3
        % bf2 = medfilt2(bf,[3,1]);
        bf2 = permute(bf,[2 3 1]);
        pex.y = 1;
    else
        % bf2 = medfilt3(bf,[3,1,1]);
        bf2 = permute(bf,[2 3 1]);
    end
    clear pedata
    pedata = bf2;
    
    assignin('base','PEdata',pedata);
    %pex.depth = pex.depth-max(TX.Delay)*PData.Lambda;
    %delete(b)
    
    if handles.save_4d.Value == 1
        clearvars -except f bScanParm pedata PData Trans TW TX pex Rcv TXArray
        f = f(1:end-10);
        file02 = [f '_4d_PE.mat'];
        fprintf('Saving 4D file...')
        save(file02);
        fprintf('Done\n')
    end
end



% --- Executes on button press in loadpe. @021
function loadpe_Callback(hObject, eventdata, handles)
[f,  p] = uigetfile(fullfile(pwd,'*4d_PE.mat'));
cd(p)
fprintf('Loading 4D Dataset...')
load([p f]);
%[~, ax] = make_axes(param,size(Xfilt),[1 2],1);
if handles.onemhz.Value == 1
    set(handles.fname,'String',file);
    fprintf('Done\n')
    assignin('base','PEdata',PEdata);
    assignin('base','param',param);
    pex.depth = pex.depth - 2.6;
    ax = pex;
    assignin('base','pex',ax);
    set(handles.active_pe,'String',num2str(size(PEdata)))
    %assignin('base','PEparam',PE);
    %set(handles.LF_chan,'String',num2str(size(LF,2)));
    set(handles.tms,'String',num2str([ax.stime(1) ax.stime(end)]));
    set(handles.tsamp,'String',num2str([1 length(ax.stime)]));
    set(handles.xmm,'String',num2str([ax.x(1) ax.x(end)]));
    set(handles.xsamp,'String',num2str([1 length(ax.x)]));
    set(handles.ymm,'String',num2str([ax.y(1) ax.y(end)]));
    set(handles.ysamp,'String',num2str([1 length(ax.y)]));
    set(handles.zmm,'String',num2str([round(ax.depth(1)) round(ax.depth(end))]));
    set(handles.zsamp,'String',num2str([1 length(ax.depth)]));
    
    if handles.reset_axes.Value == 1
        set(handles.xR,'String', num2str([ax.x(1) ax.x(end)]));
        set(handles.yR,'String', num2str([ax.y(1) ax.y(end)]));
        set(handles.zR,'String', num2str([ax.depth(1) floor(ax.depth(end))]));
        set(handles.tR,'String', num2str([ax.stime(1) ax.stime(end)]));
    end
    
else
    set(handles.fname,'String',file02);
    fprintf('Done\n')
    assignin('base','PEdata',pedata);
    assignin('base','pex',pex);
    assignin('base','param',bScanParm);
    assignin('base','Trans',Trans)
    assignin('base','PData',PData)
    assignin('base','TW',TW)
    assignin('base','TX',TX)
    if exist('TXArray')
    assignin('base','TXArray',TXArray);
    end
    set(handles.tP,'String','0');
end


% --- Executes on button press in usepe. @022
function usepe_Callback(hObject, eventdata, handles)
X = evalin('base','PEdata');
ax = evalin('base','pex');

if length(size(X)) < 4
    X(:,:,:,2) = X(:,:,:);
    ax.stime = [0 1];
end
assignin('base','Xfilt',X);
assignin('base','ax',ax);
set(handles.active_xfilt,'String','PE')
if ~isempty(handles.hfchans.String)
    set(handles.active_chan,'String',handles.hfchans.String);
end


% --- Executes on button press in realize.
function realize_Callback(hObject, eventdata, handles)
if handles.use_chop.Value == 0
    X = evalin('base','Xfilt');
    param = evalin('base','param');
    X = real(X);
    %Xcom = imag(X);
    % assignin('base','X_complex',Xcom)
    assignin('base','Xfilt',X)
else
    X = evalin('base','X_c');
    param = evalin('base','param');
    X = real(X);
    %Xcom = imag(X);
    %assignin('base','X_complex',Xcom)
    assignin('base','X_c',X)
end


% --- Executes on button press in signed_env.
function signed_env_Callback(hObject, eventdata, handles)
% hObject    handle to signed_env (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of signed_env


% --- Executes on button press in figfolder.
function figfolder_Callback(hObject, eventdata, handles)
path = uigetdir;
set(handles.savefolder,'String',path);

% --- Executes on button press in savefig. @023
function savefig_Callback(hObject, eventdata, handles)
fname = handles.savefigname.String;
if ~isempty(handles.fignum.String)
    num = str2double(handles.fignum.String);
    if handles.borderless.Value == 1
        figure(num);
        a = gca;
        a.Visible = 'off';
    end
    %     if fname(end-2:end) == 'png'
    %         fname = [fname ' -transparent'];
    %         figure(num)
    %         set(gca,'Color','none');
    %     else
    %     end
    if ~isempty(handles.savefolder.String)
        pname = handles.savefolder.String;
        sname = [pname '\' fname];
    else
        sname = fname;
    end
    
    figure(num)
    set(gca,'Color','none')
    %  eval([ 'export_fig ' sname])
    export_fig(sname,'-transparent');
else
    if ~isempty(handles.savefolder.String)
        pname = handles.savefolder.String;
        sname = [pname '\' fname];
    else
        sname = fname;
    end
    export_fig(sname)
    
    
end





function savefolder_Callback(hObject, eventdata, handles)
% hObject    handle to savefolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of savefolder as text
%        str2double(get(hObject,'String')) returns contents of savefolder as a double


% --- Executes during object creation, after setting all properties.
function savefolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savefolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function savefigname_Callback(hObject, eventdata, handles)
% hObject    handle to savefigname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of savefigname as text
%        str2double(get(hObject,'String')) returns contents of savefigname as a double


% --- Executes during object creation, after setting all properties.
function savefigname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to savefigname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fignum_Callback(hObject, eventdata, handles)
% hObject    handle to fignum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fignum as text
%        str2double(get(hObject,'String')) returns contents of fignum as a double


% --- Executes during object creation, after setting all properties.
function fignum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fignum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bbdb.
function bbdb_Callback(hObject, eventdata, handles)
% hObject    handle to bbdb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bbdb


% --- Executes on button press in large_box.
function large_box_Callback(hObject, eventdata, handles)
% hObject    handle to large_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of large_box


% --- Executes on button press in tope. @024
function tope_Callback(hObject, eventdata, handles)
if handles.use_chop.Value == 1
    X = evalin('base','X_c');
    ax = evalin('base','ax_c');
else
    X = evalin('base','Xfilt');
    ax = evalin('base','ax');
end
assignin('base','PEdata',X);
assignin('base','pex',ax);
set(handles.active_pe,'String',num2str(size(X)))



function alph_Callback(hObject, eventdata, handles)
% hObject    handle to alph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alph as text
%        str2double(get(hObject,'String')) returns contents of alph as a double


% --- Executes during object creation, after setting all properties.
function alph_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function dshift_Callback(hObject, eventdata, handles)
% hObject    handle to dshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dshift as text
%        str2double(get(hObject,'String')) returns contents of dshift as a double


% --- Executes during object creation, after setting all properties.
function dshift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dshift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function peR_Callback(hObject, eventdata, handles)
% hObject    handle to peR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of peR as text
%        str2double(get(hObject,'String')) returns contents of peR as a double


% --- Executes during object creation, after setting all properties.
function peR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to peR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot4. @025
function plot4_Callback(hObject, eventdata, handles)
param = evalin('base','param');
if handles.use_chop.Value == 1
    Xfilt = evalin('base','X_c');
    ax = evalin('base','ax_c');
else
    Xfilt = evalin('base','Xfilt');
    ax = evalin('base','ax');
end

Xfilt = real(Xfilt);
%set(handles.axes1,'ButtonDownFcn',@Plot4OnClickXZ)

%plot(handles.axes3,ax.x,'ButtonDownFcn',@Plot4OnClickXZ)
xP = str2double(handles.xP.String);
yP = str2double(handles.yP.String);
zP = str2double(handles.zP.String);
tP = str2double(handles.tP.String);

xR = str2num(handles.xR.String);
if length(xR) == 1
    xR = [xR xR xR];
else
    xR = [xR(1) xR(2) xP];
end
yR = str2num(handles.yR.String);
if length(yR) == 1
    yR = [yR yR yR];
else
    yR = [yR(1) yR(2) yP];
end
zR = str2num(handles.zR.String);
if length(zR) == 1
    zR = [zR zR zR];
else
    zR = [zR(1) zR(2) zP];
end
tR = str2num(handles.tR.String);
if length(tR) == 1
    tR = [tR tR tR];
else
    tR = [tR(1) tR(2) tP];
end
aeR = str2num(handles.aeR.String);

%ax.current = [xR(3) yR(3) zR(3) tR(3)];


dims = size(Xfilt);
%[~,ax] = make_axes(param,dims,[1 2],1);
q.x = 1:dims(1);
q.y = 1:dims(2);
q.z = 1:dims(3);
if length(size(Xfilt)) > 3
    q.t = 1:dims(4);
end
xInd = q.x(find(ax.x >= xR(1)):find(ax.x >= xR(2)));
yInd = q.y(find(ax.y >= yR(1)):find(ax.y >= yR(2)));
zInd = q.z(find(ax.depth >= zR(1)):find(ax.depth >= zR(2)));
if length(size(Xfilt)) > 3
    tInd = q.t(find(ax.stime >= tR(1)):find(ax.stime >= tR(2)));
end
if length(xR) < 3 || length(yR) < 3 || length(zR) <3 || length(tR) <3
    errordlg('All 4 dimensions need 3 values ([range1, range2, point])')
    return
end
px = q.x(find(ax.x >=xR(3),1));
py = q.y(find(ax.y >=yR(3),1));
pz = q.z(find(ax.depth >=zR(3),1));
if length(size(Xfilt)) > 3
    pt = q.t(find(ax.stime >=tR(3),1));
end
if length(size(Xfilt)) < 4
    Yxy = squeeze(Xfilt(xInd,yInd,pz));
    Yxz = squeeze(Xfilt(xInd,py,zInd));
    Yyz = squeeze(Xfilt(px,yInd,zInd));
    Yzt = squeeze(Xfilt(px,py,zInd));
else
    Yxy = squeeze(Xfilt(xInd,yInd,pz,pt));
    Yxz = squeeze(Xfilt(xInd,py,zInd,pt));
    Yyz = squeeze(Xfilt(px,yInd,zInd,pt));
    Yzt = squeeze(Xfilt(px,py,zInd,tInd));
end

if handles.hotcold.Value == 1
    h = hotcoldDB;
elseif handles.graybox.Value == 1
    h = 'gray';
else
    h = 'hot';
end

if length(size(Yxy)) > 2
    errordlg('Too many dimensions; check ranges')
    return
end
% if handles.med_box.Value == 1
%     Yxy = medfilt2(Yxt,[3 3]);
%     Yxz = medfilt2(Yxz,[3 3]);
%     Yyz = medfilt2(Yyz,[3 3]);
%     Yzt = medfilt2(Yzt,[3 3]);
% end
if handles.use_ext_fig.Value == 0
    axes(handles.axes2)
    if size(Yxy,2) == 1
        Yxy = Yxy;
    end
    imagesc(ax.x(xInd),ax.y(yInd),(Yxy'),'ButtonDownFcn',{@Plot4OnClickXY,handles})
    colormap(gca,h)
    if ~isempty(aeR)
        caxis(aeR)
    end
    handles.axes2.XLabel.String = 'Lateral (mm)';
    handles.axes2.YLabel.String = 'Elevational (mm)';
    
    axes(handles.axes1)
    if size(Yxz,2) == 1
        Yxz = Yxz';
    end
    imagesc(ax.x(xInd),ax.depth(zInd),(Yxz'),'ButtonDownFcn',{@Plot4OnClickXZ,handles})
    colormap(gca,h)
    if ~isempty(aeR)
        caxis(aeR)
    end
    handles.axes1.XLabel.String = 'Lateral (mm)';
    handles.axes1.YLabel.String = 'Depth (mm)';
    
    axes(handles.axes3)
    if size(Yyz,2) == 1
        Yyz = Yyz';
    end
    imagesc(ax.y(yInd),ax.depth(zInd),(Yyz'),'ButtonDownFcn',{@Plot4OnClickYZ,handles})
    colormap(gca,h)
    if ~isempty(aeR)
        caxis(aeR)
    end
    handles.axes3.XLabel.String = 'Elevational (mm)';
    handles.axes3.YLabel.String = 'Depth (mm)';
    
    axes(handles.axes4)
    if size(Yzt,2) == 1
        Yzt = Yzt';
    end
    if length(size(Xfilt)) > 3
        imagesc(ax.stime(tInd),ax.depth(zInd),Yzt,'ButtonDownFcn',{@Plot4OnClickTZ,handles})
    else
        imagesc(1,ax.depth(zInd),Yzt,'ButtonDownFcn',{@Plot4OnClickTZ,handles})
    end
    colormap(gca,h)
    if ~isempty(aeR)
        caxis(aeR)
    end
    handles.axes4.XLabel.String = 'Time (ms)';
    handles.axes4.YLabel.String = 'Depth (mm)';
    
else
    figure(51)
    imshow(Yxz')
    colormap(gca,h)
    if ~isempty(aeR)
        caxis(aeR)
    end
    xlabel('Lateral (mm)');
    ylabel('Depth (mm)');
    figure(52)
    imshow(Yyz)
    colormap(gca,h)
    if ~isempty(aeR)
        caxis(aeR)
    end
    xlabel('Elevational (mm)');
    ylabel('Depth (mm)');
    figure(53)
    imshow(Yxy')
    colormap(gca,h)
    if ~isempty(aeR)
        caxis(aeR)
    end
    xlabel('Lateral (mm)');
    ylabel('Elevational (mm)');
    figure(54)
    imshow(Yzt)
    colormap(gca,h)
    if ~isempty(aeR)
        caxis(aeR)
    end
    xlim([tInd(1) tInd(end)])
    ylim([zInd(1) zInd(end)])
    xlabel('Time (ms)');
    ylabel('Depth (mm)');
end


% --- Executes on button press in hotcold.
function hotcold_Callback(hObject, eventdata, handles)
% hObject    handle to hotcold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hotcold



function hfchans_Callback(hObject, eventdata, handles)
% hObject    handle to hfchans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hfchans as text
%        str2double(get(hObject,'String')) returns contents of hfchans as a double


% --- Executes during object creation, after setting all properties.
function hfchans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hfchans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MergeHF. @026
function MergeHF_Callback(hObject, eventdata, handles)
n = str2double(handles.numhf.String);
b = waitbar(0,'Loading');
for i = 1:n
    [file, path] = uigetfile(fullfile(pwd,'*4d_data.mat'));
    
    load([path file]);
    %assignin('base',['X' num2str(i)],Xfilt)
    X{i} = Xfilt;
    waitbar(i/n,b,'Merging 4d datasets')
end
delete(b)
assignin('base','Xmerged',X);
assignin('base','fpath',[path file]);
assignin('base','param',param);
assignin('base','ax',ax);
assignin('base','LF',LF);
%assignin('base','PEparam',PE);
set(handles.LF_chan,'String',num2str(size(LF,2)));
set(handles.tms,'String',num2str([ax.stime(1) ax.stime(end)]));
set(handles.tsamp,'String',num2str([1 length(ax.stime)]));
set(handles.xmm,'String',num2str([ax.x(1) ax.x(end)]));
set(handles.xsamp,'String',num2str([1 length(ax.x)]));
set(handles.ymm,'String',num2str([ax.y(1) ax.y(end)]));
set(handles.ysamp,'String',num2str([1 length(ax.y)]));
set(handles.zmm,'String',num2str([ax.depth(1) round(ax.depth(end))]));
set(handles.zsamp,'String',num2str([1 length(ax.depth)]));

if handles.reset_axes.Value == 1
    set(handles.xR,'String', num2str([ax.x(1) ax.x(end)]));
    set(handles.yR,'String', num2str([ax.y(1) ax.y(end)]));
    set(handles.zR,'String', num2str([ax.depth(1) floor(ax.depth(end))]));
    set(handles.tR,'String', num2str([ax.stime(1) ax.stime(end)]));
    
end




function numhf_Callback(hObject, eventdata, handles)
% hObject    handle to numhf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numhf as text
%        str2double(get(hObject,'String')) returns contents of numhf as a double


% --- Executes during object creation, after setting all properties.
function numhf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numhf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function channel_Callback(hObject, eventdata, handles)
% hObject    handle to channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of channel as text
%        str2double(get(hObject,'String')) returns contents of channel as a double


% --- Executes during object creation, after setting all properties.
function channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in usechan. @027
function usechan_Callback(hObject, eventdata, handles)
m = str2double(handles.channel.String);
Xmerged = evalin('base','Xmerged');
if m > length(Xmerged)
    errordlg('Your selected value exceeds the number of channels');
end
X = Xmerged{m};
assignin('base','Xfilt',X)
set(handles.active_xfilt,'String','AE');
if ~isempty(handles.hfchans.String)
    set(handles.active_chan,'String',handles.hfchans.String);
end



function framerate_Callback(hObject, eventdata, handles)
% hObject    handle to framerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of framerate as text
%        str2double(get(hObject,'String')) returns contents of framerate as a double


% --- Executes during object creation, after setting all properties.
function framerate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to framerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in graybox.
function graybox_Callback(hObject, eventdata, handles)
% hObject    handle to graybox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of graybox


% --- Executes on button press in lfmovie. @028
function lfmovie_Callback(hObject, eventdata, handles) %Low Frequency Movie
% hObject    handle to lfmovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
param = evalin('base','param');
LF = evalin('base','LF');
LF = LF(:,str2double(handles.LF_chan.String));
LF = LF/str2double(handles.lfgain.String)*1000;
tR = str2num(handles.tR.String);
if length(tR) == 1
    LF_butt_Callback(hObject, eventdata, handles)
else
    q.t = 1:param.daq.LFdaq.pts;
    ax.lf = linspace(0,param.daq.duration_ms,param.daq.LFdaq.pts);
    tInd = q.t(find(ax.lf >= tR(1),1):find(ax.lf >= tR(2),1));
    b1 = min(LF); b2 = max(LF);
    b = param.daq.LFdaq.pts/param.daq.HFdaq.NoBurstTriggers;
    
    
    if handles.LF_FFT.Value == 1
        lf = fft(LF);
        x = linspace(0,param.daq.LFdaq.fs_Hz,length(lf));
        if handles.use_ext_fig.Value == 1
            figure(33)
            plot(x,abs(lf))
        else
            axes(handles.axes3)
            plot(x,abs(lf))
        end
        if ~isempty(handles.xlims3.String)
            xlim(str2num(handles.xlims3.String));
        end
        if ~isempty(handles.ylims3.String)
            ylim(str2num(handles.ylims3.String));
        end
        
    end
    x = linspace(0,param.daq.HFdaq.duration_ms,length(LF));
    if handles.save_fig.Value == 1
        v = VideoWriter('LFdata.avi');
        v.FrameRate = str2double(handles.framerate.String);
        open(v)
    end
    for i = round(tInd(1):b:tInd(end))
        if handles.use_ext_fig.Value == 1
            figure(3)
            if ~isempty(handles.xlims.String)
                xlim(str2num(handles.xlims.String));
            end
            if ~isempty(handles.ylims.String)
                ylim(str2num(handles.ylims.String));
            end
            plot(x,LF,'k')
            hold on
            plot(x(i),LF(i),'ro','MarkerFaceColor','r')
            hold off
            set(gca,'Color','none')
            if handles.save_fig.Value == 1
                frame = getframe;
                writeVideo(v,frame);
            else
                drawnow
            end
        else
            axes(handles.axes2)
            if ~isempty(handles.xlims.String)
                xlim(str2num(handles.xlims.String));
            end
            if ~isempty(handles.ylims.String)
                ylim(str2num(handles.ylims.String));
            end
            plot(x,LF,'k')
            hold on
            plot(x(i),LF(i),'ro','MarkerFaceColor','r')
            hold off
            drawnow;
        end
        ylabel('mA');
        xlabel('ms');
        
    end
    if handles.save_fig.Value == 1
        close(v)
    end
end


% --- Executes on button press in pushbutton36.
function pushbutton36_Callback(hObject, eventdata, handles)
if handles.use_chop.Value == 1
    X = evalin('base','X_c');
    Xmerged = evalin('base','Xmerged');
    Xmerged{str2double(handles.channel.String)} = X;
    assignin('base','Xmerged',Xmerged);
else
    X = evalin('base','Xfilt');
    Xmerged = evalin('base','Xmerged');
    Xmerged{str2double(handles.channel.String)} = X;
    assignin('base','Xmerged',Xmerged);
end


% --- Executes on button press in showlf.
function showlf_Callback(hObject, eventdata, handles)
% hObject    handle to showlf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showlf

function Plot4OnClickXZ(hObject,eventdata,handles) %@029
pt = get(gca,'currentpoint');


set(handles.xP,'String',pt(1,1))
set(handles.zP,'String',pt(1,2))
if handles.show_fft.Value == 1
    param = evalin('base','param');
    if handles.use_chop.Value == 1
        ax = evalin('base','ax_c');
    else
        ax = evalin('base','ax');
    end
    yloc = find(ax.depth >= pt(1,2),1);
    xloc = find(ax.x >= pt(1,1),1);
    [X,xaxis,yaxis] = plot_fft(param,hObject.CData,1);
    axes(handles.axes2)
    plot(yaxis,X(1:length(yaxis),xloc));
    xlabel('Depth Freq MHz')
    axes(handles.axes3)
    plot(xaxis,X(yloc,1:length(xaxis)))
    xlabel('Lateral Spatial Freq kHz')
else
    plot4_Callback(hObject, eventdata, handles);
end

function Plot4OnClickYZ(hObject,eventdata,handles)
pt = get(gca,'currentpoint');
set(handles.yP,'String',pt(1,1))
set(handles.zP,'String',pt(1,2))
plot4_Callback(hObject, eventdata, handles);


function Plot4OnClickXY(hObject,eventdata,handles)
pt = get(gca,'currentpoint');
set(handles.xP,'String',pt(1,1))
set(handles.yP,'String',pt(1,2))
plot4_Callback(hObject, eventdata, handles);


function Plot4OnClickTZ(hObject,eventdata,handles)
pt = get(gca,'currentpoint');
set(handles.tP,'String',pt(1,1))
set(handles.zP,'String',pt(1,2))
if handles.show_fft.Value == 1
    param = evalin('base','param');
    if handles.use_chop.Value == 1
        ax = evalin('base','ax_c');
    else
        ax = evalin('base','ax');
    end
    yloc = find(ax.depth >= pt(1,2),1);
    xloc = find(ax.stime >= pt(1,1),1);
    [X,xaxis,yaxis] = plot_fft(param,hObject.CData,4);
    axes(handles.axes2)
    plot(yaxis,X(1:length(yaxis),xloc));
    xlabel('Depth Freq MHz')
    axes(handles.axes3)
    plot(xaxis,X(yloc,1:length(xaxis)))
    xlabel('SlowTime Freq Hz')
else

plot4_Callback(hObject, eventdata, handles);
end



function xP_Callback(hObject, eventdata, handles)
% hObject    handle to xP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xP as text
%        str2double(get(hObject,'String')) returns contents of xP as a double


% --- Executes during object creation, after setting all properties.
function xP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yP_Callback(hObject, eventdata, handles)
% hObject    handle to yP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yP as text
%        str2double(get(hObject,'String')) returns contents of yP as a double


% --- Executes during object creation, after setting all properties.
function yP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zP_Callback(hObject, eventdata, handles)
% hObject    handle to zP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zP as text
%        str2double(get(hObject,'String')) returns contents of zP as a double


% --- Executes during object creation, after setting all properties.
function zP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tP_Callback(hObject, eventdata, handles)
% hObject    handle to tP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tP as text
%        str2double(get(hObject,'String')) returns contents of tP as a double


% --- Executes during object creation, after setting all properties.
function tP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in all_movie.
function all_movie_Callback(hObject, eventdata, handles)
% hObject    handle to all_movie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of all_movie



function med_x_Callback(hObject, eventdata, handles)
% hObject    handle to med_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of med_x as text
%        str2double(get(hObject,'String')) returns contents of med_x as a double


% --- Executes during object creation, after setting all properties.
function med_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to med_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function med_y_Callback(hObject, eventdata, handles)
% hObject    handle to med_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of med_y as text
%        str2double(get(hObject,'String')) returns contents of med_y as a double


% --- Executes during object creation, after setting all properties.
function med_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to med_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function med_z_Callback(hObject, eventdata, handles)
% hObject    handle to med_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of med_z as text
%        str2double(get(hObject,'String')) returns contents of med_z as a double


% --- Executes during object creation, after setting all properties.
function med_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to med_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bb_win_num_Callback(hObject, eventdata, handles)
% hObject    handle to bb_win_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bb_win_num as text
%        str2double(get(hObject,'String')) returns contents of bb_win_num as a double


% --- Executes during object creation, after setting all properties.
function bb_win_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bb_win_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bbvar_Callback(hObject, eventdata, handles)
% hObject    handle to bbvar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bbvar as text
%        str2double(get(hObject,'String')) returns contents of bbvar as a double


% --- Executes during object creation, after setting all properties.
function bbvar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bbvar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in bb_win.
function bb_win_Callback(hObject, eventdata, handles)
% hObject    handle to bb_win (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bb_win


% --- Executes on button press in bsq.
function bsq_Callback(hObject, eventdata, handles)
% hObject    handle to bsq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bsq


% --- Executes on slider movement. @030
function tshifter_Callback(hObject, eventdata, handles)
if handles.use_chop.Value
    ax = evalin('base','ax_c');
else
    ax = evalin('base','ax');
end
%set(hObject,'Value',str2double(handles.tP.String));
val = get(hObject,'Value');
if val < ax.stime(1) || val > ax.stime(end)
    errordlg('T value outside of range');
    return
end
if length(ax.stime) > 1
    inc = ax.stime(2)-ax.stime(1);
    handles.tshifter.SliderStep = [inc/30, inc/3];
end
val = round(val,1);
val = num2str(val);
set(handles.tP,'String',val);
set(hObject,'Value',str2double(val));
set(hObject,'Max',ax.stime(end));
set(hObject,'Min',ax.stime(1));
plot_ae_Callback(hObject, eventdata, handles)
% hObject    handle to tshifter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function tshifter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tshifter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement. @031
function zshifter_Callback(hObject, eventdata, handles)
if handles.use_chop.Value
    ax = evalin('base','ax_c');
else
    ax = evalin('base','ax');
end
%set(hObject,'Value',str2double(handles.zP.String));
val = get(hObject,'Value');
if val < ax.depth(1) || val > ax.depth(end)
    errordlg('Z value outside of range');
    return
end
if length(ax.depth) > 1
    inc = ax.depth(2)-ax.depth(1);
    handles.zshifter.SliderStep = [inc/30, inc/3];
end
val = round(val,2);
val = num2str(val);
set(handles.zP,'String',val);
set(hObject,'Value',str2double(val));
set(hObject,'Max',ax.depth(end));
set(hObject,'Min',ax.depth(1));
plot_ae_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function zshifter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zshifter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement. @032
function yshifter_Callback(hObject, eventdata, handles)
if handles.use_chop.Value
    ax = evalin('base','ax_c');
else
    ax = evalin('base','ax');
end
val = get(hObject,'Value');
if val < ax.y(1) || val > ax.y(end)
    errordlg('Y value outside of range');
    return
end
if length(ax.y) > 1
    inc = ax.y(2)-ax.y(1);
    handles.y.SliderStep = [inc/10, inc];
end
set(handles.yP,'String',num2str(val));
set(hObject,'Max',ax.y(end));
set(hObject,'Min',ax.y(1));
plot_ae_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function yshifter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yshifter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in stitch_gui. @033
function stitch_gui_Callback(hObject, eventdata, handles)
Stitch;
% hObject    handle to stitch_gui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in borderless.
function borderless_Callback(hObject, eventdata, handles)
% hObject    handle to borderless (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of borderless



function output7_Callback(hObject, eventdata, handles)
% hObject    handle to output7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output7 as text
%        str2double(get(hObject,'String')) returns contents of output7 as a double


% --- Executes during object creation, after setting all properties.
function output7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function output8_Callback(hObject, eventdata, handles)
% hObject    handle to output8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of output8 as text
%        str2double(get(hObject,'String')) returns contents of output8 as a double


% --- Executes during object creation, after setting all properties.
function output8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tc.
function tc_Callback(hObject, eventdata, handles)
% hObject    handle to tc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tc


% --- Executes on button press in TC_param.
function TC_param_Callback(hObject, eventdata, handles)
tc_param;
% hObject    handle to TC_param (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in squarify_box.
function squarify_box_Callback(hObject, eventdata, handles)
% hObject    handle to squarify_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of squarify_box


% --- Executes on selection change in trans.
function trans_Callback(hObject, eventdata, handles)
% hObject    handle to trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.trans.Value == 1
    set(handles.fast_cut1,'String',1);
    set(handles.fast_cut2,'String',4);
elseif handles.trans.Value == 2
    set(handles.fast_cut1,'String',2);
    set(handles.fast_cut2,'String',4);
elseif handles.trans.Value == 3
    set(handles.fast_cut1,'String',0.3);
    set(handles.fast_cut2,'String',0.9);
elseif handles.trans.Value == 4
    set(handles.fast_cut1,'String',0.6);
    set(handles.fast_cut2,'String',1.4);
end

% Hints: contents = cellstr(get(hObject,'String')) returns trans contents as cell array
%        contents{get(hObject,'Value')} returns selected item from trans


% --- Executes during object creation, after setting all properties.
function trans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in newaescan.
function newaescan_Callback(hObject, eventdata, handles)
% hObject    handle to newaescan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of newaescan


% --- Executes on button press in ft_on.
function ft_on_Callback(hObject, eventdata, handles)
% hObject    handle to ft_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ft_on


% --- Executes on button press in st_on.
function st_on_Callback(hObject, eventdata, handles)
% hObject    handle to st_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of st_on


% --- Executes on button press in show_fft.
function show_fft_Callback(hObject, eventdata, handles)
% hObject    handle to show_fft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of show_fft


% --- Executes on button press in noleadfield.
function noleadfield_Callback(hObject, eventdata, handles)
% hObject    handle to noleadfield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noleadfield


% --- Executes on button press in leadfield.
function leadfield_Callback(hObject, eventdata, handles)
% hObject    handle to leadfield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[f,  p] = uigetfile(fullfile(pwd,'*.mat'));
cd(p)
fprintf('Loading 3D Pressure Field...')
load([p f]);
fprintf('Done\n')
assignin('base','L',LEAD_FIELD_Matrix);


% --- Executes on button press in pressurefield.
function pressurefield_Callback(hObject, eventdata, handles)
% hObject    handle to pressurefield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[f,  p] = uigetfile(fullfile(pwd,'*.mat'));
cd(p)
fprintf('Loading 3D Pressure Field...')
load([p f]);
fprintf('Done\n')
assignin('base','Pressure',Pressure);


% --- Executes on button press in J_mag.
function J_mag_Callback(hObject, eventdata, handles)
% hObject    handle to J_mag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of J_mag


% --- Executes on button press in calcJ.
function calcJ_Callback(hObject, eventdata, handles)
% hObject    handle to calcJ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% if ~exist('base','Pressure')
%     errordlg('Need to load pressure matrix first');
% else
Pressure = evalin('base','Pressure');
% end
if handles.use_chop == 0
V = evalin('base','Xfilt');
ax = evalin('base','ax');
else
V = evalin('base','X_c');
ax = evalin('base','ax_c');
end
sz = size(V);

t = str2double(get(handles.tP,'String'));
if length(ax.stime) > 1
    t0 = find(ax.stime >= t,1);
    V = V(:,:,:,t0);
end
    

V = squeeze(V);
Pressure = squeeze(Pressure);
%Temp Code
%Pressure = Pressure(:,:,:,467:734);
%
%Calculates J matrix
dec = floor(size(Pressure,4)/size(Pressure,3));
for i = 1:size(Pressure,3)
    Pressure(:,:,i,1) = Pressure(:,:,i,i*dec);
    multiWaitbar('Compressing Pressure Matrix',i/size(Pressure,3));
end
Pressure = Pressure(:,:,:,1);

if handles.J_mag.Value == 1 %taking only magnitude and not direction of pressure
    if length(size(V)) == 2
        P2 = envelope(Pressure);
    else
        for i = 1:size(Pressure,3)
            P2(:,:,i) = envelope(Pressure(:,:,i)); % Gets env of pressure
        end
    end
else
    P2 = Pressure;
end

% figure;
% for i = 1:63
%     imagesc(squeeze(Pressure(:,:,i)))
%       title(i);
%     drawnow;
%     pause(0.1)
% end
    


dims = size(P2);
in = length(P2)/length(V);

[x, y, z] = meshgrid(1:dims(2),1:dims(1),1:dims(3));
[x2, y2, z2] = meshgrid(1:(1/in):dims(2),1:(1/in):dims(1),1:(1/in):dims(3));
P2(:,:,:) = interp3(x,y,z,P2(:,:,:),x2,y2,z2);

keep = find(max(max(P2)) ~= 0);
P3 = P2(:,:,keep);


J = zeros(size(V));
if handles.noleadfield.Value == 1 % assuming little variation in lead field matrix
%     for i = 1:size(V,3)
%         if max(max(abs(P2(:,:,i)))) == 0
%             J(:,:,i) = zeros(size(J,1),size(J,2));
%         else
        J(:,:,:) = deconvlucy(V(:,:,:),P3(:,:,:)); %deconvolues pressure and AE voltage
       
%         multiWaitbar('Deconvolving',i/size(Pressure,3));
%     end
   
     assignin('base','J',J);
      multiWaitbar('CLOSEALL');
else
    J1 = deconvlucy(V,P2);
    L = evalin('base','L');
    L1 = inv(L);
    J = L1*J;
    assignin('base','J',J)
end


% --- Executes on button press in useJ.
function useJ_Callback(hObject, eventdata, handles)
% hObject    handle to useJ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
J = evalin('base','J');
if length(size(J)) == 2
    J = permute(J,[1 3 2]);
end
if handles.use_chop.Value == 1
    assignin('base','X_c',J);
else
    assignin('base','Xfilt',J);
end


% --- Executes on button press in keep_box.
function keep_box_Callback(hObject, eventdata, handles)
% hObject    handle to keep_box (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of keep_box
