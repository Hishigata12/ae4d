function varargout = ae4d(varargin)
% AE4D MATLAB code for ae4d.fig
%      AE4D, by itself, creates a new AE4D or raises the existing
%      singleton*.
%
%      H = AE4D returns the handle to a new AE4D or the handle to
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

% Last Modified by GUIDE v2.5 27-May-2018 11:34:01

% Begin initialization code - DO NOT EDIT
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
function ae4d_OpeningFcn(hObject, eventdata, handles, varargin)
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
set(handles.aeR, 'String','-9 0');

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
function plot_ae_Callback(hObject, eventdata, handles)
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

xInd = q.x(find(ax.x >= xR(1)):find(ax.x >= xR(2)));
yInd = q.y(find(ax.y >= yR(1)):find(ax.y >= yR(2)));
zInd = q.z(find(ax.depth >= zR(1)):find(ax.depth >= zR(2)));

if length(size(Xfilt)) == 3
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
if handles.med_box.Value == 1
    Y = medfilt2(Y,[3 3]);
end
if handles.use_ext_fig.Value == 0
    axes(handles.axes1)
    if handles.plotbox1.Value == 1
        imagesc(ax.stime(tInd),ax.x(xInd),(Y));
        colormap('hot')
        if ~isempty(aeR)
            caxis(aeR)
        end
        handles.axes1.XLabel.String = 'Time (ms)';
        handles.axes1.YLabel.String = 'Lateral (mm)';
    end
    if handles.plotbox1.Value == 2
        imagesc(ax.x(xInd),ax.y(yInd),(Y'))
        colormap('hot')
        if ~isempty(aeR)
            caxis(aeR)
        end
        handles.axes1.XLabel.String = 'Lateral (mm)';
        handles.axes1.YLabel.String = 'Elevational (mm)';
    end
    if handles.plotbox1.Value == 3
        imagesc(ax.x(xInd),ax.depth(zInd),(Y'))
        colormap('hot')
        if ~isempty(aeR)
            caxis(aeR)
        end
        handles.axes1.XLabel.String = 'Lateral (mm)';
        handles.axes1.YLabel.String = 'Depth (mm)';
    end
    if handles.plotbox1.Value == 4
        imagesc(ax.stime(tInd),ax.y(yInd),(Y))
        colormap('hot')
        if ~isempty(aeR)
            caxis(aeR)
        end
        handles.axes1.XLabel.String = 'Time (ms)';
        handles.axes1.YLabel.String = 'Elevational (mm)';
    end
    if handles.plotbox1.Value == 5
        imagesc(ax.y(yInd),ax.depth(zInd),(Y'))
        colormap('hot')
        if ~isempty(aeR)
            caxis(aeR)
        end
        handles.axes1.XLabel.String = 'Elevational (mm)';
        handles.axes1.YLabel.String = 'Depth (mm)';
    end
    if handles.plotbox1.Value == 6
        imagesc(ax.stime(tInd),ax.depth(zInd),Y)
        colormap('hot')
        if ~isempty(aeR)
            caxis(aeR)
        end
        handles.axes1.XLabel.String = 'Time (ms)';
        handles.axes1.YLabel.String = 'Depth (mm)';
    end
else
    figure(2)
    imshow(Y')
    colormap(gca,'hot')
    if ~isempty(aeR)
        caxis(aeR)
    end
end

assignin('base','ax',ax);
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
function load_data_Callback(hObject, eventdata, handles)
% hObject    handle to load_data (see GCBO)
[f,  p] = uigetfile(fullfile(pwd,'*4d_data.mat'));
cd(p)
fprintf('Loading 4D Dataset...')
load([p f]);
[~, ax] = make_axes(param,size(Xfilt),[1 2],1);
set(handles.fname,'String',file);
fprintf('Done\n')
assignin('base','Xfilt',Xfilt);
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



% --- Executes during object creation, after setting all properties.
function load_data_CreateFcn(hObject, eventdata, handles)
% hObject    handle to load_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in create_4d.
function create_4d_Callback(hObject, eventdata, handles)
% hObject    handle to create_4d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, path] = uigetfile(fullfile(pwd,'*_info.dat')); %Gets file location
param = read_ucsdi_info([path file]); %Gets scan parameters
[~,~,LF] = read_ucsdi_data([path file],1); %Gets input current waveform
cd(path);
if handles.match_box.Value == 1
[file2, path2] = uigetfile(fullfile(pwd,'*PEParm.mat')); %gets US pulse waveform
PE = open([path2 file2]);
US = PE.TW.Wvfm1Wy;
end
ax.HFfreq = linspace(0,param.daq.HFdaq.fs_MHz,param.daq.HFdaq.pts); %Creates fast frequency axis
ax.LFfreq = linspace(0,param.daq.HFdaq.pulseRepRate_Hz,param.daq.HFdaq.NoBurstTriggers); %creates slow frequency axis


%**************************************************************************
%Builds 4D Matrix***************~~~~~~~~~~~~~~~~**************
[~, HF1] = full_signal([path file],param,2); %Gets the raw data
if ~isempty(handles.slow_cut2.String) && ~isempty(handles.slow_cut1.String) && handles.slow_box.Value == 0
    X = w_slow_filt2(param,HF1,LF,handles.slow_box.Value,[str2double(handles.slow_cut1.String) str2double(handles.slow_cut2.String)]); %Filters in slow time 0 is match, 1 uses cutoffs
elseif handles.slow_box.Value == 1
    X = w_slow_filt2(param, HF1,LF,handles.slow_box.Value);
end
if handles.match_box.Value == 1
X = w_ae_filt2(param,X,US,1); %Filters in fast time; 0 is match, 1 uses cutoffs
end
if handles.match_box.Value == 0
  X = w_ae_filt2(param,X,1,0,[str2double(handles.fast_cut1.String) str2double(handles.fast_cut2.String)]); %Filters in fast time; 0 is match, 1 uses cutoffs  
end

%%%%%%%%%%%%%%%%%%%% Gets new axes for z and t %%%%%%%%%%%%%%%%%%
b = waitbar(0,'Matrix Conversion');
waitbar(0,b,'Enveloping and converting to 4D matrix')
HF = zeros(size(X,1),size(X,2),size(X{1},1),size(X{1},2));
for i = 1:param.velmex.XNStep
    for j = 1:param.velmex.YNStep
         %HF(i,j,:,:) = envelope(real(X{i,j})); %Converts cell array to double
          HF(i,j,:,:) = X{i,j}; %Converts cell array to double
    end
    waitbar(i/param.velmex.XNStep,b,'Converting to 4D matrix');
end

%%%%%%%%
if isempty(num2str(handles.depR.String))
    qq = [10 round(1.48*param.daq.HFdaq.pts/param.daq.HFdaq.fs_MHz-10)];
else
    qq = str2num(handles.depR.String);
end

dims = size(HF);
[M, ax] = make_axes(param,dims,qq, 12.3); %selects range for dB calculation exlcuding the 10mm around each border

% XdB = real(20*log10(real(HF)./max(max(max(max(real(HF(:,:,M.xT,:))))))));
% Xfilt = filts2D(XdB,[0 1 12],[0 2 2]);
Xfilt = HF;
delete(b)

if handles.keep.Value == 1
    assignin('base','Xfilt',Xfilt);
    assignin('base','fpath',[path file]);
    assignin('base','param',param);
    assignin('base','ax',ax);
     assignin('base','PE',PE);
    assignin('base','LF',LF);
    set(handles.fname,'String',[path file]);
end
if handles.save_4d.Value == 1
    clearvars -except Xfilt file path param ax LF PE
    f = file(1:end-4);
    f2 = [f '_4d_data.mat'];
    fprintf('Saving 4D file...')
    eval([ 'save ' f2 ' -v7.3']);
    fprintf('Done\n')
end

    


% --- Executes on button press in keep.
function keep_Callback(hObject, eventdata, handles)
% hObject    handle to keep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of keep


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
function movie_button_Callback(hObject, eventdata, handles)
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

if handles.use_ext_fig.Value == 0
    axes(handles.axes2)
else
    figure(1);
end
if handles.plotbox2.Value == 1
    if handles.save_fig.Value == 0
        for k = tInd %Mod loop
            if handles.med_box.Value == 1
                J = medfilt2((squeeze(Xfilt(xInd,yInd,zInd,k)))',[5 5]);
            else
                J = squeeze(Xfilt(xInd,yInd,zInd,k))';
            end

            if handles.use_ext_fig.Value == 0
                imagesc(ax.x(xInd),ax.depth(zInd),J) % mod plots
                h = hotcoldDB;
                colormap(h)
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
        vidwrite(Xfilt,ax.depth,ax.x,[zInd(1) zInd(end)],[xInd(1) xInd(end)],[tInd(1) tInd(end)],aeR)
    end
    
elseif handles.plotbox2.Value == 2
    if handles.save_fig.Value == 0
        for k = tInd %Mod loop
            if handles.med_box.Value == 1
                J = medfilt2((squeeze(Xfilt(xInd,yInd,zInd,k)))',[5 5]);
            else
                J = squeeze(Xfilt(xInd,yInd,zInd,k))';
            end
      
            if handles.use_ext_fig.Value == 0
                imagesc(ax.y(yInd),ax.depth(zInd),J) % mod plots
                h = hotcoldDB;
                colormap(h)
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
        vidwrite(Xfilt,ax.depth,ax.y,[zInd(1) zInd(end)],[yInd(1) yInd(end)],[tInd(1) tInd(end)],aeR)
    end
    
elseif handles.plotbox2.Value == 3
    if handles.save_fig.Value == 0
        for k = zInd %Mod loop
            if handles.med_box.Value == 1
                J = medfilt2((squeeze(Xfilt(xInd,yInd,k,tInd)))',[5 5]);
            else
                J = squeeze(Xfilt(xInd,yInd,k,tInd))';
            end
            
            if handles.use_ext_fig.Value == 0
                imagesc(ax.x(xInd),ax.y(yInd),J) % mod plots
                h = hotcoldDB;
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
        vidwrite(Xfilt,ax.y,ax.x,[yInd(1) yInd(end)],[xInd(1) xInd(end)],[zInd(1) zInd(end)],aeR)
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
function Enhance_Sig_Callback(hObject, eventdata, handles)
% hObject    handle to Enhance_Sig (see GCBO)
if handles.use_chop.Value == 0
param = evalin('base','param');
Xfilt = evalin('base','Xfilt');
m = [handles.mean_box.Value str2double(handles.mean_x.String) str2double(handles.mean_y.String) str2double(handles.mean_z.String)];
n = [handles.int_box.Value str2double(handles.int_x.String) str2double(handles.int_y.String) str2double(handles.int_z.String)];
Xfilt = filts3D(Xfilt,m,n,param);
[~,ax] = make_axes(param,size(Xfilt));
assignin('base','ax',ax);
assignin('base','Xfilt',Xfilt)
else
    param = evalin('base','param');
Xfilt = evalin('base','X_c');
ax = evalin('base','ax_c');
m = [handles.mean_box.Value str2double(handles.mean_x.String) str2double(handles.mean_y.String) str2double(handles.mean_z.String)];
n = [handles.int_box.Value str2double(handles.int_x.String) str2double(handles.int_y.String) str2double(handles.int_z.String)];
Xfilt = filts3D(Xfilt,m,n,param);

dims = size(Xfilt);
xR = [ax.x(1) ax.x(end)];
yR = [ax.y(1) ax.y(end)];
zR = [ax.depth(1) ax.depth(end)];
tR = [ax.stime(1) ax.stime(end)];
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

assignin('base','ax_c',ax);
assignin('base','X_c',Xfilt)
end
    


% --- Executes on button press in chop.
function chop_Callback(hObject, eventdata, handles)
% hObject    handle to chop (see GCBO)
Xfilt = evalin('base','Xfilt');
param = evalin('base','param');
clear X_c
xR = str2num(handles.xR.String);
if length(xR) == 1
    xR(2) = xR(1);
end
yR = str2num(handles.yR.String);
if length(yR) == 1
    yR(2) = yR(1);
end
zR = str2num(handles.zR.String);
if length(zR) == 1
    zR(2) = zR(1);
end
tR = str2num(handles.tR.String);
if length(tR) == 1
    tR(2) = tR(1);
end
aeR = str2num(handles.aeR.String);
dims = size(Xfilt);
[~,ax] = make_axes(param,dims,[1 2],1);
q.x = 1:dims(1);
q.y = 1:dims(2);
q.z = 1:dims(3);
q.t = 1:dims(4);
xInd = q.x(find(ax.x >= xR(1)):find(ax.x >= xR(2)));
yInd = q.y(find(ax.y >= yR(1)):find(ax.y >= yR(2)));
zInd = q.z(find(ax.depth >= zR(1)):find(ax.depth >= zR(2)));
tInd = q.t(find(ax.stime >= tR(1)):find(ax.stime >= tR(2)));
X = Xfilt(xInd,yInd,zInd,tInd);
if length(tInd) == 1
    X = permute(X,[1 2 3 4]);
end
%[~,ax] = make_axes(param,size(X));
dims = size(X);
if length(dims) < 4
    dims = [dims(1) dims(2) dims(3) 1];
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


assignin('base','ax_c',ax);
assignin('base','X_c',X);


% --- Executes on button press in max_box.
function max_box_Callback(hObject, eventdata, handles)
% hObject    handle to max_box (see GCBO)
if get(hObject,'Value') == 1
    if handles.use_chop.Value == 0
    ax = evalin('base','ax');
set(handles.xR,'String', num2str([ax.x(1) ax.x(end)]));
set(handles.yR,'String', num2str([ax.y(1) ax.y(end)]));
set(handles.zR,'String', num2str([ax.depth(1) floor(ax.depth(end))]));
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
function env_button_Callback(hObject, eventdata, handles)
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
function irad_Callback(hObject, eventdata, handles)
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
if length(size(HF)) == 3
    q.t = 1:size(HF,3);
    tInd =  q.t(find(ax.stime >= t(1)):find(ax.stime >= t(2)));
    HF = permute(HF,[2 1 3]);
    b = waitbar(0);
    
    for i = tInd
        R(:,:,i) = iradon(HF(:,:,i),str2double(handles.dtheta.String),'None');
        waitbar(i/tInd(end),b,'Computing inverse radon transform');
    end
    
    R = permute(R,[1 4 2 3]);
    delete(b);
    dims = size(R); %gets no dimensions of reconstructed data
    ax.depth = linspace(0,ax.depth(end),dims(3));
    ax.x = linspace(ax.x(1),ax.x(end),dims(1));
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
            waitbar((j-1)/yInd + i/yInd(end),b,'Computing inverse radon transform');
        end
    end
end

R = permute(R,[2 3 1 4]);
delete(b);
dims = size(R); %gets no dimensions of reconstructed data
ax.depth = linspace(0,ax.depth(end),dims(3));
ax.x = linspace(ax.x(1),ax.x(end),dims(1));
ax.y = linspace(ax.y(1),ax.y(end),dims(2));
    

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
function dB_Button_Callback(hObject, eventdata, handles)
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
function LF_butt_Callback(hObject, eventdata, handles)
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
function sense_button_Callback(hObject, eventdata, handles)
% hObject    handle to sense_button (see GCBO)
%datacursormode on
%dcm = datacursormode(gcf);
[x,y] = ginput(1);
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
        yR = str2num(handles.yR.String);
        if length(yR) > 1
            yR = yR(1);
        end
        sx = find(ax.x>=x,1);
        sz = find(ax.depth>=y,1);
        Sae = squeeze(Xfilt(sx,yR,sz,tInd));
        Sae2 = resample(Sae,length(Lae),length(Sae));  %This might need to be interp1
        %Sae2 = interp1(linspace(0,dims(4),dims(4)),Sae,linspace(0,dims(4),length(Lae)));
        
    end
    if handles.plotbox1.Value == 2
        zR = str2num(handles.zR.String);
        if length(zR) > 1
            zR = zR(1);
        end
        sx = find(ax.x>=x,1);
        sy = find(ax.y>=y,1);
        Sae = squeeze(Xfilt(sx,sy,zR,tInd));
        Sae2 = resample(Sae,length(Lae),length(S));
        
    end
    if handles.plotbox1.Value == 5
        xR = str2num(handles.xR.String);
        if length(xR) > 1
            xR = xR(1);
        end
        sx = find(ax.y>=x,1);
        sz = find(ax.depth>=y,1);
        Sae = squeeze(Xfilt(xR,sx,sz,tInd));
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

if handles.use_ext_fig.Value == 1
    figure(6);
    hold off;
    plot(0)
   % scatter(abs(Lae),abs(Sae2));
    scatter(Lae,Sae2,13,'r','filled')
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
     scatter(Lae,Sae2,13,'k','filled')
     ylabel('\muV')
    xlabel('mA')
    axes(handles.axes3)
    hold on
    plot(T_axis,Lnorm,'k')
    plot(T_axis,Snorm,'r')
    title(['R^2 = ' num2str(R)]);
      xlabel('ms')
    hold off
end

%m = round(mean(abs(Sae2)./abs(Lae)),2);
m = round(mean(Sae2./Lae),2);
m2 = (max(abs(Sae2))-min(abs(Sae2)))/(max(abs(Lae))-min(abs(Lae)));
ave = round(mean(abs(Sae2)),2);
dev = round(std((Sae2)),2);
if ~isempty(handles.output5.String)
    pres = str2num(handles.output5.String);
else
    pres = 1;   
end
%m3 = m/pres;
m3 = m2/pres;

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
function fwhm_button_Callback(hObject, eventdata, handles)
[x,y] = ginput(1);

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


set(handles.param1,'String','peak')
set(handles.output1,'String',num2str(P));
set(handles.param2,'String','minX')
set(handles.param3,'String','maxX')
set(handles.output2,'String',num2str(cut1));
set(handles.output3,'String',num2str(cut2));

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


f1 = find(S>=ydb(1),1);
f2 = find(S(f1:end) <= ydb(1),1) +f1;
cut1 = round(yax(f1),3);
cut2 = round(yax(f2),3);

% cut1 = round(yax(find(S>=ydb,1)),3);
% S2 = fliplr(S);
% cutt = (find(S2>=ydb,1)); %Finish this later
% cut2 = round(yax(length(yax)-cutt));

set(handles.param5,'String','minZ')
set(handles.param6,'String','maxZ')
set(handles.output5,'String',num2str(cut1));
set(handles.output6,'String',num2str(cut2));


if ~isempty(handles.xlims3.String)
    xlim(str2num(handles.xlims3.String));
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
function reset_button_Callback(hObject, eventdata, handles)
% hObject    handle to reset_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes2)
hold off
plot(0);
axes(handles.axes3)
hold off
plot(0)
axes(handles.axes4)
hold off
plot(0)


% --- Executes during object creation, after setting all properties.
function text6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in modify_button.
function modify_button_Callback(hObject, eventdata, handles)
if handles.use_chop.Value == 0
    X = evalin('base','Xfilt');
    param = evalin('base','param');
    X = circshift(X,str2double(handles.tshift.String),4);
    if str2double(handles.baseb.String) > 0
        if handles.bbdb.Value == 1
            X = 20*log10(mean(mean(mean(mean(X))))/X);
        end
        X = baseband2(X,str2double(handles.baseb.String),param.daq.HFdaq.fs_MHz);
     
        if handles.signed_env.Value == 1
            S = real(sign(X));
            dims = size(X);
            b = waitbar(0);
            for i = 1:dims(1)
                for j = 1:dims(2)                
                    X(i,j,:,:) = squeeze(S(i,j,:,:)).*envelope(squeeze(real(X(i,j,:,:))));          
                end
                waitbar(i/dims(1),b,'Basebanding');
            end
            delete(b)
         %   X = S.*envelope(real(X));
        end
    end
    if handles.invertbox.Value == 1
        X = X*(-1);
    end
    assignin('base','Xfilt',X)
    
else
    Xfilt = evalin('base','X_c');
    param = evalin('base','param');
    X = Xfilt;
    X = circshift(X,str2double(handles.tshift.String),4);
    if str2double(handles.baseb.String) > 0
      %  X = baseband2(X,str2double(handles.baseb.String),param.daq.HFdaq.fs_MHz);
       
        X = baseband_russ3(X,param.daq.HFdaq.fs_MHz,str2double(handles.baseb.String));
        if handles.signed_env.Value == 1
            S = sign(real(X));
            dims = size(Xfilt);
            b = waitbar(0);
            for i = 1:dims(1)
                for j = 1:dims(2)   
                    Xfilt(i,j,:,:) = envelope(real(squeeze(Xfilt(i,j,:,:))));
                   % X(i,j,:,:) = squeeze(S(i,j,:,:)).*envelope(squeeze(real(X(i,j,:,:))));      
                end
                waitbar(i/dims(1),b,'Basebanding');
            end
            if handles.bbdb.Value == 1
                Xfilt = 20*log10(Xfilt./max(max(max(max(Xfilt)))));
            end
            X = S.*abs(Xfilt);
            delete(b)
        end
    end
    if handles.invertbox.Value == 1
        X = X*(-1);
    end 
    assignin('base','X_c',X)
end



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
function overlay_Callback(hObject, eventdata, handles)
% hObject    handle to overlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



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
function noise_button_Callback(hObject, eventdata, handles)
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
function create_pe_Callback(hObject, eventdata, handles)
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
for i = 1:param.velmex.XNStep
    for j = 1:param.velmex.YNStep
         %HF(i,j,:,:) = envelope(real(X{i,j})); %Converts cell array to double
          HF(i,j,:,:) = X{i,j}; %Converts cell array to double
    end
    waitbar(i/param.velmex.XNStep,b,'Converting to 4D matrix');
end
delete(b)
PEdata = HF(:,:,:,1:2);
[~, ax] = make_axes(param,size(HF));
ax.pe = linspace(0,1.48*param.daq.HFdaq.pts/param.daq.HFdaq.fs_MHz/2,size(PEdata,3));
if handles.keep.Value == 1
    assignin('base','PEdata',PEdata);
    assignin('base','fpath',[path file]);
    assignin('base','param',param);
    assignin('base','ax',ax);
    set(handles.fname,'String',[path file]);
end
if handles.save_4d.Value == 1
    clearvars -except file path param ax PEdata
    f = file(1:end-4);
    f2 = [f '_4d_PE.mat'];
    fprintf('Saving 4D file...')
    save(f2);
    fprintf('Done\n')
end

else
    [f, p]  = uigetfile(fullfile(pwd,'*_PEParm.mat'));
load([p f]);
cd(p)
    [PEbsqFile, p] = uigetfile(fullfile(pwd,'*.bsq'));
    fid = fopen(PEbsqFile,'rb');
if fid > 0
    
    n = fread(fid,1,'int32');
    dsize = fread(fid,[1,n],'int32');
    nOffset = (n+1)*4;
    fclose(fid);
    
    PEImage = multibandread(PEbsqFile,[dsize(1:2),prod(dsize(3:end))],...
        'single',nOffset,'bsq','ieee-le',{'Band','Direct',bScanParm.nScanPt});
   
   
    
x = 3;
end
end


% --- Executes on button press in loadpe.
function loadpe_Callback(hObject, eventdata, handles)
[f,  p] = uigetfile(fullfile(pwd,'*4d_PE.mat'));
cd(p)
fprintf('Loading 4D Dataset...')
load([p f]);
%[~, ax] = make_axes(param,size(Xfilt),[1 2],1);
set(handles.fname,'String',file);
fprintf('Done\n')
assignin('base','PEdata',PEdata);
assignin('base','param',param);
assignin('base','ax',ax);
%assignin('base','PEparam',PE);
%set(handles.LF_chan,'String',num2str(size(LF,2)));
set(handles.tms,'String',num2str([ax.stime(1) ax.stime(end)]));
set(handles.tsamp,'String',num2str([1 length(ax.stime)]));
set(handles.xmm,'String',num2str([ax.x(1) ax.x(end)]));
set(handles.xsamp,'String',num2str([1 length(ax.x)]));
set(handles.ymm,'String',num2str([ax.y(1) ax.y(end)]));
set(handles.ysamp,'String',num2str([1 length(ax.y)]));
set(handles.zmm,'String',num2str([ax.pe(1) round(ax.pe(end))]));
set(handles.zsamp,'String',num2str([1 length(ax.depth)]));

if handles.reset_axes.Value == 1    
set(handles.xR,'String', num2str([ax.x(1) ax.x(end)]));
set(handles.yR,'String', num2str([ax.y(1) ax.y(end)]));
set(handles.zR,'String', num2str([ax.pe(1) floor(ax.pe(end))]));
set(handles.tR,'String', num2str([ax.stime(1) ax.stime(end)]));

end


% --- Executes on button press in usepe.
function usepe_Callback(hObject, eventdata, handles)
X = evalin('base','PEdata');
assignin('base','Xfilt',X);


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

% --- Executes on button press in savefig.
function savefig_Callback(hObject, eventdata, handles)
fname = handles.savefigname.String;
if ~isempty(handles.fignum.String)
    num = str2double(handles.fignum.String);
    if fname(end-2:end) == 'png'
        fname = [fname ' -transparent'];
        figure(num)
        set(gca,'Color','none');
    else
    end
    if ~isempty(handles.savefolder.String)
        pname = handles.savefolder.String;
        sname = [pname '\' fname];
    else
        sname = fname;
    end
    
    figure(num)
    eval([ 'export_fig ' sname])
else
    if ~isempty(handles.savefolder.String)
        pname = handles.savefolder.String;
        sname = [pname '\' fname];
    else
        sname = fname;
    end
    eval([ 'export_fig ' sname])
    
    
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
