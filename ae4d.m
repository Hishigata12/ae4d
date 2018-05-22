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

% Last Modified by GUIDE v2.5 22-May-2018 15:47:23

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
q.t = 1:dims(4);
xInd = q.x(find(ax.x >= xR(1)):find(ax.x >= xR(2)));
yInd = q.y(find(ax.y >= yR(1)):find(ax.y >= yR(2)));
zInd = q.z(find(ax.depth >= zR(1)):find(ax.depth >= zR(2)));
tInd = q.t(find(ax.stime >= tR(1)):find(ax.stime >= tR(2)));
Y = squeeze(Xfilt(xInd,yInd,zInd,tInd));
if length(size(Y)) > 2
    errordlg('Too many dimensions; check ranges')
    return
end
Y = medfilt2(Y,[3 3]);
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

fprintf('Loading 4D Dataset...')
load(uigetfile(pwd,'*4D_data.mat'));
[~, ax] = make_axes(param,size(Xfilt),[1 2],1);
set(handles.fname,'String',file);
fprintf('Done\n')
assignin('base','Xfilt',Xfilt);
assignin('base','fpath',[path file]);
assignin('base','param',param);
assignin('base','ax',ax);
assignin('base','LF',LF);
assignin('base','PEparam',PE);
set(handles.LF_chan,'String',num2str(size(LF,2)));

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
       % HF((i-1)*sL+j,:,:) = HF1{i,j}; %converts cells to pseudo 4-D array
         HF(i,j,:,:) = envelope(real(X{i,j})); %Converts cell array to double
        %X2(i,j,:,:) = envelope(squeeze(real(X(i,j,:,:))));
      % HF_bb(i,j,:,:) = X_bb{i,j};
    end
    waitbar(i/param.velmex.XNStep,b,'Converting to 4D matrix');
end
%%%%TESTING THIS OUT
% for i = 1:param.velmex.XNStep
%     for j = 1:param.velmex.YNStep
%         for k = 1:size(HF,4)
%             HF(i,j,:,k) = envelope(squeeze(HF(i,j,:,k)));
%         end
%     end
%     waitbar(.5+i/param.velmex.XNStep/3,b,'Fast Time Envelope');
% end


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
    save(f2);
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
yR = str2num(handles.yR.String);
if length(yR) == 1
    yR = [yR yR];
end
zR = str2num(handles.zR.String);
tR = str2num(handles.tR.String);
aeR = str2num(handles.aeR.String);
dims = size(Xfilt);
%[~,ax] = make_axes(param,dims,[1 2],1);
q.x = 1:dims(1);
q.y = 1:dims(2);
q.z = 1:dims(3);
q.t = 1:dims(4);
xInd = q.x(find(ax.x >= xR(1)):find(ax.x >= xR(2)));
yInd = q.y(find(ax.y >= yR(1)):find(ax.y >= yR(2)));
zInd = q.z(find(ax.depth >= zR(1)):find(ax.depth >= zR(2)));
tInd = q.t(find(ax.stime >= tR(1)):find(ax.stime >= tR(2)));

if handles.save_fig.Value == 0
    if handles.use_ext_fig.Value == 0
    axes(handles.axes2)
    else 
        figure(1);
    end
    for k = tInd
        J = medfilt2((squeeze(Xfilt(xInd,yInd,zInd,k)))',[3 3]);
      %  I = insertText(J,[0 200],['t = ' num2str(ax.stime(tInd))],'FontSize',14,...
           % 'BoxColor','green','TextColor','black'); 
       % I = rgb2gray(I);
        if handles.use_ext_fig.Value == 0
        imagesc(ax.x(xInd),ax.depth(zInd),J)
        colormap('hot')
        else 
            imshow(J)
            colormap(gca,'hot')
        end
        if ~isempty(aeR)
        caxis(aeR)
        end
        title(['t = ' num2str(ax.stime(k))]);
        handles.axes2.XLabel.String = 'Lateral (mm)';
        handles.axes2.YLabel.String = 'Depth (mm)';
        drawnow
    end
end



if handles.save_fig.Value == 1
    vidwrite(Xfilt,ax.depth,ax.x,[zInd(1) zInd(end)],[xInd(1) xInd(end)],[tInd(1) tInd(end)],aeR)
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
xR = str2num(handles.xR.String);
yR = str2num(handles.yR.String);
zR = str2num(handles.zR.String);
tR = str2num(handles.tR.String);
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
%[~,ax] = make_axes(param,size(X));
dims = size(X);

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
if isempty(num2str(handles.depR.String))
    qq = [10 60];
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
if handles.LF_FFT.Value == 1
    lf = fft(LF);
    x = linspace(0,param.daq.LFdaq.fs_Hz,length(lf));
    if handles.use_ext_fig.Value == 1
        figure(3)
        plot(x,abs(lf))
    else
        axes(handles.axes2)
        plot(x,abs(lf))         
    end
else
    x = linspace(0,param.daq.HFdaq.duration_ms,length(LF));
      if handles.use_ext_fig.Value == 1
        figure(3)
        plot(x,abs(LF))
    else
        axes(handles.axes2)
        plot(x,abs(LF))
      end
end
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


