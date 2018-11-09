function varargout = Jay(varargin)
% JAY MATLAB code for Jay.fig
%      JAY, by itself, creates a new JAY or raises the existing
%      singleton*.
%
%      H = JAY returns the handle to a new JAY or the handle to
%      the existing singleton*.
%
%      JAY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in JAY.M with the given input arguments.
%
%      JAY('Property','Value',...) creates a new JAY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Jay_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Jay_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Jay

% Last Modified by GUIDE v2.5 08-Nov-2018 13:32:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Jay_OpeningFcn, ...
                   'gui_OutputFcn',  @Jay_OutputFcn, ...
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


% --- Executes just before Jay is made visible.
function Jay_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Jay (see VARARGIN)

% Choose default command line output for Jay
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Jay wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Jay_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
X = evalin('base','Xmerged');
assignin('base','X',X);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, path] = uigetfile('*.mat');
P = load([path file],'Pressure');
Pressure = P.Pressure;
assignin('base','Pressure',Pressure);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, path] = uigetfile('*.mat');
L2 = load([path file],'L');
L = L2.L;
assignin('base','L',L);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, path] = uigetfile('*.mat');
J2 = load([path file],'J');
J = J2.J;
assignin('base','J',J);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
J = evalin('base','J');
if handles.three.Value
    if handles.usepressure.Value
        Pressure = evalin('base','Pressure');
        %Method 1
        dimsP = round(size(Pressure)/2);
        q = round(str2double(get(handles.psize,'String'))/2);
        P = Pressure(dimsP(1)-q:dimsP(1)+q,dimsP(2)-q:dimsP(2)+q,dimsP(2)-q:dimsP(2)+q);
        Pressure = P;
        %%%%%%%%%%%%
        dimsP = size(Pressure);
    end
     dimsJ = size(J);
    if handles.useleadfield.Value
        L = evalin('base','L');
          V = zeros(dimsJ(1),dimsJ(2),dimsJ(3),size(L,5));
    else
        V = zeros(dimsJ(1),dimsJ(2),dimsJ(3));
    end
    F = get(handles.fullvol,'Value');
    if F
        J0 = padarray(J,[q q q],'symmetric','both');
        L0 = padarray(L,[q q q],'symmetric','both');
        clear L J
        J = J0;
        L = L0;
        clear L0 J0
        dimsJ = size(J);
    end
  
    if handles.useleadfield.Value
        Y = zeros(dimsJ(1),dimsJ(2),dimsJ(3),size(L,5));
        for m = 1:size(L,5)
            Y(:,:,:,m) = dot(J,L(:,:,:,:,m),4);
        end
        for m = 1:size(L,5)
            for i = floor(dimsP(1)/2)+1:dimsJ(1)-floor(dimsP(1)/2)
                for j = floor(dimsP(2)/2)+1:dimsJ(2)-floor(dimsP(2)/2)
                    for k = floor(dimsP(3)/2)+1:dimsJ(3)-floor(dimsP(3)/2)
                        %                         Y = dot(J(i-floor(dimsP(1)/2):i+floor(dimsP(1)/2),...
                        %                             j-floor(dimsP(2)/2):j+floor(dimsP(2)/2),k-floor(dimsP(3)/2):k+floor(dimsP(3)/2),:),...
                        %                             L(i-floor(dimsP(1)/2):i+floor(dimsP(1)/2),j-floor(dimsP(2)/2):j+floor(dimsP(2)/2),...
                        %                             k-floor(dimsP(3)/2):k+floor(dimsP(3)/2),:,m),4);
                        Y2 = Y(i-floor(dimsP(1)/2):i+floor(dimsP(1)/2),...
                            j-floor(dimsP(2)/2):j+floor(dimsP(2)/2),k-floor(dimsP(3)/2):k+floor(dimsP(3)/2),m);
                        Y2 = reshape(Y2,[numel(Y2),1]);
                        P2 = reshape(Pressure,[numel(Pressure),1]);
                        V(i,j,k,m) = dot(Y2,P2);
                        % V(i,j,k,m) = dot(dot(J(i-floor(dimsP(1)/2):i+floor(dimsP(1)/2),j-floor(dimsP(2)/2):j+floor(dimsP(2)/2),k-floor(dimsP(3)/2):j+floor(dimsP(3)/2),:),...
                        % L(i-floor(dimsP(1)/2):i+floor(dimsP(1)/2),j-floor(dimsP(2)/2):j+floor(dimsP(2)/2),k-floor(dimsP(3)/2):j+floor(dimsP(3)/2),:,m),4),Pressure);
                    end
                end
                multiWaitbar(['Solving for Vae: Channel ' num2str(m) ' of ' num2str(size(L,5))],i/(dimsJ(1)-dimsP(1)));
            end   
        end
          V(dimsJ(1)-dimsP(1)+1:dimsJ(1),dimsJ(2)-dimsP(2)+1:dimsJ(2),dimsJ(3)-dimsP(3)+1:dimsJ(3),1:size(L,5)) = 0;
    else
        Y = zeros(dimsJ(1),dimsJ(2),dimsJ(3),size(J,4));
        for m = 1:size(J,4)
            for i = floor(dimsP(1)/2)+1:dimsJ(1)-floor(dimsP(1)/2)
                for j = floor(dimsP(2)/2)+1:dimsJ(2)-floor(dimsP(2)/2)
                    for k = floor(dimsP(3)/2)+1:dimsJ(3)-floor(dimsP(3)/2)
                        
                        Y = J(i-floor(dimsP(1)/2):i+floor(dimsP(1)/2),...
                        j-floor(dimsP(2)/2):j+floor(dimsP(2)/2),k-floor(dimsP(3)/2):k+floor(dimsP(3)/2),m);
                        Y2 = reshape(Y,[numel(Y),1]);
                        P2 = reshape(Pressure,[numel(Pressure),1]);
                        V(i,j,k,m) = dot(Y2,P2);
                        
                        % V(i,j,k,m) = dot(dot(J(i-floor(dimsP(1)/2):i+floor(dimsP(1)/2),j-floor(dimsP(2)/2):j+floor(dimsP(2)/2),k-floor(dimsP(3)/2):j+floor(dimsP(3)/2),:),...
                        % L(i-floor(dimsP(1)/2):i+floor(dimsP(1)/2),j-floor(dimsP(2)/2):j+floor(dimsP(2)/2),k-floor(dimsP(3)/2):j+floor(dimsP(3)/2),:,m),4),Pressure);
                    end
                end
                multiWaitbar(['Solving for Vae: Component ' num2str(m) ' of ' num2str(size(J,4))],i/(dimsJ(1)-dimsP(1)));
            end 
        end
          V(dimsJ(1)-dimsP(1)+1:dimsJ(1),dimsJ(2)-dimsP(2)+1:dimsJ(2),dimsJ(3)-dimsP(3)+1:dimsJ(3),1:size(J,4)) = 0;
    end
  
    assignin('base','V',V)
else
end
multiWaitbar('CLOSEALL');
    

% --- Executes on button press in usepressure.
function usepressure_Callback(hObject, eventdata, handles)
% hObject    handle to usepressure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of usepressure


% --- Executes on button press in useleadfield.
function useleadfield_Callback(hObject, eventdata, handles)
% hObject    handle to useleadfield (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of useleadfield


% --- Executes on button press in three.
function three_Callback(hObject, eventdata, handles)
% hObject    handle to three (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of three


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.paramenu.Value == 1
    data = evalin('base',V);
    tit1 = 'V calc ';
elseif handles.paramenu.Value == 2
      data = evalin('base',X);
      tit1 = 'V exp ';
elseif handles.paramenu.Value == 3
      data = evalin('base',J);
      tit1 = 'J ';
elseif handles.paramenu.Value == 4
      data = evalin('base',L);
      tit1 = 'L ';
elseif handles.paramenu.Value == 5
      data = evalin('base',Pressure);
      tit1 = 'Pressure ';
end
comp = get(handles.compmenu,'Value');
chan = get(handles.chanmenu,'Value');
px = str2double(handles.px.String);
py = str2double(handles.py.String);
pz = str2double(handles.pz.String);
pt = str2double(handles.pt.String);
xl = str2num(handles.xlims.String);
yl = str2num(handles.ylims.String);
if handles.plotmenu.Value == 1
    if handles.paramenu.Value == 4
    data = squeeze(data(:,py,:,comp,chan));
    elseif handles.paramenu.Value == 3
       data = squeeze(data(:,py,:,comp));
    elseif handles.paramenu.Value == 2
       data = squeeze(data(:,py,:,pt,chan));
    elseif handles.paramenu.Value == 1
       data = squeeze(data(:,py,:,chan));
    elseif handles.paramenu.Value == 5
       data = squeeze(data(:,py,:));
    end
    data = data';
    tit2 = 'X-Z plot';
    xlab = 'Lateral';
    ylab = 'Depth';
end
axes(handles.axes1);
imagesc(data);
set(handles.axes1,'XLabel',xlab);
set(handles.axes1,'YLabel',ylab);
set(handles.axes1,'Title',[tit1 tit2]);
set(handles.axes1,'XLim',xl);
set(handles.axes1,'YLim',yl);

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.exportv.Value
    V = evalin('base','V');
    assignin('base','X_c',V);
elseif handles.exportj.Value
    J = evalin('base','J');
    assignin('base','X_c',J);
else 
    errordlg('Select either J or V to export to X_c');
end


% --- Executes on button press in exportv.
function exportv_Callback(hObject, eventdata, handles)
% hObject    handle to exportv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of exportv


% --- Executes on button press in exportj.
function exportj_Callback(hObject, eventdata, handles)
% hObject    handle to exportj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of exportj


% --- Executes on selection change in parammenu.
function parammenu_Callback(hObject, eventdata, handles)
% hObject    handle to parammenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.paramenu.Value == 1
    handles.chanmenu.String = 1:size(V,4);
elseif handles.paramenu.Value == 2
    handles.chanmenu.String = 1:size(X,5);
elseif handles.paramenu.Value == 3
    handles.chanmenu.String = [];
elseif handles.paramenu.Value == 4
    handles.chanmenu.String = 1:size(L,5);
elseif handles.paramenu.Value == 5
    handles.chanmenu.String = [];
end
    

% Hints: contents = cellstr(get(hObject,'String')) returns parammenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from parammenu


% --- Executes during object creation, after setting all properties.
function parammenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parammenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in chanmenu.
function chanmenu_Callback(hObject, eventdata, handles)
% hObject    handle to chanmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns chanmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from chanmenu


% --- Executes during object creation, after setting all properties.
function chanmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chanmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in compmenu.
function compmenu_Callback(hObject, eventdata, handles)
% hObject    handle to compmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns compmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from compmenu


% --- Executes during object creation, after setting all properties.
function compmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to compmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plotmenu.
function plotmenu_Callback(hObject, eventdata, handles)
% hObject    handle to plotmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plotmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotmenu


% --- Executes during object creation, after setting all properties.
function plotmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px_Callback(hObject, eventdata, handles)
% hObject    handle to px (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px as text
%        str2double(get(hObject,'String')) returns contents of px as a double


% --- Executes during object creation, after setting all properties.
function px_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function py_Callback(hObject, eventdata, handles)
% hObject    handle to py (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of py as text
%        str2double(get(hObject,'String')) returns contents of py as a double


% --- Executes during object creation, after setting all properties.
function py_CreateFcn(hObject, eventdata, handles)
% hObject    handle to py (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pz_Callback(hObject, eventdata, handles)
% hObject    handle to pz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pz as text
%        str2double(get(hObject,'String')) returns contents of pz as a double


% --- Executes during object creation, after setting all properties.
function pz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function pt_Callback(hObject, eventdata, handles)
% hObject    handle to pt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pt as text
%        str2double(get(hObject,'String')) returns contents of pt as a double


% --- Executes during object creation, after setting all properties.
function pt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function xli_Callback(hObject, eventdata, handles)
% hObject    handle to xli (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xli as text
%        str2double(get(hObject,'String')) returns contents of xli as a double


% --- Executes during object creation, after setting all properties.
function xli_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xli (see GCBO)
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
