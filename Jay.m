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

% Last Modified by GUIDE v2.5 15-Nov-2018 12:00:13

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
function pushbutton5_Callback(hObject, eventdata, handles) %Inverse
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
V = evalin('base','V');
L = evalin('base','L');
%P = evalin('base','Pressure');
if handles.usepressure.Value
    Pressure = evalin('base','Pressure');
    %This limits the size of the pressure matrix to P Size
    %input
    
    %Step one: deconvolve P from V
    %Method 3 FFT division
    if handles.deconvmethod.Value == 3
        Vs = size(V);
        Ps = size(Pressure);
        [x,y,z] = meshgrid(linspace(1,Ps(1),Ps(1)),linspace(1,Ps(2),Ps(2)),linspace(1,Ps(3),Ps(3)));
        [x1,y1,z1] = meshgrid(linspace(1,Ps(1),Vs(1)),linspace(1,Ps(2),Vs(2)),linspace(1,Ps(3),Vs(3)));
        P2 = interp3(x,y,z,Pressure,x1,y1,z1);
        HPF(1:101) = linspace(1,0,101);
        HPF(102:201) = linspace(0,1,100);
        bpFilt = designfilt('bandpassfir','FilterOrder',20, ...
            'CutoffFrequency1',2000,'CutoffFrequency2',4000, ...
            'SampleRate',20000);
        for m = 1%1:size(V,4)
            for i = 1:size(V,1)
                for j = 1:size(V,2)
                    Vf = fft2(squeeze(V(i,j,:,m)));
                    Pf = fft2(squeeze(P2(101,101,:)))+1e5;
                    %Pf2 =
                    %interp1(linspace(1,length(Pf),length(Pf)),Pf,linspace(1,length(Pf),length(Vf)));
                    T1 = ifft(Vf./Pf);
                    % T2 = T1.*HPF';
                    J1(i,j,:,m) = filtfilt(bpFilt,T1);
                    %                     if mod(j,10)==0
                    %                     plot(abs(T2));
                    %                     drawnow;
                    %                     end
                    % J1(i,j,:,m) = real(ifft(T2));
                end
                multiWaitbar(['Solving for Vae: Channel ' num2str(m) ' of ' num2str(size(V,4))],i/(size(V,1)));
            end
        end
        cut = round(size(X,3)*.12);
        J1(:,:,1:cut,:) = 0;
        J1(:,:,size(X,3)-cut:size(X,3),:) = 0;
        
        %Method 2 Limits size of Pressure then manually deconvolves in space
    elseif handles.deconvmethod.Value == 2
        dimsP = round(size(Pressure)/2);
        q = round(str2double(get(handles.psize,'String'))/2);
        P = Pressure(dimsP(1)-q:dimsP(1)+q,dimsP(2)-q:dimsP(2)+q,dimsP(2)-q:dimsP(2)+q);
        Pressure = P;
        %%%%%%%%%%%%
        dimsP = size(Pressure);
        dimsJ = size(V);
        J1 = zeros(size(V,1),size(V,2),size(V,3),size(V,4));
        for m = 1%1:size(V,4)
            for i = floor(dimsP(1)/2)+1:dimsJ(1)-floor(dimsP(1)/2)
                for j = floor(dimsP(2)/2)+1:dimsJ(2)-floor(dimsP(2)/2)
                    for k = floor(dimsP(3)/2)+1:dimsJ(3)-floor(dimsP(3)/2)
                        
                        Y = V(i-floor(dimsP(1)/2):i+floor(dimsP(1)/2),...
                            j-floor(dimsP(2)/2):j+floor(dimsP(2)/2),k-floor(dimsP(3)/2):k+floor(dimsP(3)/2),m);
                        Y2 = reshape(Y,[numel(Y),1]);
                        P2 = reshape(Pressure,[numel(Pressure),1]);
                        %                         noisef = min(abs(P2));
                        %                         P2 = P2+noisef/2;
                        namp = max(P2)/1e4;
                        noyse = P2 + namp*rand(size(P2))+1;
                        %                         P4 = awgn(P2,30,'measured');
                        P3 = 1./noyse;
                        %                         J1(i,j,k,m) = dot(Y2,P3);
                        J1(i,j,k,m) = sum(Y2./noyse);
                        
                        % V(i,j,k,m) = dot(dot(J(i-floor(dimsP(1)/2):i+floor(dimsP(1)/2),j-floor(dimsP(2)/2):j+floor(dimsP(2)/2),k-floor(dimsP(3)/2):j+floor(dimsP(3)/2),:),...
                        % L(i-floor(dimsP(1)/2):i+floor(dimsP(1)/2),j-floor(dimsP(2)/2):j+floor(dimsP(2)/2),k-floor(dimsP(3)/2):j+floor(dimsP(3)/2),:,m),4),Pressure);
                    end
                end
                multiWaitbar(['Solving for Vae: Channel ' num2str(m) ' of ' num2str(size(V,4))],i/(dimsJ(1)-dimsP(1)));
            end
        end
        J1(dimsJ(1)-dimsP(1)+1:dimsJ(1),dimsJ(2)-dimsP(2)+1:dimsJ(2),dimsJ(3)-dimsP(3)+1:dimsJ(3),1:size(V,4)) = 0;
        
        
        %Method 1 built in deconv
    elseif handles.deconvmethod.Value == 1
        if length(size(V)) == 3
            J1 = deconvlucy(V,Pressure);
        elseif length(size(V)) == 4
            for i = 1:size(V,4)
                J1(:,:,:,i) = deconvlucy(V(:,:,:,i),Pressure);
                multiWaitbar('Deconvolving Pressure',i/size(V,4));
            end
        elseif length(size(V)) == 5
            for i = 1:size(V,4)
                for j = 1:size(V,5)
                    J1(:,:,:,i,j) = deconvlucy(V(:,:,:,i,j),Pressure);
                end
                multiWaitbar('Deconvolving Pressure',i/size(V,4));
            end
        end
        %Method 4 basebanding to pressure freq
    elseif handles.deconvmethod.Value == 4
        Ps = size(Pressure);
        p1 = round(Ps(1)/2);
        p2 = round(Ps(2)/2);
        Pf = fft(squeeze(Pressure(p1,p2,:)));
        [val,loc] = max(Pf);
        wc = loc/str2double(handles.ws.String);
        if isempty(handles.wc.String)
            J2 = baseband2(V(:,:,:,1),wc,str2double(handles.ws.String),2,4); %change bandpass for other Trans
        else
            J2 = baseband2(V(:,:,:,1),str2double(handles.wc.String),str2double(handles.ws.String),2,4);
            %Need to scale this by pressure
        end
        P = abs(hilbert(Pressure));
        for i = 1:size(J2,4)
            J3(:,:,:,i) = deconvlucy(real(J2(:,:,:,i)),P);
        end
        J1 = real(J3);
    end
else
    J1 = V;
end


if handles.useleadfield.Value
    Vs = size(J1);
    Ls = size(L);
    dif = Vs(1:3)-Ls(1:3);
    if dif ~= 0
        J2 = J1;
        clear J1
        [x,y,z] = meshgrid(linspace(1,Vs(1),Vs(1)),linspace(1,Vs(2),Vs(2)),linspace(1,Vs(3),Vs(3)));
        [x1,y1,z1] = meshgrid(linspace(1,Vs(1),Ls(1)),linspace(1,Vs(2),Ls(2)),linspace(1,Vs(3),Ls(3)));
        J1 = interp3(x,y,z,J2,x1,y1,z1);
    end
    
    for i = 1:size(J1,1)
        for j = 1:size(J1,2)
            for k = 1:size(J1,3)
                for m = 1:size(L,4)
                %J(i,j,k,:) = squeeze(L(i,j,k,:,:))\squeeze(V(i,j,k,:));
                %                     J2(i,j,k,m) = L/J1
                J(i,j,k,m) = dot(pinv(squeeze(L(i,j,k,m,:))),squeeze(J1(i,j,k,:)));
                end
            end
        end
        multiWaitbar('Calculating J(x,y,z)',i/size(J1,1));
    end
else
    J = J1;
end

multiWaitbar('CLOSEALL');
assignin('base','Jrecon',J);


% --- Executes on button press in pushbutton6. 
function pushbutton6_Callback(hObject, eventdata, handles) %Forward
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
        Y = zeros(dimsJ(1),dimsJ(2),dimsJ(3));
       % for m = 1:size(J,4)
       for i = 1:size(J,1)
           for j = 1:size(J,2)
               for k = 1:size(J,3)
                   Jm(i,j,k) = sqrt(J(i,j,k,1)^2+J(i,j,k,2)^2+J(i,j,k,3)^2);
               end
           end
           multiWaitbar('Converting to magnitude of J',i/size(J,1))
       end
            for i = floor(dimsP(1)/2)+1:dimsJ(1)-floor(dimsP(1)/2)
                for j = floor(dimsP(2)/2)+1:dimsJ(2)-floor(dimsP(2)/2)
                    for k = floor(dimsP(3)/2)+1:dimsJ(3)-floor(dimsP(3)/2)
                        Y = Jm(i-floor(dimsP(1)/2):i+floor(dimsP(1)/2),...
                        j-floor(dimsP(2)/2):j+floor(dimsP(2)/2),k-floor(dimsP(3)/2):k+floor(dimsP(3)/2));
                        Y2 = reshape(Y,[numel(Y),1]);
                        P2 = reshape(Pressure,[numel(Pressure),1]);
                        V(i,j,k) = dot(Y2,P2);
                        
                        % V(i,j,k,m) = dot(dot(J(i-floor(dimsP(1)/2):i+floor(dimsP(1)/2),j-floor(dimsP(2)/2):j+floor(dimsP(2)/2),k-floor(dimsP(3)/2):j+floor(dimsP(3)/2),:),...
                        % L(i-floor(dimsP(1)/2):i+floor(dimsP(1)/2),j-floor(dimsP(2)/2):j+floor(dimsP(2)/2),k-floor(dimsP(3)/2):j+floor(dimsP(3)/2),:,m),4),Pressure);
                    end
                end
                multiWaitbar('Solving for Vae',i/(dimsJ(1)-dimsP(1)));
            end 
        %end
          V(dimsJ(1)-dimsP(1)+1:dimsJ(1),dimsJ(2)-dimsP(2)+1:dimsJ(2),dimsJ(3)-dimsP(3)+1:dimsJ(3)) = 0;
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



function psize_Callback(hObject, eventdata, handles)
% hObject    handle to psize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of psize as text
%        str2double(get(hObject,'String')) returns contents of psize as a double


% --- Executes during object creation, after setting all properties.
function psize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to psize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fullvol.
function fullvol_Callback(hObject, eventdata, handles)
% hObject    handle to fullvol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fullvol


% --- Executes on selection change in deconvmethod.
function deconvmethod_Callback(hObject, eventdata, handles)
% hObject    handle to deconvmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns deconvmethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from deconvmethod


% --- Executes during object creation, after setting all properties.
function deconvmethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deconvmethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function wc_Callback(hObject, eventdata, handles)
% hObject    handle to wc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of wc as text
%        str2double(get(hObject,'String')) returns contents of wc as a double


% --- Executes during object creation, after setting all properties.
function wc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ws_Callback(hObject, eventdata, handles)
% hObject    handle to ws (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ws as text
%        str2double(get(hObject,'String')) returns contents of ws as a double


% --- Executes during object creation, after setting all properties.
function ws_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ws (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
