function varargout = er_UI_GLMplot(varargin)
% ER_UI_GLMPLOT MATLAB code for er_UI_GLMplot.fig
%      ER_UI_GLMPLOT, by itself, creates a new ER_UI_GLMPLOT or raises the existing
%      singleton*.
%
%      H = ER_UI_GLMPLOT returns the handle to a new ER_UI_GLMPLOT or the handle to
%      the existing singleton*.
%
%      ER_UI_GLMPLOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ER_UI_GLMPLOT.M with the given input arguments.
%
%      ER_UI_GLMPLOT('Property','Value',...) creates a new ER_UI_GLMPLOT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before er_UI_GLMplot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to er_UI_GLMplot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help er_UI_GLMplot

% Last Modified by GUIDE v2.5 18-Mar-2019 22:47:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @er_UI_GLMplot_OpeningFcn, ...
                   'gui_OutputFcn',  @er_UI_GLMplot_OutputFcn, ...
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


% --- Executes just before er_UI_GLMplot is made visible.
function er_UI_GLMplot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to er_UI_GLMplot (see VARARGIN)
global sel_cont SelY

sel_cont = 1;
handles.STATS = varargin{1};
handles.Y = varargin{2};

% defines contrast list
c_number = size(handles.STATS.model.C,1);
str = cell(c_number,1);
for l =1:c_number
    str{l} = num2str(handles.STATS.model.C(l,:));
end
set(handles.popupmenu1,'String',str);
set(handles.popupmenu1,'Value',1);

%axes(handles.axes1);
axis_ctrl(handles,'axes1',handles.STATS.CB,sel_cont,[],'Contrast*B')
axis_ctrl(handles,'axes2',handles.STATS.T,sel_cont, [],'T-value')
axis_ctrl(handles,'axes3',handles.STATS.P,sel_cont,[0 0.05],'P-value')

set(handles.radiobutton1,'Value',1);set(handles.radiobutton2,'Value',0);set(handles.radiobutton3,'Value',0)
handles.Pmode = 1;


%set slider and edit box for single variable plot
maxSlX = size(handles.Y,2)-1;
maxSlY = size(handles.Y,3)-1;
if maxSlX > 0
    set(handles.SlX,'min',0,'max',maxSlX,'sliderstep',[1/maxSlX 10/maxSlX]);
else
    set(handles.SlX,'min',0,'max',maxSlX,'sliderstep',[0 1]);
end
if maxSlY > 0
    set(handles.SlY,'min',0,'max',maxSlY,'sliderstep',[1/maxSlY 10/maxSlY]);
else
    set(handles.SlY,'min',0,'max',maxSlY,'sliderstep',[0 1]);
end



if not(isfield(SelY,'row')) || (isfield(SelY,'row') && SelY.row > (maxSlX+1) )
    SelY.row = 1;
    set(handles.edCoX,'String','1');
else
    set(handles.edCoX,'String',num2str(SelY.row));
end
if not(isfield(SelY,'col')) || (isfield(SelY,'col') && SelY.col > (maxSlY+1) )
    SelY.col = 1;
    set(handles.edCoY,'String','1');
else
    set(handles.edCoY,'String',num2str(SelY.col));
end
plot_contrast_fit(handles,SelY.row,SelY.col);


    % 
% set(handles.SlY,'min',0);
% set(handles.SlX,'max',size(handles.Y,2)-1);
% set(handles.SlY,'max',size(handles.Y,3)-1);


% Choose default command line output for er_UI_GLMplot
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes er_UI_GLMplot wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function axis_ctrl(handles,axes_,Y,Cindx,limits,title_)


switch axes_
    case 'axes1'
        ceckbox =get(handles.checkbox1,'Value');
        edn1 = '2';
        edn2 = '3';           
    case 'axes2'
        ceckbox =get(handles.checkbox2,'Value');
        edn1 = '4';
        edn2 = '5'; 
    case 'axes3'
        checkbox = 0;
end

eval(['axes(handles.',axes_,')'])
Y = squeeze(Y(Cindx,:,:));
if isempty(limits)
    if ceckbox
        maxx = max(Y(Y ~= inf));
        minn = min(Y(Y ~= inf));
        % for simmetric colorbar
        if abs(maxx) > abs(minn)
            minn = -abs(maxx);
            maxx = abs(maxx);
        else
            minn = -abs(minn);
            maxx = abs(minn);
        end
        eval(['set(handles.edit',edn1,',''String'',num2str(minn));']);
        eval(['set(handles.edit',edn2,',''String'',num2str(minn));']);
    else
        eval(['minn = str2double(get(handles.edit',edn1,',''String''));']);
        eval(['maxx = str2double(get(handles.edit',edn2,',''String''));']);
    end
else
    minn = limits(1);
    maxx = limits(2);
end
imagesc(Y,[minn maxx]);colorbar; colormap jet
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);
title(title_);




% --- Outputs from this function are returned to the command line.
function varargout = er_UI_GLMplot_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sel_cont
sel_cont = get(hObject,'Value');
axis_ctrl(handles,'axes1',handles.STATS.CB,sel_cont,[],'Contrast*B')
axis_ctrl(handles,'axes2',handles.STATS.T,sel_cont,[],'T-value')
axis_ctrl(handles,'axes3',handles.STATS.P,sel_cont,[0 0.05],'P-value')

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value') && handles.Pmode~=1
    set(handles.radiobutton2,'Value',0); set(handles.radiobutton3,'Value',0);
    axis_ctrl(handles,'axes3',handles.STATS.P,1,[0 0.05],'P-value');
    %set(hObject,'Value',1)
    handles.Pmode = 1;
end
guidata(hObject, handles);

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value') && handles.Pmode~=2
    set(handles.radiobutton1,'Value',0); set(handles.radiobutton3,'Value',0);
    axis_ctrl(handles,'axes3',handles.STATS.Pfdr,1,[0 0.05],'P-value');
    %set(hObject,'Value',1);
    handles.Pmode = 2;
end
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value') && handles.Pmode~=3
    set(handles.radiobutton1,'Value',0); set(handles.radiobutton2,'Value',0);
    axis_ctrl(handles,'axes3',handles.STATS.P*numel(handles.STATS.P(1,:,:)),1,[0 0.05],'P-value');
    %set(hObject,'Value',1);
    handles.Pmode = 3;
end
guidata(hObject, handles);
% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sel_cont
if get(hObject,'Value')
    set(handles.edit2,'Enable','off');
    set(handles.edit3,'Enable','off');
    axis_ctrl(handles,'axes1',handles.STATS.CB,sel_cont,[],'Contrast*B')
else
    set(handles.edit2,'Enable','on');
    set(handles.edit3,'Enable','on');
end

% Hint: get(hObject,'Value') returns toggle state of checkbox1



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes1);
lim = caxis;
caxis([str2double(get(hObject,'String')),lim(2)]);
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
lim = caxis;
caxis([lim(1), str2double(get(hObject,'String'))]);
% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sel_cont
if get(hObject,'Value')
    set(handles.edit4,'Enable','off');
    set(handles.edit5,'Enable','off');
    axis_ctrl(handles,'axes2',handles.STATS.T,sel_cont,[],'T-value')
else
    set(handles.edit4,'Enable','on');
    set(handles.edit5,'Enable','on');
end

% Hint: get(hObject,'Value') returns toggle state of checkbox2



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes2);
lim = caxis;
caxis([str2double(get(hObject,'String')),lim(2)]);
% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes2);
lim = caxis;
caxis([lim(1), str2double(get(hObject,'String'))]);
% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edCoX_Callback(hObject, eventdata, handles)
% hObject    handle to edCoX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edCoX as text
%        str2double(get(hObject,'String')) returns contents of edCoX as a double


% --- Executes during object creation, after setting all properties.
function edCoX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edCoX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edCoY_Callback(hObject, eventdata, handles)
% hObject    handle to edCoY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edCoY as text
%        str2double(get(hObject,'String')) returns contents of edCoY as a double


% --- Executes during object creation, after setting all properties.
function edCoY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edCoY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% % --- Executes on button press in pb_plot.
% function pb_plot_Callback(hObject, eventdata, handles)
% % hObject    handle to pb_plot (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes4


% 
% function edit8_Callback(hObject, eventdata, handles)
% % hObject    handle to edCoX (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of edCoX as text
% %        str2double(get(hObject,'String')) returns contents of edCoX as a double
% 
% 
% % --- Executes during object creation, after setting all properties.
% function edit8_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edCoX (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end
% 
% 
% 
% function edit9_Callback(hObject, eventdata, handles)
% % hObject    handle to edCoY (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'String') returns contents of edCoY as text
% %        str2double(get(hObject,'String')) returns contents of edCoY as a double
% 
% 
% % --- Executes during object creation, after setting all properties.
% function edit9_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to edCoY (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: edit controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


% --- Executes on button press in PuBoPlot.
function PuBoPlot_Callback(hObject, eventdata, handles)
% hObject    handle to PuBoPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global SelY
SelY.row = str2double(get(handles.edCoX,'String'));
SelY.col = str2double(get(handles.edCoY,'String'));
plot_contrast_fit(handles,SelY.row,SelY.col);

function plot_contrast_fit(handles,x,y)
global sel_cont

X = handles.STATS.model.X;
C = handles.STATS.model.C(sel_cont,:);
selY = handles.Y(:,x,y);

if strcmp(handles.STATS.model.ynan,'remove')
    %check for NaNs
    NaNy = isnan(selY);
    if sum(NaNy) > 0
        selY(NaNy) = [];
        X(NaNy) = [];
    end
    if handles.STATS.model.zscore 
        selY = zscore(selY);
        X = zscore(X);
    end
end
        
predictor = X*C';

axes(handles.axes6);

scatter(predictor,selY,'MarkerEdgeColor',[0.3 0.3 0.3],'MarkerFaceColor',[0.3 0.3 0.3]);
ylabel('response variable (Y)');
xlabel('predictor');
if handles.STATS.model.zscore 
    xlim([-4 4]);
    ylim([-4 4]);
    box on;
    hold on
    xL = xlim; yL = ylim;
    line([0 0], yL,'LineStyle','--','color','k');
    line(xL,[0 0], 'LineStyle','--','color','k');
    hold off;
end



% --- Executes on slider movement.
function SlX_Callback(hObject, eventdata, handles)
% hObject    handle to SlX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global SelY
SelY.row = (int32(get(hObject,'Value'))+1);
set(handles.edCoX,'String',num2str(SelY.row));
plot_contrast_fit(handles,SelY.row,SelY.col)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function SlX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SlX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function SlY_Callback(hObject, eventdata, handles)
% hObject    handle to SlY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global SelY
SelY.col = (int32(get(hObject,'Value'))+1);
set(handles.edCoY,'String',num2str(SelY.col));
plot_contrast_fit(handles,SelY.row,SelY.col)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function SlY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SlY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes3);
lim = caxis;
caxis([lim(1), str2double(get(hObject,'String'))]);

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
