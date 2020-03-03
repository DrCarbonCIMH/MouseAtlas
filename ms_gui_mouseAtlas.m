function varargout = ms_gui_mouseAtlas(varargin)
% MS_GUI_MOUSEATLAS MATLAB code for ms_gui_mouseAtlas.fig
%      MS_GUI_MOUSEATLAS, by itself, creates a new MS_GUI_MOUSEATLAS or raises the existing
%      singleton*.
%
%      H = MS_GUI_MOUSEATLAS returns the handle to a new MS_GUI_MOUSEATLAS or the handle to
%      the existing singleton*.
%
%      MS_GUI_MOUSEATLAS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MS_GUI_MOUSEATLAS.M with the given input arguments.
%
%      MS_GUI_MOUSEATLAS('Property','Value',...) creates a new MS_GUI_MOUSEATLAS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ms_gui_mouseAtlas_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ms_gui_mouseAtlas_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ms_gui_mouseAtlas

% Last Modified by GUIDE v2.5 02-Mar-2018 06:07:55

% Version 0.1; main purpose of this application is to create nifti files
% of selected ROIs for later usage as 'binary' masks. The visualization of
% the cooresponding ROIs is just for illustration purposes. Please note,
% that the use of isosurface does not result always in completely correct
% illustrations!
% Author: Markus Sack; markus.sack@zi-mannheim.de
% The software comes with no warranty. Always double check your data!

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ms_gui_mouseAtlas_OpeningFcn, ...
                   'gui_OutputFcn',  @ms_gui_mouseAtlas_OutputFcn, ...
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


% --- Executes just before ms_gui_mouseAtlas is made visible.
function ms_gui_mouseAtlas_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ms_gui_mouseAtlas (see VARARGIN)

% Choose default command line output for ms_gui_mouseAtlas
handles.output = hObject;

% load the annotation matrix
% handles.atlas='/adata1/markus/Allen_mouse_atlas/version3/CCFV3/annotation_50.nii';
d=fileparts(which('ms_gui_mouseAtlas'));
load([d filesep 'ms_AllenStuff.mat']);
if ~isdeployed
    if isempty(which('spm')) || ~isempty(strfind(which('spm'),'spm8'))
        errordlg('Please, first add your SPM 12 path and start the application again','No SPM found');
    end
    if isempty(which('ms_allen2pax_batch'))
%         addpath([d filesep '20171030_Allen2Pax']);
        addpath([d filesep 'transJobs']);
    end
end

handles.Vavg=Vavg; % that is for the avg template
handles.Mtxavg=Mtxavg; % that is for the avg template

switchMatrix=[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
handles.RotMat=V.mat * switchMatrix; % for later transformation of ROIs
handles.Vatlas = V;
handles.Mtx=Mtx; % this is the annotation info
% create a simple brain; just the outlines of the annotation matrix
tmp = handles.Mtx; tmp(tmp>0)=1; tmp(isnan(tmp))=0;
handles.brain=tmp; clear tmp

handles.shape = '<html> <font color="green"><b>'; % that defines the look of selected ROIs
% import java function for the nodes
import javax.swing.*
import javax.swing.tree.*
% import com.mathworks.mwswing.checkboxtree.* % that is somewhat too complicated!

% create the node tree
root=uitreenode('v0', (regions{1}.id),[regions{1}.name],[],false);
[tree, container]=uitree('v0','Root',root,'ExpandFcn',{@build_tree,regions},'SelectionChangeFcn',{@mySelectFcn, handles});
set(container,'Parent',handles.uipanel);
set(container,'Units','normalized','Position',[0 0 1 1]);
tree.setMultipleSelectionEnabled(true);
handles.tree=tree;
handles.regions=regions{1};

% % try the checkbox version
% jtree = tree.getTree;
% jCheckBoxTree = CheckBoxTree(jtree.getModel);
% jScrollPane = com.mathworks.mwswing.MJScrollPane(jCheckBoxTree);
% [jComp,hc] = javacomponent(jScrollPane,[150,10,120,110],gcf);

% try another thing; works quite nice

jtree = handle(tree.getTree,'CallbackProperties');
set(jtree, 'MousePressedCallback', {@mousePressedCallback, jtree, tree, handles});

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ms_gui_mouseAtlas wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ms_gui_mouseAtlas_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function mySelectFcn(tree,value,handles)
jtree = tree.getTree;
% path = jtree.getSelectionPath;
node=jtree.getSelectionPath.getLastPathComponent;
p=findall(handles.axes1,'Tag',strrep(char(node.getName()),handles.shape,''));
if ~isempty(p)
    handles.trans_slider.Value = p.FaceAlpha;
end

function nodes=build_tree(tree,value,reg)
% for more info: http://undocumentedmatlab.com/blog/customizing-uitree-nodes
fprintf('try to get childs for %d\n',value)
tarReg=getRegionById( value, reg{1} );
if ~isempty(tarReg)
    for ix=1:length(tarReg.children)
        if isempty(tarReg.children{ix}.children)
            nodes(ix)=uitreenode('v0',tarReg.children{ix}.id,[tarReg.children{ix}.name],[],true);
        else
            nodes(ix)=uitreenode('v0',tarReg.children{ix}.id,[tarReg.children{ix}.name],[],false);
        end
        rgb=reshape(sscanf(tarReg.children{ix}.color_hex_triplet.','%2x'),3,[]);
        I = uint8(zeros([16 16 3])); I(:,:,1) = rgb(1); I(:,:,2) = rgb(2); I(:,:,3) = rgb(3);
        nodes(ix).setIcon(im2java(I));
    end
end

function [ target ] = getRegionById( id, reg )
% finds a region by its ID
% reg is the 'big structure' containing all regions
fprintf('trying to find region by its ID.. ');
child=1; tar=0;
target='';
while child
    for i=1:length(reg)
        if reg(i).id==id 
%             fprintf('found id\n'); 
            target=reg(i); 
            tar=1; 
            break; 
        end
    end
    if tar; break; end;
    addreg=[];
    for i=1:length(reg)
        if ~isempty(reg(i).children)
              addreg=[addreg reg(i).children{:}];
        end
    end
    reg=addreg;
    for i=1:length(reg)
        if reg(i).id==id 
%             fprintf('found id\n'); 
            target=reg(i); 
            tar=1; 
            break; 
        end
    end
    child=0;
    for i=1:length(reg)
        if ~isempty(reg(i).children)
            child =1;
        end
    end
end
fprintf('found\n')


% --- Executes on button press in pushbutton_update.
function pushbutton_update_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text_updateROIs,'Visible','On'); pause(0.01)
nodes=ms_getSelectedNodes(handles);
plotNodes(nodes,handles,hObject);

function nodes=ms_getSelectedNodes(handles)
tree=handles.tree;
%nodes=tree.getSelectedNodes;
root=tree.getRoot;
nodetmp = root.getNextNode;
ix=1;
while ~isempty(nodetmp)
%     disp(nodetmp.getNextNode)
    if strfind(char(nodetmp.getName),handles.shape)
        nodes(ix)=nodetmp;
        ix =ix+1;
    end
    nodetmp = nodetmp.getNextNode; % this is important for the "parent separated by child" thing! 
end

function plotNodes(nodes,handles,hObject)
% get the axes and delete all child
h=handles.axes1;
% for ix=1:length(nodes)
%     toPreserve{ix}=findall(h,'Type','Patch','Tag',strrep(char(nodes(ix).getName()),handles.shape,''));
% end
allPatches = findall(h,'Type','Patch');
if ~isempty(allPatches)
    PatchesToKeep = zeros(length(allPatches),1);
    for ix=1:length(nodes)
        PatchesToKeep = PatchesToKeep + ~cellfun(@isempty, strfind(cellstr(char(allPatches.Tag)),strrep(char(nodes(ix).getName()),handles.shape,'')));
        allNodesNames{ix}=strrep(char(nodes(ix).getName()),handles.shape,'');
    end
    if isfield(handles,'ForceUpdate'); PatchesToKeep = zeros(length(allPatches),1); handles=rmfield(handles,'ForceUpdate'); end % quick and dirty
    delete(allPatches(~PatchesToKeep)); allPatches(~PatchesToKeep)=[];
    NodesToRemove=zeros(length(nodes),1); % we cannot delete nodes!!
    for ix=1:length(nodes)
        if sum(~cellfun(@isempty, strfind(cellstr(char(allPatches.Tag)),allNodesNames{ix})))
%         NodesToRemove = NodesToRemove + ~cellfun(@isempty, strfind(cellstr(char(allPatches.Tag)),allNodesNames{ix}));
            NodesToRemove(ix)=1;
        end
    end
    tmp=nodes; clear nodes; nodes=tmp(~NodesToRemove);
end

% delete(allchild(h));
% camlight; lighting(h, 'gouraud');

delete(findall(h,'Type','Light'))
% plot a brain if needed
if handles.checkbox_showBrain.Value
    if isempty(findall(h,'Type','Patch','Tag','brain'))
        fprintf('plotting brain.. ')
        fvbrain=isosurface(smooth3(padarray(handles.brain,[2 2 2]),'box',3),0.2);
        % get it into voxelspace
        tmp=handles.RotMat * [fvbrain.vertices'; ones(1,size(fvbrain.vertices',2))];
        tmp(end,:)=[]; fvbrain.vertices=tmp'; clear tmp
        
        hpatchBrain=patch(fvbrain,'Parent',h);
        hpatchBrain.Tag='brain';
        reducepatch(hpatchBrain,0.1,'fast')
        set(hpatchBrain,'FaceColor',[0.7 0.7 0.7])
        set(hpatchBrain,'EdgeColor','none')
        set(hpatchBrain,'FaceAlpha',0.2);
        fprintf('done\n')
    end
end

if ~isempty(nodes)
% prepare the ROIs
[mtx, tarReg] = getMtxAndTarRegs(nodes, handles);
% here we 'save' the mtx and tarRegs before we do the processing stuff for
% plotting
handles.mtx=mtx;
handles.tarRegs=tarReg;

% create isosurface
% parfor ix=1:length(tarReg)
for ix=1:length(tarReg)
    p=findall(handles.axes1,'Tag',tarReg(ix).name);
    if isempty(p)
        fprintf('create isosurface for %s.. ',tarReg(ix).name);
        tic
        mtx{ix}(mtx{ix}>0)=1;
        mtx{ix}=smooth3(mtx{ix},'box',3); %smoothing just looks nice
        fv{ix}=isosurface(mtx{ix});
        fprintf('done\n');
        toc
    end
end
if exist('fv','var')
ToRem=~logical(size(fv)); % if a region is too small for illustration
for ix=1:length(fv)
    % transfer to worldspace
    try
        tmp=handles.RotMat * [fv{ix}.vertices'; ones(1,size(fv{ix}.vertices',2))];
        tmp(end,:)=[]; fv{ix}.vertices=tmp'; clear tmp
    catch
        warndlg(sprintf('Region %s is too small for illustration!\nBut it will be in the final atlas.',tarReg(ix).name));
%         (fprintf('Region %s is too small for illustration!\nBut it will be in the final atlas.',tarReg(ix).name));
        ToRem(ix)=1; %warnflag=1;
    end
end
% remove too small regions (if needed)
fv(ToRem) =[]; tarReg(ToRem)=[];
% patch the ROIs
% for ix=1:length(tarReg)
for ix=1:length(fv)
    fprintf('create patch for %s.. ',tarReg(ix).name);
%     handles.text_updateROIs.String = sprintf('Patching %s..',tarReg(ix).name); pause(0.01)
    tic
    rgb=reshape(sscanf(tarReg(ix).color_hex_triplet.','%2x'),3,[]).'/255;
    hpatch{ix}=patch(fv{ix},'Parent',h, 'FaceColor',rgb,'EdgeColor','none');
    hpatch{ix}.Tag = tarReg(ix).name;
    reducepatch(hpatch{ix},0.1,'fast') % that makes patch less performance hungry
%     isonormals(mtx{ix},hpatch{ix});
    fprintf('done\n');
    toc
end
end
end
axes(h); camlight; lighting(h, 'gouraud');
axis equal
h.XLabel.String='X'; h.YLabel.String='Y'; h.ZLabel.String='Z';
fprintf('plot complete\n');
set(handles.text_updateROIs,'Visible','Off');
handles.text_updateROIs.String = 'Updating ROIs..';
if exist('warnflag','var')
    warndlg('At least one selected ROI had to be removed! See command window')
end

guidata(hObject,handles);

function [mtx, tarReg] = getMtxAndTarRegs(nodes, handles)

for ix=1:length(nodes) % get the selected nodes
%     handles.text_updateROIs.String = sprintf('Working on %s..',strrep(char(nodes(ix).getName), handles.shape,'')); pause(0.01)
    id(ix)=nodes(ix).getValue;
    tarReg(ix)=getRegionById( id(ix), handles.regions); % find the selected regions by their ID
    tarRegList{ix}=get_list_subregionsNum(tarReg(ix)); % get the all needed IDs -> that means the number ID of every child and the selected region itself!
    tarRegList{ix}(isnan(tarRegList{ix}))=[];
end
for ix=1:length(tarReg)
    tarReg(ix).tarRegList=tarRegList{ix};
end
for ix=1:length(tarReg) % create the actual matrices for later visualization
    mtx{ix}=zeros(size(handles.Mtx));
    fprintf('create matrix for %s...',tarReg(ix).name);
    for iy=1:length( tarReg(ix).tarRegList)
        mtx{ix}(handles.Mtx== tarReg(ix).tarRegList(iy))=tarReg(ix).id;
    end
    if isempty(find(mtx{ix},1))
        warning(sprintf('\nCould not find any entries for %s. \n ROI will not be taken into account!',tarReg(ix).name))
%         warnflag=1;
        warndlg('At least one selected ROI had to be removed! See command window')
    end
    fprintf(' done\n');
end

% control for "empty" mtx and remove them
rem=logical(false(1,length(tarReg)));
for ix=1:length(tarReg)
    if isempty(find(mtx{ix},1))
        rem(ix)=1;
    end
end
mtx(rem)=[]; tarReg(rem)=[]; clear rem; clear iy;
% catch the case that nothing is to plot
if isempty(tarReg)
    warndlg('There is nothing to plot!')
    return
end

% take care about the hemispheres
for ix=1:length(mtx)
    if handles.HemiLeft_radiobutton.Value
        mtx{ix}(:,:,1:size(mtx{ix},3)/2)=0; % left hemisphere
    elseif handles.HemiRight_radiobutton.Value
        mtx{ix}(:,:,size(mtx{ix},3)/2+1:end)=0; % right hemisphere
    end
end



function numList=get_list_subregionsNum(reg)
%returns a list of numbers of the lowest regions level which are within the
%given region 'reg'

%!! reg has to be a structure! !!
fprintf('get all children of region %s..', reg.name')
child=1;
n=1;
while child
    for i=1:length(reg)
        if ~isempty(reg(i).children)
            reg(end+1:end+length(reg(i).children))=[reg(i).children{:}];
            tmplist(n)=reg(i).id; n=n+1; % it seems that the lowest "region" isn't always the end of annotation!! e.g. hypothalamus -> so in the end we need all numbers!
            reg(i)=[];
        end
    end
    child=0;
    for i=1:length(reg)
        if ~isempty(reg(i).children)
            child =1;
        end
    end
end

numList=[reg.id];
if exist('tmplist','var')
    numList=[numList unique(tmplist)];
end
fprintf('done\n');


% --- Executes on button press in pushbutton_createAtlas.
function pushbutton_createAtlas_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_createAtlas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'SavePath')
    [FileName, PathName] = uiputfile('*.nii', 'Save Atlas', handles.SavePath);
else
    [FileName, PathName] =uiputfile('*.nii', 'Save Atlas');
end
AtlasName = [PathName FileName];
if ~isnumeric(PathName)
    handles.SavePath = PathName;
    h=waitbar(0,['create ' FileName]);
    % we load all needed again; so we don't have to update/reload everything!
    nodes=ms_getSelectedNodes(handles);
    [mtx, tarReg] = getMtxAndTarRegs(nodes, handles);
    handles.mtx=mtx; handles.tarRegs=tarReg;
    
    disp(['create atlas.. ' AtlasName]);
    Vout=handles.Vatlas;
    Vout.fname=AtlasName;
%     Mtx=zeros(size(handles.mtx{1}));
    Mtx=nan(size(handles.mtx{1})); % are NaNs better?
    for ix=1:length(handles.mtx) % here it is that the parent is separated by its (possibly) selected children; but just because we do this "nextNode" thing to find out which nodes are selected at all!!
        Mtx(handles.mtx{ix}>0)=max(unique(handles.mtx{ix}));
        waitbar(ix/length(handles.mtx),h)
    end
    if handles.BinaryMask_checkbox.Value
        Mtx(Mtx>0)=1;
    end
    if handles.sepHem_radiobutton.Value
        warndlg('You have chosen seperate hemispheres. The resulting data cannot be loaded then. I recommend also saving an atlas with "both" hemispheres!');
        newMtx=nan(size(handles.mtx{1})); newMtxR=newMtx;
        for ix=1:length(handles.tarRegs)
            newMtx(Mtx==handles.tarRegs(ix).id)=ix;
            newMtxR(Mtx==handles.tarRegs(ix).id)=length(handles.tarRegs)+ix;
            newtarRegs(ix) = handles.tarRegs(ix); newtarRegs(ix).name=[handles.tarRegs(ix).name ' leftHem']; newtarRegs(ix).id=ix;
            newtarRegs(length(handles.tarRegs)+ix) = handles.tarRegs(ix); newtarRegs(length(handles.tarRegs)+ix).name=[handles.tarRegs(ix).name ' rightHem']; newtarRegs(length(handles.tarRegs)+ix).id=length(handles.tarRegs)+ix;
        end
        Mtx = newMtx; Mtx(:,:,size(Mtx,3)/2+1:end)=newMtxR(:,:,size(Mtx,3)/2+1:end);
        handles.tarRegs = newtarRegs;
    end
    V(1)=spm_write_vol(Vout,Mtx); % this is the new atlas
    
    spm_jobman('initcfg');
    % transform to Dorr in Pax space
    matlabbatch = ms_allen2pax_batch(V(1).fname);
    spm_jobman( 'run', matlabbatch );
    delete([PathName matlabbatch{1}.spm.util.reorient.prefix FileName]) % cleanup
    movefile([PathName 'p' matlabbatch{1}.spm.util.reorient.prefix FileName], [PathName strrep(FileName, '.nii', '_inPax.nii')])  
    
    % also provide the allen average image
    Voutavg=handles.Vavg;
    Voutavg.fname = [PathName Voutavg.fname];
    V(2)=spm_write_vol(Voutavg,handles.Mtxavg);
    % make some zero padding
    Vzp=spm_vol([PathName strrep(FileName, '.nii', '_inPax.nii')]);
    Mtxzp=spm_read_vols(Vzp); Mtxzp=padarray(Mtxzp,[8 8 8]); Vzp.dim=size(Mtxzp);
    Vzp.mat(:,end)=Vzp.mat(:,end)-diag(Vzp.mat).*[8 8 8 0]';
    spm_write_vol(Vzp, Mtxzp);
    
    % show the atlas and the average template
    try
        spm_check_registration(V, {FileName, 'Avg. Template'}) % setting the interpolation would be nice!
    catch
        disp('cannot open spm_check_registration, sorry!');
    end
    % create a txt file
    fid=fopen(strrep(AtlasName,'.nii','.txt'),'w');
    if handles.BinaryMask_checkbox.Value; fprintf(fid, '# the mask included these regions\n'); end
    t=struct2cell(handles.tarRegs); % questionable if we should sort
    [~,idx]=sort(t(5,1,:));
    % rename the tarregs to make clear if it is left or right
    for ix=1:length(handles.tarRegs)
        if handles.HemiLeft_radiobutton.Value
            handles.tarRegs(ix).name = [handles.tarRegs(ix).name ' leftHem']; % left hemisphere
        elseif handles.HemiRight_radiobutton.Value
            handles.tarRegs(ix).name = [handles.tarRegs(ix).name ' rightHem']; % right hemisphere
        end
    end
    for ix=1:length(handles.tarRegs)
        if ix < length(handles.tarRegs)
            fprintf(fid, '%s\t %d\n',strrep(handles.tarRegs(idx(ix)).name,' ','_'), handles.tarRegs(idx(ix)).id);
        else
            fprintf(fid, '%s\t %d',strrep(handles.tarRegs(idx(ix)).name,' ','_'), handles.tarRegs(idx(ix)).id);
        end
    end
    fclose(fid);
    fprintf('done\n');
    close(h)
    guidata(hObject, handles);
end


% --- Executes on button press in checkbox_showBrain.
function checkbox_showBrain_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_showBrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_showBrain

function mousePressedCallback(hTree, eventData, jtree, mtree, handles) %,additionalVar)
  % if eventData.isMetaDown % right-click is like a Meta-button
  % if eventData.getClickCount==2 % how to detect double clicks
  % there is a eventData.ControlDown bool!
%  get(eventData) 
  % Get the clicked node
    clickX = eventData.getX;
    clickY = eventData.getY;
    treePath = jtree.getPathForLocation(clickX, clickY);
    % check if a node was clicked
    if ~isempty(treePath) && eventData.isMetaDown
%         disp('something was selected')
        node = treePath.getLastPathComponent;
%         char(node.getName)
%         handles.shape = '<html> <font color="green"><b>';
        if strfind(char(node.getName), handles.shape)
            nam = char(node.getName);
            node.setName(strrep(nam,handles.shape,''));
        else
            node.setName([handles.shape char(node.getName)])
        end
        model=mtree.Model; model.nodeChanged(node); % that does the trick!
%         jtree.collapsePath(treePath)
%         jtree.expandPath(treePath)
%         mtree.reloadNode(node);
%         mtree.repaint;
%         jtree.treeDidChange();
    end


% --- Executes on button press in SetLight_pushbutton.
function SetLight_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to SetLight_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% camlight(handles.axes1, 'Position', 'headlight')
% l=light('Position',handles.axes1.CameraPosition, 'Style','local');
l=findobj(handles.axes1,'Type','Light');
l.Position=[handles.axes1.CameraPosition(1) handles.axes1.CameraPosition(2) handles.axes1.CameraPosition(3)];


% --------------------------------------------------------------------
function save_menu_Callback(hObject, eventdata, handles)
% hObject    handle to save_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% you cannot save a uitree! so we don't even try!

% --------------------------------------------------------------------
function load_menu_Callback(hObject, eventdata, handles)
% hObject    handle to load_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% we use the annotation txt file
if isfield(handles,'SavePath')
    [f, p] = uigetfile('*.txt', 'Pick a annotations txt-file', handles.SavePath);
else
    [f, p] = uigetfile('*.txt', 'Pick a annotations txt-file');
end
if ~isequal(f,0)
    handles.pushbutton_update.Enable='off';
    handles.pushbutton_createAtlas.Enable='off';
    fid = fopen([p f]);
    C = textscan(fid, '%s%f');
    fclose(fid);
    % create a new node tree
    delete(handles.uipanel.Children); regions=handles.regions;
    root=uitreenode('v0', (regions.id),[regions.name],[],false);
%     [tree, container]=uitree('v0','Root',root,'ExpandFcn',{@build_tree,{regions}});
    [tree, container]=uitree('v0','Root',root,'ExpandFcn',{@build_tree,{regions}},'SelectionChangeFcn',{@mySelectFcn, handles});
    set(container,'Parent',handles.uipanel);
    set(container,'Units','normalized','Position',[0 0 1 1]);
    tree.setMultipleSelectionEnabled(true);
    handles.tree=tree;
    jtree = handle(tree.getTree,'CallbackProperties');
    set(jtree, 'MousePressedCallback', {@mousePressedCallback, jtree, tree, handles});
    
    model=handles.tree.Model;
    ids = num2cell(C{2});
    for ix=1:length(ids)
        % maybe it works to read the txt file and set/expand everything
        [ target ] = getRegionById( ids{ix}, handles.regions ); clear wayUp;
        wayUp(1) = ids{ix};
        while ~isempty(target.parent_structure_id)
            wayUp(end+1)=target.parent_structure_id;
            target = getRegionById( wayUp(end), handles.regions );
        end
        wayUp = wrev(wayUp); wayUp = wayUp(2:end); % root is already given
        tree=handles.tree;
        root=tree.getRoot;
        handles.tree.expand(root) % sth we always have to do
        pause(0.1); % we need some time to expand that stuff
        tmp=root;
        endflag=0; % shape = '<html> <font color="green"><b>';
        while ~endflag % that expands everything till we are in the correct "branch"
            tmp=tmp.getNextNode;
            if length(wayUp) > 1
                if tmp.getValue == wayUp(1)
                    handles.tree.expand(tmp);
                    wayUp = wayUp(2:end); 
                    while ~tmp.getChildCount; pause(0.01); disp('tmp has no child'); end % works
%                     pause(0.1); % we need some time to expand that stuff
                end
            else
                if tmp.getValue == wayUp; 
%                     disp(tmp); 
                    pause(0.1); 
                    tmp.setName([handles.shape char(tmp.getName)]); 
%                     handles.tree.repaint; % not enough
                    model.nodeChanged(tmp); 
                    endflag=1;  
                end
            end
            
        end
    end
    guidata(hObject, handles);
    handles.pushbutton_update.Enable='on';
    handles.pushbutton_createAtlas.Enable='on';
end


% --- Executes on button press in BinaryMask_checkbox.
function BinaryMask_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to BinaryMask_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BinaryMask_checkbox


% --- Executes on slider movement.
function trans_slider_Callback(hObject, eventdata, handles)
% hObject    handle to trans_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% disp(get(hObject,'Value'))
tree=handles.tree;
jtree = tree.getTree;
paths = jtree.getSelectionPaths;
for ix=1:length(paths)
    node=paths(ix).getLastPathComponent;
    % node=tree.getSelectedNodes;
    % disp(node.getName());
    p=findall(handles.axes1,'Tag',strrep(char(node.getName()),handles.shape,''));
    if ~isempty(p)
        p.FaceAlpha=get(hObject,'Value');
    end
end
% --- Executes during object creation, after setting all properties.
function trans_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trans_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in HemiBoth_radiobutton.
function HemiBoth_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to HemiBoth_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ForceUpdate=1;
guidata(hObject,handles);


% --- Executes on button press in HemiLeft_radiobutton.
function HemiLeft_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to HemiLeft_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ForceUpdate=1;
guidata(hObject,handles);


% --- Executes on button press in HemiRight_radiobutton.
function HemiRight_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to HemiRight_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ForceUpdate=1;
guidata(hObject,handles);


% --- Executes on button press in sepHem_radiobutton.
function sepHem_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to sepHem_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.ForceUpdate=1;
guidata(hObject,handles);
