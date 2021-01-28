classdef audioClip < handle
    %AUDIOCLIP represent a segment of audiodata
    %   Detailed explanation goes here
    
    properties
        %         + info      - a struct containig info on the source with fields:
        %               * ratNum   - rat serial number/id
        %               * paradigm - experimental paradigm
        %               * dayNum   - day number of experiment
        info = struct('fileName','','filePath','','ratNum',0,'paradigm','','dayNum',0);
        
        vec = [];               % the audio vector, if this is empty, no plot is made
        times_ = [];
        fs = 2.5e5;             % sampling rate of the audio vector
        window = 512;           % fft window length in bins
        overlap = 128;          % fft overlap length in bins
        fft = 512;              % number of DFT point also bins
        threshold = -inf;       % threshold - minimal value for spectrogram
        clim = [];              % color map limits, auto if empty
        climMode = 'auto'       % if clim is on auto, clim will be disregarded
        stft = [];              % short-time Fourier transform of the input signal
        cyclicalF = [];         % cyclicalF   - cyclical frequencies
        timeVec = [];           % vector of time instants
        psd = [];               % psd vector
        %         roiList = [];           % table of roi elements and their types
        
        emptyFlag = true;        % empty flag, true if the object is empty
        roiTable
        triggerTime = 0;
        sourceType = 'wav'      %type of souce vector - wav or mat
        vecSource = '';        %source of the vector, file of workspace
        ylims;
        absTime = 0;            %If the clip was spliced from a bigger clip, this  property saves the time
    end
    properties (GetAccess = 'public', SetAccess = 'public')
        addedTypeList = table('Size',[0,3],'VariableNames',{'Name','Color','Comments'},'VariableTypes',{'string','double','string'});
    end
    properties (Constant)
        roiTypeList = cell2table(...
            {"USV",[1,0,1],"";...
            'Low',[1,1,0],"";...
            "AC",[0,1,1],"";...
            "USV-Auto",[0,1,0],"AutoDetected USV";...
            "Noise",[1,0,0],"Noise segments identified as USV";...
            "FN",[0,0,0],"Missed syllables";...
            "USV-DS",[1,1,0.8],"AutoDetected USV by DeepSqueak"},...
            'VariableNames',{'Name','Color','Comments'})
        
        roiTableTemplate = table('Size',[0,8],'VariableNames',{'Label','TimeStart','TimeEnd','FrLow','FrHigh','Duration',...
            'SourcePath','Comments'},'VariableTypes',{'string','double','double','double','double','double','string','string'});
    end
    properties(Dependent)
        summaryTable
        audioLen
        allROIs
        path
        times
        
    end
    
    methods
        
        function this = audioClip(S,dataIn)
            
            this.roiTable = audioClip.roiTableTemplate;
            
            if nargin==0
                %AUDIOCLIP Construct an instance of this class
                %   Detailed explanation goes here
                return
            else
                if isstruct(S)
                    this.fs            = S.fs           ;
                    this.roiTable      = S.roiTable     ;
                    this.times_        = S.times_       ;
                    this.triggerTime   = S.triggerTime    ;
                    this.ylims         = S.ylims          ;
                    this.absTime       = S.absTime        ;
                    this.sourceType    = S.sourceType     ;
                    this.clim          = S.clim           ;
                    this.climMode      = S.climMode       ;
                    this.window        = S.window       ;
                    this.fft           = S.fft          ;
                    this.overlap       = S.overlap      ;
                    this.addedTypeList = S.addedTypeList;
                    this.vec           = S.vec            ;
                    this.info          = S.info         ;
                else
                    if isvector(S)
                        this.vec = S;
                    end
                    
                    if exist('dataIn','var')&&~isempty(dataIn)
                        this.info = dataIn;
                    end
                end
                this.emptyFlag = false;
            end
        end
        
        function set.climMode(this,val)
            %Set mode of clim
            %Clim must be a char or string array with two possible values:
            %'auto' or 'manual'
            if ischar(val)||isstring(val)
                
                switch val
                    case {'auto','manual'}
                        this.climMode = val;
                    otherwise
                        error('clim mode must be auto or manual')
                end
            else
                error('invalid clim mode')
            end
        end
        function len = get.audioLen(this)
            
            if isempty(this.vec)
                len = 0;
            else
                len = numel(this.vec)/this.fs;
            end
        end
        function [obj,sh,psd] = drawSpectrogram(obj,axSpect, axP)
            %Draw a spectrogram and a normlized pds
            %   Draw plots and update plots accordingly
            
            %{
Inputs:
    axSpect - handel for axis on which to plot the spectrogram
    axP     - handel for axis on which to plot the power spectral density
    dataIn  - an instanse of audioClip class
Outputs:
    sh        - handel for surf plot of spectrogram
    psd       - handel for the line plot psd
    dataOut   - an instanse of audioClip class, updated
       
    
            %}
            if nargin==1
                % Create figure
                figure1 = figure;
                
                % Create axes
                axSpect = axes('Parent',figure1,...
                    'Position',[0.13 0.3 0.775 0.7]);
                
                % Set the remaining axes properties
                set(axSpect,'XTick',zeros(1,0));
                % Create axes
                axP = axes('Parent',figure1,...
                    'Position',[0.13 0.07 0.775 0.23]);
            end
            
            if isempty(obj.vec)
                error ('Empty audio vector')
            end
            if obj.times(1)<0||obj.times(2)>obj.audioLen
                error('Time range must be between %f and %f',0,obj.audioLen)
            end
            vectorInds = (1+obj.fs*obj.times(1)):(obj.fs*obj.times(2));
            %Perform fft on audio vector
            if ~isempty(obj.vec)
                [s,f,t,ps] = spectrogram(obj.vec(vectorInds),obj.window,obj.overlap,obj.fft,obj.fs,'yaxis','MinThreshold',obj.threshold);
            else
                sh = [];
                psd = [];
                
                return
            end
            
            %Create plot with set threshold and clim
            sh = surf(axSpect,log10(abs(s)),'EdgeColor','none');
            view(axSpect,2);
            if ~isempty(obj.clim)
                axSpect.CLim = obj.clim;
            end
            
            %Create psd plot
            psSum = smooth(sum(ps),0.01);
            psSum = psSum./max(psSum);
            psd = plot(axP,t,psSum);
            
            %Make axes tight
            axis(axSpect,'tight')
            axis(axP,'tight')
            
            %Update output values
            obj.stft = s;
            obj.cyclicalF = f;
            obj.timeVec = t;
            obj.psd = ps;
        end
        function tab = get.allROIs(obj)
            added = obj.addedTypeList;
            basicTypes = obj.roiTypeList;
            if ~isempty(added)
                tab = [basicTypes;added];
            else
                tab = basicTypes;
            end
            
        end
        function ret = get.times(obj)
            maxLen = inf;
            if isempty(obj.times_)
                ret = [0,min(obj.audioLen,10)];
            else
                tempRet = obj.times_;
                ret(1) = max(0,tempRet(1));
                ret(2) = min(tempRet(2),obj.audioLen);
                %Make sure the length doen't go over maxLen seconds
                ret(2) = min(ret(2),ret(1)+maxLen);
            end
        end
        function set.times_(obj,t)
            %             maxDisp = 20;
            if isempty(t)
                return
            end
            
            if isnumeric(t)&&(all(size(t) == [1,2]))
                %                 if t(end)<=obj.audioLen&&t(1)>=0&&diff(t)<=maxDisp %#ok
                %                     obj.times_ = t;
                %                 end
                %             else
                %                 error('Invalid time range')
                obj.times_ = t;
            end
        end
        function set.roiTable(obj,input)
            try
                [obj.roiTableTemplate;input]; %#ok
            catch
                error('Invalid roi table input')
            end
            obj.roiTable = input;
        end
        function ret = isempty(obj)
            try
                ret = obj.emptyFlag ;
            catch
                ret = true;
            end
        end
        function set.ylims(obj,vals)
            if all(size(vals)~=[1,2])||~isnumeric(vals)
%                 warning('ylims must be a 1x2 vector')
                obj.ylims = [0,obj.fs/2]; %#ok
            else
                obj.ylims = vals;
            end
            
        end
        function ret = get.path(obj)
            ret = fullfile(obj.info.filePath,obj.info.fileName);
        end
        function set.addedTypeList(app,input)
            t = app.addedTypeList;
            if istable(input)||iscell(input)
                %If the new ori types allready exist in the variable,
                %replace them with new comment and new color
                newNames = input{:,1};
                oldNames = t{:,1};
                [~,ia,ib] = union(newNames,oldNames);
                t = [t(ib,:);input(ia,:)];
                app.addedTypeList = t;
            else
                error('Invalid input')
            end
            
        end
        function tab = roi2table(obj,roilist)
            
            if ~exist('roilist','var')||isempty(roilist)
                tab = obj.roiTableTemplate;
                return
            end
            
            if isa(roilist,'images.roi.Rectangle')
                %                 if isvalid(roilist)
                %
                %                 end
                Label = string(get(roilist,'Label'));
                positions = get(roilist,'Position');
                if iscell(positions)
                    positions = cell2mat(positions);
                end
            elseif isstruct(roilist)
                Label = roilist.Label;
                positions = roilist.Positions;
            else
                tab = obj.roiTableTemplate;
                return
            end
            
            
            TimeStart = positions(:,1);
            TimeEnd = TimeStart + positions(:,3);
            FrLow = positions(:,2);
            FrHigh = FrLow + positions(:,4);
            Duration = TimeEnd - TimeStart;
            %
            %I removed audioVector as a variable
            %             ind = @(s,e,fs,vec) vec(round((s*fs+1):(e*fs)),1);
            %             fs_ = obj.fs;
            %             vec_ = obj.vec;
            %             AudioVector = cellfun(ind, num2cell(TimeStart),num2cell(TimeEnd),repmat({fs_},size(TimeStart)),repmat({vec_},size(TimeStart)),'UniformOutput',false);
            %
            SourcePath = repmat({obj.path},size(TimeStart));
            Comments = repmat("",size(TimeStart));
            %             tab = table(Label,TimeStart,TimeEnd,FrLow,FrHigh,Duration,AudioVector,SourcePath,Comments);
            tab = table(Label,TimeStart,TimeEnd,FrLow,FrHigh,Duration,SourcePath,Comments);
            tab = sortrows(tab,{'TimeStart','Label'});
        end
        function obj = replaceRowsInTable(obj,tabIn,deleteAllFlag)
            if ~isa(tabIn,'table')
                tabIn = obj.roi2table(tabIn);
            end
            if ~exist('deleteAllFlag','var')
                deleteAllFlag = false;
                
            end
            if deleteAllFlag
                lables = unique(obj.roiTable.Label);
            else
                lables = unique(tabIn.Label);
            end
            
            obj = obj.deleteRowsFromTableByRange(obj.times,lables);
            obj = obj.addRowsToTable(tabIn);
            %             end
        end
        function obj = removeRowsFromTable(obj,lablesToRemove,rangesToRemove)
            %Inputs : lables - labels of rows to remove
            %         ranges - time ranges of rows to remove
            szLabels = size(lablesToRemove,1);
            szRanges = size(rangesToRemove,1);
            if szLabels~=szRanges
                error('Sizes of lables and ranges must be equel')
            end
            if isempty(obj.roiTable)
                return
            end
            
            lablesInObj = obj.roiTable.Label;
            rangesInObj = obj.roiTable{:,{'TimeStart','TimeEnd'}};
            
            indLabels = ismember(lablesInObj, lablesToRemove); %Find inds of tags
            indStarts = ismember(rangesInObj(:,1), rangesToRemove(:,1));
            indEndingss = ismember(rangesInObj(:,2), rangesToRemove(:,2));
            
            %Remove rows where all identifiers are the same
            indToRemove = all([indLabels,indStarts,indEndingss]);
            obj.roiTable(indToRemove,:) = [];
        end
        
        function obj = addRowsToTable(obj,tabIn)
            %Inputs : tabIn - table with data or roi object
            if isempty(tabIn)
                return
            end
            
            if isa(tabIn,'images.roi.Rectangle')||isa(tabIn,'matlab.graphics.GraphicsPlaceholder')
                tabIn = obj.roi2table(tabIn);
            end
            if isempty(obj.roiTable)
                obj.roiTable = tabIn;
            else
                obj.roiTable = [obj.roiTable;tabIn];
                obj.roiTable = sortrows(obj.roiTable,{'TimeStart','Label'});
            end
            
        end
        function obj = deleteRowsFromTableByRange(obj,rangeToDelete,lablesToDelete)
            
            %There is nothing to delete
            if isempty(lablesToDelete)
                return
            end
            
            %There is nothing to delete
            if isempty(obj.roiTable)
                return;
            end
            
            %Out of range
            if ~isempty(rangeToDelete) %If the range is empty, delete every thing
            if rangeToDelete(1)<0||rangeToDelete(2)>obj.audioLen
                return
            end
            else
               rangeToDelete = [0 obj.audioLen] ;
            end
            
            tab = obj.roiTable;
            %To delet a row is has to have the correct lable and be in the time range
            indsToDelete =  find(any(tab.Label == lablesToDelete',2)&(tab.TimeStart>rangeToDelete(1) & tab.TimeEnd<rangeToDelete(2)));
            
            
            if ~isempty(indsToDelete)
                tab(indsToDelete,:)=[];
            end
            obj.roiTable = tab;
            
        end
        function explore(obj)
            %If the is no source data, remmber that this was loaded from
            %workspace
            if isempty(obj.path)
                obj.vecSource = 'workspace';
                obj.sourceType = 'mat';
            end
            audioBrowser(obj,audioClip);
        end
        
        function [newVecOut, newFsOut] = playSegment(obj,fromFs,toFs,volume,speed,timeRange,vecIn)
           % playSegment(obj,fromFs,toFs,volume,speed,timeRange)
           
            if ~exist('fromFs','var')||isempty(fromFs)
                fromFs = 60000;
            end
            
            if ~exist('toFs','var')||isempty(toFs)
               toFs = 5000; 
            end
            
           
            
            if ~exist('volume','var')||isempty(volume)
                volume = 1;
            end
            
             if ~exist('speed','var')||isempty(speed)
               speed = 1; 
            end
            
             if ~exist('timeRange','var')||isempty(timeRange)
                timeRange = obj.times;
             end
                          
             if ~exist('vecIn','var')||isempty(vecIn)
                vecSegment = obj.vec(round((1+timeRange(1)*obj.fs):(timeRange(2)*obj.fs)));
             else
                 vecSegment = vecIn;
             end
             vecSegment = vecSegment./max(abs(vecSegment));
            fsratio = round(fromFs/toFs);
            lengthRatio = ceil(fsratio*speed);
            
            newFs = obj.fs/fsratio;
            newVec = (downsample(vecSegment,lengthRatio)).*volume;     %Down sample the vector so it will have the desired speed and multiply by volume       
            newVec = highpass(newVec,500,newFs);
            if nargout==0
                if newFs>100000
                    disp('Reducing sampling rate')
                    downSampleRatio = 5;
                    newVec = (downsample(vecSegment,downSampleRatio)).*volume;     %Down sample but dont change pitch
                    newFs = ceil(newFs/downSampleRatio);
                end
                soundview(newVec,newFs)
                %sound(newVec,newFs)
            else
                newVecOut = newVec;
                newFsOut = newFs;
            end
        end
        function vecseg = vecSegment(obj,timeRange)
            if ~exist('timeRange','var')
               timeRange = obj.times; 
            end
            vecseg = obj.vec(round((1+obj.fs*timeRange(1)):(obj.fs*timeRange(2))));
        end
        function edit(app)
            splitterEditorApp(app) ;
        end
        function h = editROIs(obj)
            h = editRoisApp(obj);
        end
        function obj1 = plus(obj1,obj2)
           %The must have the same path!
           if strcmp(obj1.path,obj2.path)
               T2 = obj2.roiTable;
               obj1.addedTypeList = obj2.addedTypeList;
               obj1.addRowsToTable(T2);
           end
        end
        function s = saveobj(this)
            s.info          = this.info         ;
            s.vec           = []                ;
            s.times_        = this.times_       ;
            s.fs            = this.fs           ;
            s.window        = this.window       ;
            s.overlap       = this.overlap      ;
            s.fft           = this.fft          ;
            s.threshold     = this.threshold    ;
            s.clim          = this.clim         ;
            s.climMode      = this.climMode     ;
            s.path          = this.path         ;
            s.stft          = []                ;
            s.cyclica       = []                ;
            s.timeVec       = []                ;
            s.psd           = []                ;
            
            s.addedTypeList = this.addedTypeList;
            s.emptyFlag     = this.emptyFlag    ;
            s.roiTable      = this.roiTable     ;
            s.sourceType    = this.sourceType   ;
            s.ylims         = this.ylims        ;
            s.absTime       = this.absTime      ;
            s.triggerTime   = this.triggerTime  ;
            
        end
    end
    methods(Static)
        function this = loadobj(S,flags)
            if ~exist('flags','var')
                flags = true;
            end
            if isstruct(S)
                path = S.path;
                if ~isfield(S,'sourceType')
                    S.sourceType = 'wav' ;
                end
                if ~isfield(S,'ylims')
                    S.ylims = [];
                end
                
                if ~isfield(S,'absTime')
                    S.absTime = 0;
                end
                
                if ~isfield(S,'clim')
                    S.clim = [];
                end
                
                if ~isfield(S,'climMode')
                    
                    S.climMode = 'auto';
                end
                if ~isfield(S,'triggerTime')
                    
                    S.triggerTime = 0;
                end
                
                
                %Make sure file name has right extention
                [~,fileName,ext] = fileparts(path);
                if isempty(ext)
                    ext = S.sourceType ;
                    path = [path,'.',ext];
                else
                    if ~strcmp(ext(2:end),S.sourceType)
                        warning('File extention (%s) and sourceType (%s) are not the same,\nChaging to extention',ext,S.sourceType)
                         S.sourceType = ext(2:end);
                    end
                end
                
                %Search for the file with full path, if it cant be found,
                %try to look for the file in current folder
                
                if ~isfile(path)
                    
                    path = [pwd,'\',fileName,ext];
                    if isfile(path)
                           info = getRecodringDetails([fileName,'.',ext],pwd,fileName);
                           S.info = info;
                    end
                end
                try
                    switch S.sourceType
                        
                        case 'wav'
                            vec = audioread(path) ;
                        case 'mat'
                            readfile = matfile(path);
                            classesInVar = whos(readfile);
                            indDouble = find(strcmpi('double',{classesInVar.class}),1,'first');
                            vec  = eval(['readfile.',classesInVar(indDouble).name]);
                    end
                    S.vec = vec;
                    this         = audioClip(S);
                    
                    
                catch
                    if flags
                        answer = questdlg('Unable to find file on path, please choose how to proceed', ...
                            'Load audio file', ...
                            'Load Manualy','Load without audio vector','Cancel','Cancel');
                        
                        %For some reason, this doesn't work with 4 options
                        %so I removed load roitable option
                        %                           answer = questdlg('Unable to find file on path, please choose how to proceed', ...
                        %                             'Load audio file', ...
                        %                             'Load Manualy','Load only roiTable','Cancel','Cancel');
                        %
                        % Handle response
                        switch answer
                            case 'Load Manualy'
                                [file,path] = uigetfile(...
                                    {'*.mat;*.wav','All Files (*.mat,*.wav)';...
                                    '*.wav','Audio File (.*wav)';...
                                    '*.mat','MATLAB File (*.mat)'});
                                
                                if file==0
                                    warning('No file was chocen, returning an empty audioClip file')
                                    this = audioClip;
                                    return
                                else
                                    [~,f,type]=fileparts(file);
                                    switch type
                                        case '.wav'
                                            vec = audioread(fullfile(path,file)) ;
                                        case '.mat'
                                            readfile = matfile(fullfile(path,file));
                                            classesInVar = whos(readfile);
                                            indDouble = find(strcmpi('double',{classesInVar.class}),1,'first');
                                            vec  = eval(['readfile.',classesInVar(indDouble).name]);
                                    end
                                    
                                    info = getRecodringDetails(file,path,f);
                                    S.info = info;
                                    S.vec = vec;
                                    this = audioClip(S);
                                    return
                                end
                                
                            case 'Load only roiTable'
                                this = S.roiTable;
                            case 'Load without audio vector'
                                S.vec = [];
                                this = audioClip(S);
                            case 'Cancel'
                                this = audioClip;
                                return
                        end
                    else
                        warning('Unable to find file, returning an empty object')
                        this = audioClip;
                        return
                    end
                end
            end
        end
    end
end

