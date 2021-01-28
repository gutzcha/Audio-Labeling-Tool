function [...
    AllClustersForFile,...
    AudioFullName_tSNE_ClusterAnalysis,...
    AllClustersForFileSyllables,...
    AllClustersDurationForFile,TriggerTime,LablesMap] = prepareTabForDeepPheno(tab)
%Convert table to cell data variables
%{
Inputs:
        tab - an roi Table or audioClip object
Outputps:
        AllClustersForFile:
           - cell array with onset of syllable fragmets. each syllable is
           split to 6 ms segments.
        AudioFullName_tSNE_ClusterAnalysis:
           - Path to the audio files
        AllClustersForFileSyllables:
          - cell array with onset of sellables, each syllable with it's own
          length
        AllClustersDurationForFile
         - cell array with duration of each syllable in sec.

If there are syllables from multiple sources, this a cells of cells
%}

%Check input
if isa(tab, 'audioClip')
    typeList = tab.allROIs;
    TriggerTime = tab.triggerTime;
    tab = tab.roiTable;
elseif isa(tab,'table')
    typeList = audioClip.roiTypeList;
    TriggerTime = 0;
else
    error('Invalid input')
end

%Get list of unique paths, each path represents a new file.
AudioFullName_tSNE_ClusterAnalysis = unique(tab.SourcePath); %This is ready
numFiles = size(AudioFullName_tSNE_ClusterAnalysis,1);
numTypes = size(typeList,1);

%Initialize other outputs
AllClustersForFile = cell(numFiles,numTypes);
AllClustersForFileSyllables = cell(numFiles,numTypes);
AllClustersDurationForFile = cell(numFiles,numTypes);

%AllClustersForFile
%Split each syllable into 6 ms fragments
segLenOne = 6e-3; %each segment is 6 ms
%We must seperate between the lines by the file name>cluster(lable)

%Crearte containter - directory between labels and index
LablesMap = (typeList.Name)';
%Function to segment syllables syllables 
segmentSyllables= @(s,d,l) s:l:(s+d);
 
%First seperate into different files #TO DO
for ifile = 1:numFiles
    for itype = 1:numTypes
        thisType = typeList{itype,1};
        inds = tab.Label == thisType;
        start = tab{inds,"TimeStart"};
        duration = tab{inds,"Duration"};
        AllClustersForFileSyllables{ifile,itype} = start;
        AllClustersDurationForFile{ifile,itype} = duration;
        
        % segmentSyllables
        segLen = ones(size(start)).*segLenOne;
        f = arrayfun (segmentSyllables,start,duration,segLen,'UniformOutput',false);
        AllClustersForFile{ifile,itype} =  [f{:,:}]';
        
    end
end