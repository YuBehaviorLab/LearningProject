clear;clc;
rng(20220819); % Set the random seed for reproducibility of the results
rpath = pwd; % rootpath
%% Initiate
path_functions = fullfile(rpath,'0_Functions');addpath(path_functions);
path_archive = fullfile(rpath,'1_Archive');

GroupName = {'Lesion','Sham'};
PeriodName = {'Pre','Post'};
prdVar = {1:8,9:16};

% sessions selected for plotting (jy)
plotRange = 2:15;

% extract data again or not?
toextract = 1;

%% Extract Data
if toextract
    tarPath = {};
    for i=1:length(GroupName) % Get each subject's data folder & grouping info
        grpName = GroupName{i};
        grpPath = fullfile(rpath,grpName);
        grpDir = dir(grpPath); % each rat has a folder. read folders' names and get the rat's data
        for j=1:length(grpDir)
            if ~grpDir(j).isdir || strcmp(grpDir(j).name,'.') || ...
                    strcmp(grpDir(j).name, '..') || strcmp(grpDir(j).name,'0_shelve')
                continue;
            end
            sbjDir = fullfile(grpPath,grpDir(j).name);
            if isfolder(sbjDir)
                tmpInfo = {sbjDir,grpDir(j).name,grpName};
                tarPath(end+1,1:3) = tmpInfo;
            end
        end
    end
 
    disp(tarPath)
 
%     tarPath =
%     8×3 cell array
%     {'C:\Users\jiani\…'}    {'Hsiao-hsien'}    {'Lesion'}
%     {'C:\Users\jiani\…'}    {'ONeal'      }    {'Lesion'}
%     {'C:\Users\jiani\…'}    {'Risotto'    }    {'Lesion'}
%     {'C:\Users\jiani\…'}    {'Walter'     }    {'Lesion'}
%     {'C:\Users\jiani\…'}    {'Kieslowski' }    {'Sham'  }
%     {'C:\Users\jiani\…'}    {'Matias'     }    {'Sham'  }
%     {'C:\Users\jiani\…'}    {'Panini'     }    {'Sham'  }
%     {'C:\Users\jiani\…'}    {'Ulysess'    }    {'Sham'  }

    btAll2d = {}; % 'b'pod 't'able Allsessions 2d(for all subjects)
    for ipath=1:size(tarPath,1) % Extract the raw data & packaging
        dataPath = tarPath{ipath,1};
        grpName = tarPath{ipath,3};
        cd(dataPath);
        % extract and processing
        FileNames = arrayfun(@(x)x.name, dir('*DSRT*.mat'), 'UniformOutput', false); % each 'FileNames' is a session
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
        disp(FileNames)
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
%         first half pre, second half post (conditions are not labeled otherwise)
%         {'Hsiao-hsien_DSRT_06_3FPs_20220129_172638.mat'}
%         {'Hsiao-hsien_DSRT_06_3FPs_20220206_190958.mat'}
%         {'Hsiao-hsien_DSRT_06_3FPs_20220208_192320.mat'}
%         {'Hsiao-hsien_DSRT_06_3FPs_20220209_194527.mat'}
%         {'Hsiao-hsien_DSRT_06_3FPs_20220210_172708.mat'}
%         {'Hsiao-hsien_DSRT_06_3FPs_20220211_131813.mat'}
%         {'Hsiao-hsien_DSRT_06_3FPs_20220212_141957.mat'}
%         {'Hsiao-hsien_DSRT_06_3FPs_20220213_200854.mat'}

%         {'Hsiao-hsien_DSRT_06_3FPs_20220225_140032.mat'}
%         {'Hsiao-hsien_DSRT_06_3FPs_20220226_134509.mat'}
%         {'Hsiao-hsien_DSRT_06_3FPs_20220228_143040.mat'}
%         {'Hsiao-hsien_DSRT_06_3FPs_20220301_160559.mat'}
%         {'Hsiao-hsien_DSRT_06_3FPs_20220302_141907.mat'}
%         {'Hsiao-hsien_DSRT_06_3FPs_20220303_143602.mat'}
%         {'Hsiao-hsien_DSRT_06_3FPs_20220304_124042.mat'}
%         {'Hsiao-hsien_DSRT_06_3FPs_20220305_141157.mat'}

        btAll_raw = cell(1,length(FileNames));
        for i=1:length(FileNames)
            % get data from a single session 
            btAll_raw{i} = DSRT_DataExtract_Block(FileNames{i},false);
        end
        % merge multiple files of one day into one file
        btAll_merge = DSRT_DataMerge_Block(btAll_raw,2);

%         Each cell is a session, organized in a table 

%         17 is the number of variables. 
%         btAll_merge =
%         1×16 cell array
%         Columns 1 through 5
%         {304×17 table}    {203×17 table}    {279×17 table}    {278×17 table}    {297×17 table}
%         Columns 6 through 10
%         {252×17 table}    {250×17 table}    {260×17 table}    {170×17 table}    {242×17 table}
%         Columns 11 through 15
%         {262×17 table}    {281×17 table}    {288×17 table}    {275×17 table}    {285×17 table}
%         Column 16
%         {280×17 table}

%          Columns 1 through 6
%          {'Subject'}    {'Date'}    {'StartTime'}    {'Task'}    {'iTrial'}    {'BlockNum'}
%          Columns 7 through 12
%          {'TrialNum'}    {'TrialType'}    {'TimeElapsed'}    {'FP'}    {'RW'}    {'DarkTry'}
%          Columns 13 through 17
%          {'ConfuseNum'}    {'Outcome'}    {'HT'}    {'RT'}    {'MT'}


        % add group information to data file
        % add  eg "Lesion-Pre"
        btAll = cell(1,length(btAll_merge));
        for i=1:length(btAll_merge)
            bt = btAll_merge{i};
            for j=1:length(prdVar)
                if ismember(i,prdVar{j})
                    dateGrp = PeriodName{j}; % Pre/Post
                    newVars = repelem(strcat(string(grpName),'-',string(dateGrp)),size(bt,1))';
                    break;
                end
            end
            btAll{i} = addvars(bt,newVars,'After','Date','NewVariableNames','Group');
        end

        % Save
        savename = 'bmixedAll_' + upper(btAll{1}.Subject(1)); %tarPath{ipath,2}
        save(savename, 'btAll');
        btAll2d(end+1,1:length(btAll)) = btAll;
    end
    save(fullfile(path_archive,'bmixedAllsbj.mat'),'btAll2d','plotRange');
    task2d = cellfun(@(in) getGroup(in),btAll2d,'UniformOutput',false); % for check
end;
%% Plot
cd(path_archive);
load('bmixedAllsbj.mat');

plotLesionAfterLearningJ(btAll2d,plotRange);
%% Functions
function out = getGroup(in)
    if ~isempty(in)
        out = char(in.Group(1));
    else
        out = '';
    end
end