%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear


% ---------Specify FIS System--------------------------------------------
FIStype = uigetfile('*.fis','Load the fis file:');   % If you want to use a different FIS, specifiy it here 

data = readtable('C:\et_al\Shared\Projects\USA\CHaMP\ResearchProjects\HabitatSuitability\wrk_Data\FISValidation\ChinookSpawner\UGR_ValidationSites.csv');
fPath = table2array(data(:,{'visit_dir'}));

for ii = 1:length(fPath)
    fileName = dir(fullfile(fPath{ii},'**/coverIndex_Inputs.csv'));
    fileName = char(strcat({fileName.folder}, filesep, {fileName.name}));
    disp(fileName);
    FuzzyCoverIndex(fileName, FIStype);
end

% % ---------Specify FIS System--------------------------------------------
% FIStype = uigetfile('*.fis','Load the fis file:');   % If you want to use a different FIS, specifiy it here 
% 
%# build a list of file names with absolute path
% % % fPath = uigetdir('.', 'Select directory containing all input files');
% % % if fPath==0, error('no folder selected'), end
% % % fileNames = dir(fullfile(fPath,'**/coverIndex_Inputs.csv'));
% % % fileNames = strcat({fileNames.folder}, filesep, {fileNames.name});
% % % 
% % % %# process each file
% % % for ii = 1:length(fileNames)
% % %     disp(fileNames{ii});
% % %     %FuzzyCoverIndex(fileNames{ii}, FIStype);
% % % end
