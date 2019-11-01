%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear


% ---------Specify FIS System--------------------------------------------
FIStype = uigetfile('*.fis','Load the fis file:');   % If you want to use a different FIS, specifiy it here 

%# specity number of inputs
numInputs = 4;

visits = readtable('C:\et_al\Shared\Projects\USA\CHaMP\ResearchProjects\HabitatSuitability\wrk_Data\FISValidation\ChinookSpawner\UGR_ValidationSites.csv');
visits_sub = visits(visits.AveBFW > 10.0,:);
fPath = table2array(visits_sub(:,{'visit_dir'}));

for ii = 1:length(fPath)
    fileName = dir(fullfile(fPath{ii},'**/FuzzyHSI_Inputs.csv'));
    fileName = char(strcat({fileName.folder}, filesep, {fileName.name}));
    disp(fileName);
    FuzzyHSI_LargeSites(fileName, FIStype, numInputs);
end

% % ---------Specify FIS System--------------------------------------------
% FIStype = uigetfile('*.fis','Load the fis file:');   % If you want to use a different FIS, specifiy it here 
% 
% %# specity number of inputs
% numInputs = 4;
% 
% % build a list of file names with absolute path
% fPath = uigetdir('.', 'Select directory containing all input files');
% if fPath==0, error('no folder selected'), end
% fileNames = dir(fullfile(fPath,'**/FuzzyHSI_Inputs.csv'));
% fileNames = strcat({fileNames.folder}, filesep, {fileNames.name});
% 
% % process each file
% for ii = 1:length(fileNames)
%     disp(fileNames{ii});
%     FuzzyHSI_LargeSites(fileNames{ii}, FIStype, numInputs);
% end

