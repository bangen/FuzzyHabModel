%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       FUZZY INFERENCE SYSTEM for Fish Cover
%        for use with Fuzzy Logic Toolbox 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This script uses a fuzzy inference system to estimate instream fish 
% cover index. The inputs are derived from a Champ crew survey.
% The script reqires inputs of LWD Count, uncdercut distance and 
% deep pool distance in CSV input file.  
%
% User Inputs: A CSV file with a header row (ignored) and 1st column is
% easting (x), 2nd column is northing (y), 3rd column is LWD Count, 4th
% column is undercut distance, 5th column is deep pool distance.

% Parameters:  None
% Outputs:     The program writes a new CSV file that is the same as the input,
% with an added column of fuzzy cover index (between 0 and 1).

function FuzzyCoverIndex(filepath, FIStype)

% Call up the *.fis file produced from the fuzzy logic toolbox
    % This will only run if the Fuzzy Logic Toolbox is loaded
    % The fuzzy rule system can be modified using the toolbox
	aFIS = readfis(FIStype);  

% ---------Read in CSV File of Inputs-------------------------------------

% Set output pathname and filename
outname = regexp(FIStype, '[.]', 'split');
outname = outname(1,1);
filepathStr = filepath;
[pathstr,name,ext] = fileparts(filepathStr);
%outpathname = strrep(pathstr, 'Inputs', 'Output');
%outfilename = strcat(outpathname, '\', outname, '.csv');
outfilename = strcat(pathstr, '\', outname, '.csv');
outfilename = outfilename{1};
disp(outfilename);

%% Use csvread to get rest of the data
data = csvread(filepath,2,0);
x=data(:,1);
y=data(:,2);
LWCount=data(:,3);
UCDist=data(:,4);
PoolDist=data(:,5);

%% Clear temporary variables
clearvars data raw R columnIndices;

disp('Done reading input file.');
    
%%------Get Rid of NonZero-----------------------------------------------    

% Find addresses of rows that meet condition
% Update array 
% Clear addresses after each find, update
    
    % Find LWCounts >= 0
    clear addF;
    addF = find(LWCount >= 0);
    x = x(addF);
    y = y(addF);  
    LWCount = LWCount(addF);
	UCDist = UCDist(addF);
	PoolDist = PoolDist(addF);
    
    % Find UCDists >= 0
    clear addF;
    addF = find(UCDist >= 0);
    x = x(addF);
    y = y(addF);  
    LWCount = LWCount(addF);
	UCDist = UCDist(addF);
	PoolDist = PoolDist(addF);
    
    % Find PoolDist >= 0
    clear addF;
    addF = find(PoolDist >= 0);
    x = x(addF);
    y = y(addF);  
    LWCount = LWCount(addF);
	UCDist = UCDist(addF);
	PoolDist = PoolDist(addF);

    % Find input values > fis max
    % Set value to just below fis max so that an error is 
    
    % Find LWCount above Max
    clear addF;
    addF = LWCount >= 400;
    LWCount(addF) = 399;
    
    % Find UCDists above Max
    clear addF;
    addF = UCDist >= 100;
    UCDist(addF) = 99; 
    
    % Find PoolDist above Max
    clear addF;
    addF = PoolDist >= 100;
    PoolDist(addF) = 99;
    
% --------Evaluate the inputs using the fuzzy inference system------------

disp('Doing FIS Calculations');

total = length(LWCount);
pct = round(total/10);
percent = 0;

CoverFIS = zeros(1, length(LWCount));

for i=1:(length(LWCount)) 
    
    %CoverFIS(1,i) = evalfis([UCDist(i) LWCount(i) PoolDist(i)], aFIS);
    CoverFIS(i) = evalfis([LWCount(i) UCDist(i) PoolDist(i)], aFIS);  
    
     if (i == pct)
         pct = pct + (round(total/10));
         percent = percent + 10;
         fprintf('%u percent done with habitat calculations.\n',percent);
     end
end

disp('Done evaluating your input file');       

%------Write the FIS Habitat Suitability Results to an CSV file format-----

fid3 = fopen(outfilename, 'wt');    %create output file to write to

fprintf(fid3, strcat('x,', 'y,', 'LWCount,', 'UCDist,', 'PoolDist,', 'FuzzyCoverIndex\n')); % write header

total = length(LWCount);
pct = round(total/10);
percent = 0;

% write out data
for j=1:(length(CoverFIS))                                                      

    fprintf(fid3,'%14.3f,%14.3f,%8.3f,%8.3f,%8.3f,%5.2f\n',  x(j), y(j), LWCount(j), UCDist(j), PoolDist(j), CoverFIS(j));
    if (j == pct)
         pct = pct + (round(total/10));
         percent = percent + 10;
         fprintf('%u percent done writing output file.\n',percent);
     end

end

fclose(fid3);

disp('Program finished.');  

end


