%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       FUZZY INFERENCE SYSTEM for Fish Habitat Suitability
%        for use with Fuzzy Logic Toolbox 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% This script uses a fuzzy inference system to estimate habitat
% suitability. The inputs are derived from Champ data (Delft model
% run output [velcity, depth] and grain size [D50 from ocular channel
% unit ocular estimates]). The script reqires inputs of Depth, Velocity  
% D50, and Cover in CSV input file.  
%
% User Inputs: A CSV file with a header row (ignored) and 1st column is
% easting (x), 2nd column is northing (y), 3rd column is velocity, 4th
% column is depth, 5th column is D50, 6th column is cover.

% Parameters:  None
% Outputs:     The program writes a new CSV file that is the same as the input,
% with an added column of fuzzy habitat suitability index (between 0 and 1).
%

function FuzzyHSI(filepath, FIStype, numInputs)

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
outpathname = strrep(pathstr, 'Inputs', 'Output');
outfilename = strcat(outpathname, '\', outname, '.csv');
outfilename = outfilename{1};
disp(outfilename);
	
%% Use csvread to get rest of the data
data = csvread(filepath,2,0);
x = data(:,1);
y = data(:,2);
Vel = data(:,3);
Depth = data(:,4);
D50 = data(:,5);
Cover = data(:,6);

%% Clear temporary variables
clearvars data raw R columnIndices;

disp('Done reading input file.');
    
%%------Get Rid of NonZero-----------------------------------------------    

% Find addresses of rows that meet condition
% Update array 
% Clear addresses after each find, update
    
    % Find Vels >= 0
    clear addF;
    addF = find(Vel >= 0);
    x = x(addF);
    y = y(addF);  
    Vel = Vel(addF);
	Depth = Depth(addF);
	D50 = D50(addF);
    Cover = Cover(addF);
    
    % Find Depths >= 0
    clear addF;
    addF = find(Depth >= 0);
    x = x(addF);
    y = y(addF);  
    Vel = Vel(addF);
	Depth = Depth(addF);
	D50 = D50(addF);
    Cover = Cover(addF);
    
    % Find D50 >= 0
    clear addF;
    addF = find(D50 >= 0);
    x = x(addF);
    y = y(addF);  
    Vel = Vel(addF);
	Depth = Depth(addF);
	D50 = D50(addF);
    Cover = Cover(addF);
    
    % Find Cover >= 0
    clear addF;
    addF = find(Cover >= 0);
    x = x(addF);
    y = y(addF);  
    Vel = Vel(addF);
	Depth = Depth(addF);
	D50 = D50(addF);
    Cover = Cover(addF);
    
    % Find input values > fis max
    % Set value to just below fis max so that an error is 
    
    % Find Vel above Max
    clear addF;
    addF = Vel >= 4;
    Vel(addF) = 3.9;
    
    % Find Depths above Max
    clear addF;
    addF = Depth >= 4;
    Depth(addF) = 3.9; 
    
    % Find D50 above Max
    clear addF;
    addF = D50 >= 4000;
    D50(addF) = 3999;
    
    % Find D50 above Max
    clear addF;
    addF = Cover >= 1;
    Cover(addF) = 0.99;
    
% --------Evaluate the inputs using the fuzzy inference system------------

disp('Doing FIS Calculations');

total = length(Vel);
pct = round(total/10);
percent = 0;

FuzzyHSI = zeros(1, length(Vel));

if numInputs == 3

	for i=1:(length(Vel)) 
		
		FuzzyHSI(i) = evalfis([Depth(i) Vel(i) D50(i)], aFIS);  
		
		 if (i == pct)
			 pct = pct + (round(total/10));
			 percent = percent + 10;
			 fprintf('%u percent done with habitat calculations.\n',percent);
		 end
	end

	disp('Done evaluating your input file');       

	%------Write the FIS Habitat Suitability Results to an CSV file format-----

	fid3 = fopen(outfilename, 'w');    %create output file to write to

	fprintf(fid3, strcat('x,', 'y,', 'Vel,', 'Depth,', 'D50,', 'FuzzyHSI\n')); % write header

	total = length(Vel);
	pct = round(total/10);
	percent = 0;

	% write out data
	for j=1:(length(FuzzyHSI))                                                      

		fprintf(fid3,'%14.3f,%14.3f,%8.3f,%8.3f,%8.3f,%5.2f\n',  x(j), y(j), Vel(j), Depth(j), D50(j), FuzzyHSI(j));
		if (j == pct)
			 pct = pct + (round(total/10));
			 percent = percent + 10;
			 fprintf('%u percent done writing output file.\n',percent);
		 end

	end

	fclose(fid3);
	
else

	for i=1:(length(Vel)) 
		
		FuzzyHSI(i) = evalfis([Depth(i) Vel(i) D50(i) Cover(i)], aFIS);  
		
		 if (i == pct)
			 pct = pct + (round(total/10));
			 percent = percent + 10;
			 fprintf('%u percent done with habitat calculations.\n',percent);
		 end
	end

	disp('Done evaluating your input file');       

	%------Write the FIS Habitat Suitability Results to an CSV file format-----

	fid3 = fopen(outfilename, 'w');    %create output file to write to

	fprintf(fid3, strcat('x,', 'y,', 'Vel,', 'Depth,', 'D50,', 'Cover,', 'FuzzyHSI\n')); % write header

	total = length(Vel);
	pct = round(total/10);
	percent = 0;

	% write out data
	for j=1:(length(FuzzyHSI))                                                      

		fprintf(fid3,'%14.3f,%14.3f,%8.3f,%8.3f,%8.3f,%8.3f,%5.2f\n',  x(j), y(j), Vel(j), Depth(j), D50(j), Cover(j), FuzzyHSI(j));
		if (j == pct)
			 pct = pct + (round(total/10));
			 percent = percent + 10;
			 fprintf('%u percent done writing output file.\n',percent);
		 end

	end

	fclose(fid3);
	
disp('Program finished.');  

end


