function [] = Change_file(input_file,output_file,va2search,delimiter,num2write)

% open input file
[fidIn,errmsg] = fopen(input_file); % ,'r','n','Macintosh')
if fidIn < 3
   disp(errmsg);
end

tline = 'dum'; index = 1; flag = 0; kk = [];

while ischar(tline)
    tline = fgets(fidIn);
    kk = strfind(tline,va2search);
    if ~isempty(kk) && ~flag,
        flag = 1;
        tlinef = tline;
        indexf = index;
    end
    data{index} = tline;
    index = index + 1;
end

if ~exist('tlinef','var')
   error(['Change_file.m: the string ',va2search,' was not found in the file: ',input_file]) 
end

tlinef_cut = strrep(tlinef,va2search,'');

split_str2 = strsplit(tlinef_cut,{delimiter,';'});
split_str2{2}= [' ',delimiter,' ',num2str(num2write,16),';'];

data{1,indexf} = strjoin([{va2search} split_str2]);

fclose(fidIn);

% write to new file
% open input file
fidOut = fopen(output_file,'w');

for l = data
    fprintf(fidOut,'%s', l{1,1});
end
fclose(fidOut);