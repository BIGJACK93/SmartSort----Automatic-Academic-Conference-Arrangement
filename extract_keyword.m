%extract keywords in the approach customised for the file%

function keywords = extract_keyword5(input_str, database_name)

mod_str = upper(input_str);
keyword_bank1 = strtrim(database_name);
keyword_cell1 = upper(keyword_bank1);
picked_word = {};
keyword_num1 = length(keyword_cell1);

for index = 1: keyword_num1 
    pick = keyword_cell1{index};
    state = strfind(mod_str, pick);
    if isempty(state) == 0
        picked_word{end+1} = pick;
    end
end

keywords = picked_word;

end