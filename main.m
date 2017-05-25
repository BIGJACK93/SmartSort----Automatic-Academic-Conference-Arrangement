%   Copyright 2017 Wei Cao, Parastoo Sadeghi

clc;
clear;
close all;


%% load xlsx file containing all paper information

init_data = readtable('isit.xlsx');
init_data.Properties.VariableNames{1} = 'Sequence';
journal_num = length(init_data.Sequence);


%% load isit arrangement file

manual_assign = readtable('manual.csv');
manual_assign.Properties.VariableNames{1} = 'Title';
manual_assign.Properties.VariableNames{2} = 'Position';
manual_assign.Properties.VariableNames{3} = 'ID';


%% getting dictionary file ready

termostat_data = load('database.mat');
keywords_list = struct2cell(termostat_data);
keywords_list = keywords_list{1};

keyword_bank = strtrim(keywords_list);
keyword_cell = upper(keyword_bank);
keyword_num = length(keyword_cell);
removal_list = []; 

for index = 1: keyword_num
    pick = keyword_cell{index};
    for  index2 = 1: keyword_num
        pick2 = keyword_cell{index2};
        state = strfind(pick2, pick);
        if isempty(state) == 0 && index ~= index2
            removal_list(end+1) = index;
        end
    end
end

removal_list = unique(removal_list, 'sorted')';
keyword_cell(removal_list) = [];


%% creat registry for information

keyword_library = {};
topic_library = {};
kw_registry = cell(journal_num, 6);
ppr_registry = cell(journal_num, 2);
tp_registry = cell(journal_num, 2);
tpc_registry = cell(journal_num, 2);


%% slice off useful information

seq_data = init_data.Sequence;
tit_data = init_data.Title;
tpc_data = init_data.PaperTopics;
abs_data = init_data.Abstract;
kw_data = init_data.Keywords;


%% topic registration

for index = 1:journal_num
    tp_registry{index, 1} = seq_data(index);
    tp_registry{index, 2} = tpc_data{index};
    topic_library{end+1} = tpc_data{index};
end

topic_library = unique(topic_library, 'sorted');


%% determine manual assignment approach from raw file

ref_titnsess = manual_assign.Title;
ref_pos = manual_assign.Position;
ref_id = manual_assign.ID;

ref_sess_raw = ref_titnsess(find(ref_id));
ref_tit_raw = ref_titnsess(find(~ref_id));
% gvn_tit = upper(tit_data);
gvn_tit = tit_data;
% ref_tit = upper(ref_tit);
ref_tit = ref_tit_raw;


ref_leng = length(ref_tit);
ref_sess_num = length(ref_sess_raw);


for index = 1:ref_sess_num
    pick = ref_sess_raw{index};
    ref_sess_raw{index} = strtrim(regexprep(pick, '\d+(?:_(?=\d))?', ''));
end


%% mark sessions on the same topic

dup_marker = zeros(ref_sess_num);


for index1 = 1: ref_sess_num
    pick1 = ref_sess_raw{index1};
    for index2 = 1: ref_sess_num
        pick2 = ref_sess_raw{index2};
        state = strfind(pick1, pick2);
        if isempty(state) == 0
            dup_marker(index1, index2) = 1;
            dup_marker(index2, index1) = 1;
        end
    end
end


overlap = ismember(ref_tit, gvn_tit);
overlap_pos = find(overlap == 1);
overlap_num = length(overlap_pos);

if overlap_num ~= journal_num
    error('Critial error. Please check your input files.')
end


%% generate numerical form of the title list from the manual assignment file

num_dummy = zeros(1, ref_leng);

for index = 1:ref_leng
    
    pick = ref_tit{index};
    bool_bar = ismember(gvn_tit, pick);
    tit_pos = find(bool_bar);
    stat = isempty(tit_pos);
    
    if stat == 1
        tit_pos = 0;
    end
    
    num_dummy(index) = tit_pos;
    
end
 

%% group the papers into sessions using the boolean tracking array

max_sess_cap = 6; %according to file
ref_cupboard = zeros(ref_sess_num, max_sess_cap);
tot_leng = length(ref_titnsess);
cursor_y = 0;
cursor_x = 0;
mapper = 1;
counter = 0;

for index = 1: tot_leng
    
    id = ref_id(index);
    
    if id == 1
        cursor_x = 0;
        cursor_y = cursor_y + 1;
        counter = counter + 1;
        index = index + 1;
        
    else
        mapper = index - counter;
        pick = num_dummy(mapper);
        cursor_x = cursor_x + 1;
        ref_cupboard(cursor_y, cursor_x) = pick;
        
        if mapper > ref_leng
            error('Critial error. Please check your code.')
        end
        
    end
end


%% keyword extraction of each sections

for index = 1:journal_num
    
    kw_registry{index, 1} = seq_data(index);
    
    temp_data_1 = tit_data{index};
    keyword_set_1 = extract_keyword(temp_data_1, keyword_cell);
    kw_registry(index, 2) = {keyword_set_1};
    for index2 = 1:length(keyword_set_1)
        keyword_library{end+1} = keyword_set_1{index2};
    end
    
    temp_data_2 = tpc_data{index};
    keyword_set_2 = extract_keyword(temp_data_2, keyword_cell);
    kw_registry(index, 3) = {keyword_set_2};
    for index2 = 1:length(keyword_set_2)
        keyword_library{end+1} = keyword_set_2{index2};
    end
    
%     temp_data_3 = abs_data{index};
%     keyword_set_3 = extract_keyword(temp_data_3, keyword_cell);
%     kw_registry(index, 4) = {keyword_set_3};
%     for i = 1:length(keyword_set_3)
%         keyword_library{end+1} = keyword_set_3{i};
%     end
    
    temp_data_4 = kw_data{index};
    keyword_set_4 = extract_keyword(temp_data_4, keyword_cell);
    kw_registry(index, 5) = {keyword_set_4};
    for index2 = 1:length(keyword_set_4)
        keyword_library{end+1} = keyword_set_4{index2};
    end
    
%     keyword_set_5 = {keyword_set_1{:}, keyword_set_2{:}, keyword_set_3{:}, keyword_set_4{:}};
    keyword_set_5 = {keyword_set_1{:}, keyword_set_2{:}, keyword_set_4{:}};
    keyword_set_5 = unique(keyword_set_5, 'sorted');
    kw_registry(index, 6) = {keyword_set_5};
    
end


%% final library of all keywords

keyword_sorted = unique(keyword_library, 'sorted');


%% registering keywords for each papaer in a cell

for index = 1:journal_num
    
  requested_words = kw_registry{index, 6};
  bool_bar = ismember(keyword_sorted, requested_words);
  indexes = find(bool_bar);
  ppr_registry{index, 1} = kw_registry{index, 1};
  ppr_registry{index, 2} = indexes;
  
end


%% registering topic for each papaer in a cell

for index = 1:journal_num
    
  requested_tpc = tp_registry{index, 2};
  bool_bar = ismember(topic_library, requested_tpc);
  indexes = find(bool_bar);
  tpc_registry{index, 1} = tp_registry{index, 1};
  tpc_registry{index, 2} = indexes;
  
end

%% start labelling keywords

kw_num = length(keyword_sorted);
sum_matrix = zeros(kw_num, journal_num);

for index = 1: journal_num
    tags = ppr_registry{index, 2};
    for index2 = 1: length(tags)
        dim_kw = tags(index2);
        sum_matrix(dim_kw, index) = 1;
    end
end


%% degree of correlations of papers are to be calculated

dim = size(sum_matrix);
mat_width = dim(2);
mat_height = dim(1);
combi_bank = combnk(1:mat_width, 2);
corr_mat = zeros(mat_width);

for index = 1: length(combi_bank)
    roll_calls = combi_bank(index, :);
    roll_call1 = roll_calls(1);
    roll_call2 = roll_calls(2);
    slice1 = sum_matrix(:, roll_call1);
    slice2 = sum_matrix(:, roll_call2);
    feeder = [slice1, slice2];
    corr_val = corr_stat(feeder);
    corr_mat(roll_call1, roll_call2) = corr_val;
    corr_mat(roll_call2, roll_call1) = corr_val;
end


%% replace all non-zero values by 1 to generate a boolean adjancy matrix

weighed_mat = corr_mat;
corr_mat(find(corr_mat>0)) = 1;


%% generate a cell containing all the matrices for greedy clique sorting

cliq_cell = cell(journal_num);
combi_bank2 = combnk(1:journal_num, 2);
comb_num = length(combi_bank2);

for index = 1:comb_num
    comb = combi_bank2(index,:);
    rc1 = comb(1);
    rc2 = comb(2);
    ext1 = ppr_registry{rc1, 2};
    ext2 = ppr_registry{rc2, 2};
    dim1 = length(ext1);
    dim2 = length(ext2);
    temp_mat = zeros(dim1, dim2);
    for index2 = 1: dim1
        for index3 = 1: dim2
            if ext1(index2) == ext2(index3)
                temp_mat(index2, index3) = 1;
            end
        end
    end
    cliq_cell{rc1, rc2} = temp_mat;
    cliq_cell{rc2, rc1} = temp_mat';
    
end

for index = 1:journal_num
    indi_kw = length(ppr_registry{index, 2});
    cliq_cell{index, index} = zeros(indi_kw);
end


%% transform cell into matrix for greedy clique selection

cliq_mat = cell2mat(cliq_cell);


%% greedy algorithm for maximal clique sorting

sess_cap = 4;
sess_num = floor(journal_num./sess_cap);
leftover = mod(journal_num, sess_cap);
employed = journal_num - leftover;
container = zeros(1, sess_cap);
cupboard = zeros(sess_num, sess_cap);
blacklist = [];
whitelist = [1:journal_num];
fulllist = [1:journal_num];
corr_column = zeros(journal_num, 1);


%% calculate degree of correlation for each individual paper

for index = 1: journal_num
    corr_slice = cliq_cell(index, :);
    merge_mat = cell2mat(corr_slice);
    corr_sum = sum(sum(merge_mat));
    corr_column(index) = corr_sum;
end


%% sort papers in terms of their degree of correlation

[corr_val,  index_column] = sort(corr_column, 'descend');
corr_registry = [index_column, corr_val];
pool = corr_registry;
photocopy = pool;

for index = 1:sess_num
    
    pool_spcs = pool(:, 2);
    pool_index = pool(:, 1);
    pick = max(pool_spcs);
    pick = pick(1);
    pos = find(pool_spcs == pick(1));
    pos = pos(1);
    index_val = pool_index(pos);
    
    bool_filter = ismember(pool_index, index_val);
    removal_pos = find(bool_filter);
    pool(removal_pos, :) = [];
    
    pool_spcs = pool(:, 2);
    pool_index = pool(:, 1);
    
    container(1) = index_val;
    
    selection_slice = weighed_mat(index_val, :);
    marked_slice = [selection_slice,
        1:journal_num];
    ava_mkt = marked_slice(:, pool_index);
    ava_sel = ava_mkt(1, :);
    mk_bar = ava_mkt(2, :);
    [ava_val, ava_pos] = sort(ava_sel, 'descend');
    ext_pos = ava_pos(1:sess_cap-1);
    cdd_num = mk_bar(ext_pos);
    container(2:sess_cap) = cdd_num;
    
    blacklist = [blacklist, container];
    leng1 = length(blacklist);
    blacklist = unique(blacklist);
    leng2 = length(blacklist);
    
    %error alert when critical error is detected
    
    if leng1 ~= leng2
        error('Critial error. Please check your code.')
    end
    
    whitelist = setxor(blacklist, fulllist);
    
    bool_filter = ismember(pool_index, container);
    removal_pos = find(bool_filter);
    pool(removal_pos, :) = [];
    
    cupboard(index, :) = container;
    container = zeros(1, sess_cap);
    
end
    
dims = size(cupboard);
leng = dims(1);
wid = dims(2)+1;
disp_cell = cell(leng, wid);

for index1 = 1:leng
    for index2 = 2:wid
        coord_val = cupboard(index1, index2-1);
        disp_cell{index1, index2} = tit_data{coord_val};
    end
end


%% asssigning title to each session

container2 = [];
title_bar = zeros(leng, 1);

for index1 = 1:leng
    slice = cupboard(index1, :);
    for index2 = 1:sess_cap
        ppr_index = cupboard(index1, index2);
        container2 = [container2, tpc_registry{ppr_index, 2}];
    end
    tpc_index = mode(container2);
    title_bar(index1) = tpc_index;
    tpc = topic_library{tpc_index};
    disp_cell{index1, 1} = tpc;
    container2 = [];
end
        

disp(disp_cell);


%% mark sessions on the same topic

dup_marker_2 = zeros(leng);

for index1 = 1:leng
    for index2 = 1:leng
        pick1 = title_bar(index1);
        pick2 = title_bar(index2);
        stat = (pick1==pick2);
        dup_marker_2(index1, index2) = stat;
        dup_marker_2(index2, index1) = stat; 
    end
end


%% evaluating the similarity of the auto assignment with the manual one

actual_size = numel(cupboard);
sima_array = zeros(actual_size, 1);
remove_marker = [];
cupboard_tr = cupboard';
cupboard_tr = cupboard_tr(:);

for index = 1:actual_size

    volume = cupboard_tr(index);
    [cup_x, cup_y] = find(cupboard == volume);
    [ref_x, ref_y] = find(ref_cupboard == volume);
    asso_slice_1 = dup_marker_2(cup_x, :);
    asso_slice_2 = dup_marker(ref_x, :);
    asso_list_1 = find(asso_slice_1);
    dims_asso = size(asso_slice_2);
    width_asso = min(dims_asso);
    if width_asso == 1
        asso_list_2 = find(asso_slice_2);
    else
        [useless_dim, asso_list_raw] = find(asso_slice_2);
        asso_list_2 = unique(asso_list_raw);
    end
    
    if isempty(cup_x)==0 && isempty(ref_x)==0
        cup_slice_0 = cupboard(asso_list_1, :);
        ref_slice_0 = ref_cupboard(asso_list_2, :);
        cup_slice = cup_slice_0(:);
        ref_slice = ref_slice_0(:);
        comm = intersect(cup_slice, ref_slice);
        num_cumm = length(comm);
        empty_slot = find(ref_slice == 0);
        ref_slice_leng = length(ref_slice) - length(empty_slot);
        divider = min(ref_slice_leng, length(cup_slice));
        sima_array(index) = num_cumm./divider;
    else
        remove_marker(end+1) = index;
    end
    
end

sima_array(remove_marker) = [];
simalarity = mean(sima_array).*100;
format_spec = 'The similarity between auto and manual assignment is: %-.2f%%.';
prt_str = sprintf(format_spec, simalarity);
disp(prt_str);


track_array = zeros(length(sima_array), 1);

for index = 1:length(sima_array)
    track_array(index) = mean(sima_array([1:index]));
end

indexes = [1:length(track_array)];

figure;
plot(indexes, track_array);
grid on;

figure;
scatter(indexes, sima_array, 'filled');
grid on;


%% test cell2uitable functions
clc
close all
[dim1, dim2] = size(disp_cell);
selection_bar = cell(dim1, dim2-1);
selection_bar(:) = {true};
disp_cell2 = cell(dim1, 2*dim2-1);
col_cord1 = [1:2:(2*dim2-3)]+1;
col_cord1 = [1, col_cord1];
disp_cell2(:, col_cord1) = disp_cell(:, [1:dim2]);
col_cord2 = [1:2:(2*dim2-3)]+2;
disp_cell2(:,col_cord2) = selection_bar(:, :);
tit_cell = cell(1, 2*dim2-1);
tit_cell{1} = ['Session Topic'];
tit_cord = col_cord1(2:end);
for index = 1:dim2-1
    pointer = tit_cord(index);
    content = ['Paper Number', num2str(index)];
    tit_cell{pointer} = content;
end

cell2uitable(disp_cell2, 'oversizecolfactor', 0.8, 'colnames', tit_cell)