% Procedure YMrCA - Platinum Pedigree Cluster Edition
% Оптимизировано для 80+ локусов без использования Excel

% --- 1. ПАРАМЕТРЫ СИМУЛЯЦИИ ---
G = 150;                % Максимальное число поколений
keep_max = 5000;        % Лимит комбинаций (C)
add_descrepancy = 15;   % Допустимое расхождение аллелей (AD)
min_prob_n0 = 1e-15;    % Минимальная вероятность
prob_decrease = 0.99;   % Коэффициент уменьшения вероятности


% --- 2. ЗАГРУЗКА ДАННЫХ ---
fprintf('Загрузка данных...\n');

% ИСПОЛЬЗУЕМ 'preserve', чтобы MATLAB не портил имена колонок с точками
haplo_table = readtable('haplo_matrix.csv', 'VariableNamingRule', 'preserve');
alleles_scale = load('alleles_scale.txt');

marker_names = haplo_table.Properties.VariableNames(2:end);
data = table2array(haplo_table(:, 2:end));
sample_names = string(haplo_table.name); % Преобразуем в массив строк

% Создаем все возможные пары (28 пар для 8 образцов)
pairs = combnk(1:size(data,1), 2);

% --- 3. ГЛАВНЫЙ ЦИКЛ ПО ПАРАМ ---
for koppel = 1:size(pairs, 1)
    idxA = pairs(koppel, 1);
    idxB = pairs(koppel, 2);
    
    fprintf('\n========================================\n');
    fprintf('Пара %d/%d: %s vs %s\n', koppel, size(pairs,1), sample_names(idxA), sample_names(idxB));
    
    % Выделяем валидные (не-NaN) локусы, общие для этой пары
    valA = data(idxA, :);
    valB = data(idxB, :);
    valid = ~isnan(valA) & ~isnan(valB);
    
    curr_markers = marker_names(valid);
    vA = valA(valid);
    vB = valB(valid);
    n_loci = numel(vA);
    
    fprintf('  Используется локусов: %d\n', n_loci);
    if n_loci < 5
        fprintf('  Слишком мало локусов, пропускаю...\n');
        continue;
    end

    % Находим индексы стартовых и целевых аллелей на шкале
    idx_start = zeros(1, n_loci);
    idx_compare = zeros(1, n_loci);
    for m = 1:n_loci
        [~, idx_start(m)] = min(abs(alleles_scale - vA(m)));
        [~, idx_compare(m)] = min(abs(alleles_scale - vB(m)));
    end

    % Базовое расхождение
    P_int_diff_marker = sum(abs(idx_start - idx_compare));
    if P_int_diff_marker < 3, P_int_diff_marker = 3; end

    % --- 4. БАЙЕСОВСКАЯ СИМУЛЯЦИЯ (АЛГОРИТМ АВТОРОВ) ---
    probgen = 1.0;
    idx_gen = idx_start;
    minprob = min_prob_n0 * prob_decrease;
    
    solution_prob = cell(G, 1);
    solution_idx = cell(G, 1);

    for n = 1:G
        fprintf('  Поколение: %d\r', n);
        
        % --- МАРКЕР 1 ---
        m1_file = ['matrices/', curr_markers{1}, '.csv'];
        M1 = readmatrix(m1_file);
        
        ind_t = M1(:, idx_gen(:,1)) >= minprob;
        [rown, coln] = find(ind_t);
        
        if isempty(rown), break; end
        
        subidx = sub2ind(size(M1), rown, idx_gen(coln, 1));
        prob = probgen(coln) .* M1(subidx);
        idx_gen_tmp = idx_gen(coln, :);
        
        possible_idx = zeros(numel(prob), n_loci, 'int8');
        possible_idx(:, 1) = int8(rown);
        
        % --- МАРКЕРЫ 2..N ---
        for i = 2:n_loci
            mi_file = ['matrices/', curr_markers{i}, '.csv'];
            Mi = readmatrix(mi_file);
            
            % Умножение вероятностей
            indtemp = (prob .* Mi(:, idx_gen_tmp(:, i))') >= minprob;
            [row_idx, col_idx] = find(indtemp);
            
            int_ind = numel(row_idx);
            if int_ind == 0, prob = []; break; end
            
            possible_idx_next = zeros(int_ind, i, 'int8');
            possible_idx_next(:, 1:i-1) = possible_idx(row_idx, 1:i-1);
            possible_idx_next(:, i) = int8(col_idx);
            
            subidx = sub2ind(size(Mi), col_idx, idx_gen_tmp(row_idx, i));
            prob = prob(row_idx) .* Mi(subidx);
            idx_gen_tmp = idx_gen_tmp(row_idx, :);
            possible_idx = possible_idx_next;
        end
        
        if isempty(prob), break; end

        % --- ГРУППИРОВКА И ОТСЕЧЕНИЕ ---
        [a, b, c] = unique(possible_idx, 'rows', 'first');
        [~, y] = sort(b);
        idx_new_gen = a(y, :);
        
        % Лексикографическая сортировка (аналог sortrows)
        [~, q] = sortrows(idx_new_gen);
        ic = q(c, :);
        
        % Адаптация порога
        if numel(prob) > (keep_max * 10)
            minprob = minprob / 0.7;
        elseif numel(prob) <= (keep_max * 5)
            minprob = minprob * prob_decrease;
        end
        
        prob_new_gen = zeros(size(idx_new_gen, 1), 1);
        for j = 1:numel(prob)
            prob_new_gen(ic(j)) = prob_new_gen(ic(j)) + prob(j);
        end
        
        % Фильтр AD
        dist_to_target = sum(abs(double(idx_new_gen) - double(idx_compare)), 2);
        keep = dist_to_target <= (P_int_diff_marker + add_descrepancy);
        prob_new_gen = prob_new_gen(keep);
        idx_new_gen = idx_new_gen(keep, :);
        
        % Лимит C
        if numel(prob_new_gen) > keep_max
            [~, I] = sort(prob_new_gen, 'descend');
            I = I(1:keep_max);
            probgen = prob_new_gen(I);
            idx_gen = double(idx_new_gen(I, :));
        else
            probgen = prob_new_gen;
            idx_gen = double(idx_new_gen);
        end
        
        solution_prob{n} = probgen;
        solution_idx{n} = idx_gen;
    end
    fprintf('\n');

    % --- 5. ФОРМИРОВАНИЕ РЕЗУЛЬТАТА ---
    probability = zeros(1, G);
    for n = 1:G
        if isempty(solution_idx{n}), continue; end
        % Ищем полное совпадение с целью
        [~, ~, ib] = intersect(idx_compare, solution_idx{n}, 'rows');
        if ~isempty(ib)
            probability(n) = solution_prob{n}(ib);
        end
    end
    
    % Сохраняем в текстовый CSV (одна строка с вероятностями)
    out_name = sprintf('couple%d.csv', koppel);
    writematrix(probability, out_name);
end

fprintf('\n*** РАСЧЕТ ВСЕЙ СЕМЬИ ЗАВЕРШЕН ***\n');
