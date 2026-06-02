%% ========= 用户输入 =========
filename = 'part_0.txt';
p = 4;
numElem = 16;                 % Number of elem
tag = 'NURBSExtraction1';
tol = 1e-12;

%% ========= 参数 =========
n = p + 1;                     % 每个提取算子的矩阵尺寸
needNum = n * n;               % 每个提取算子的数字个数
totalNeed = needNum * numElem; % 总数字个数

%% ========= 打开文件 =========
fid = fopen(filename, 'r');
if fid == -1
    error('无法打开文件: %s', filename);
end

%% ========= 查找标签并读取全部提取系数 =========
found = false;
vals = [];

while ~feof(fid)
    line = strtrim(fgetl(fid));

    if strcmp(line, tag)
        found = true;

        % 从下一行开始持续读取，直到凑够 totalNeed 个数
        while ~feof(fid) && numel(vals) < totalNeed
            nextLine = strtrim(fgetl(fid));

            if isempty(nextLine)
                continue;
            end

            nums = sscanf(nextLine, '%f').';

            % 如果这一行读不到数字，就跳过或停止
            if isempty(nums)
                continue;
            end

            vals = [vals, nums]; %#ok<AGROW>
        end

        break;
    end
end

fclose(fid);

%% ========= 检查 =========
if ~found
    error('没有找到标签 %s', tag);
end

fprintf('读到的 extraction 系数总数: %d\n', numel(vals));
fprintf('理论应为: %d\n', totalNeed);

if numel(vals) < totalNeed
    error('提取系数数量不足：需要 %d 个，只读到 %d 个。', totalNeed, numel(vals));
elseif numel(vals) > totalNeed
    warning('提取系数数量多于预期：需要 %d 个，实际读到 %d 个。只取前 %d 个。', ...
        totalNeed, numel(vals), totalNeed);
    vals = vals(1:totalNeed);
end

%% ========= 拆分成每个单元的提取算子矩阵 =========
C_all = cell(numElem, 1);

for e = 1:numElem
    idx1 = (e-1)*needNum + 1;
    idx2 = e*needNum;

    block = vals(idx1:idx2);

    % 假设文件中按“行展开”存储
    C_all{e} = reshape(block, n, n).';
end

fprintf('成功恢复出 %d 个提取算子矩阵，每个大小为 %d x %d。\n', numElem, n, n);

%% ========= 显示前几个提取算子 =========
showNum = min(numElem, 3);
for e = 1:showNum
    fprintf('\n第 %d 个提取算子矩阵:\n', e);
    disp(C_all{e});
end

%% ========= 检查相邻提取算子是否相同 =========
fprintf('\n开始检查相邻提取算子是否相同...\n');

same_adjacent = false(numElem-1, 1);

for e = 1:numElem-1
    maxDiff = max(abs(C_all{e} - C_all{e+1}), [], 'all');

    if maxDiff < tol
        same_adjacent(e) = true;
        fprintf('第 %d 个与第 %d 个提取算子相同，maxDiff = %.3e\n', ...
            e, e+1, maxDiff);
    else
        fprintf('第 %d 个与第 %d 个提取算子不同，maxDiff = %.3e\n', ...
            e, e+1, maxDiff);
    end
end

%% ========= 汇总 =========
idx_same = find(same_adjacent);

fprintf('\n========== 汇总 ==========\n');
fprintf('相同的相邻提取算子对数: %d\n', numel(idx_same));

if isempty(idx_same)
    fprintf('没有相同的相邻提取算子。\n');
else
    fprintf('相同的相邻对如下：\n');
    for k = 1:numel(idx_same)
        fprintf('(%d, %d)\n', idx_same(k), idx_same(k)+1);
    end
end