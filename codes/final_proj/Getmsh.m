function [ID, IEN, numNodes] = Getmsh(filename)
    % 输出:
    %   ID       - 节点坐标矩阵 [节点编号, x 坐标, y 坐标]
    %   IEN      - 单元连接矩阵 [单元编号, 节点编号1, 节点编号2, 节点编号3, 节点编号4]
    %   numNodes - 节点总数

    % 初始化输出变量
    ID = [];
    IEN = [];
    numNodes = 0;

    try
        % 打开文件
        fid = fopen(filename, 'r');
        if fid == -1
            error('无法打开文件: %s', filename);
        end

        % 定义标志位，用于跟踪当前读取的部分
        isNodeSection = false;
        isElementSection = false;

        % 按行读取文件内容
        while ~feof(fid)
            line = strtrim(fgetl(fid)); % 读取当前行并去掉多余的空格

            % 如果当前行是 $Nodes，开始读取节点部分
            if contains(line, '$Nodes')
                isNodeSection = true;
                numNodes = str2double(fgetl(fid)); % 从下一行读取节点总数
                continue;
            end

            % 如果当前行是 $EndNodes，结束读取节点部分
            if contains(line, '$EndNodes')
                isNodeSection = false;
                continue;
            end

            % 如果当前行是 $Elements，开始读取单元部分
            if contains(line, '$Elements')
                isElementSection = true;
                fgetl(fid); % 跳过下一行（单元部分的头信息）
                continue;
            end

            % 如果当前行是 $EndElements，结束读取单元部分
            if contains(line, '$EndElements')
                isElementSection = false;
                continue;
            end

            % 处理节点部分的数据
            if isNodeSection
                % 按格式解析节点数据 [节点编号, x, y, z]
                nodeData = sscanf(line, '%d %f %f %f');
                if numel(nodeData) == 4
                    % 存储节点编号和 x, y 坐标 (忽略 z 坐标)
                    ID = [ID; nodeData(1), nodeData(2:3)];
                end
            end

            % 处理单元部分的数据
            if isElementSection
                % 按格式解析单元数据 [单元编号, 类型, 标签, 节点编号...]
                elementData = sscanf(line, '%d %d %d %d %d %d %d');
                if numel(elementData) >= 5
                    % 获取单元类型
                    elementType = elementData(2);
                    if elementType == 1
                        % 如果是线单元，可以选择跳过（本例中仅处理四边形单元）
                        continue;
                    elseif elementType == 4
                        % 如果是四边形单元，提取单元编号和节点编号
                        elementID = elementData(1); % 单元编号
                        nodeIDs = elementData(6:end); % 关联的节点编号
                        % 存储到 IEN 中
                        IEN = [IEN; elementID, nodeIDs'];
                    end
                end
            end
        end

        % 关闭文件
        fclose(fid);

    catch ME
        % 如果发生错误，关闭所有已打开的文件
        fclose('all');
        rethrow(ME); % 抛出错误信息
    end
end
