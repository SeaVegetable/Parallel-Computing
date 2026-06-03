# Matrix-Free Jacobian 预计算提示词

这个提示词根据当前需求自动整理，用来约束后续实现和迭代方向。

## 目标

在 `matrix_free` CPU 版本中，基于现有 `ElementMF` 的几何映射实现，预计算每个分区内：

1. 每个单元。
2. 每个积分点。
3. 对应的 `2x2` Jacobian 逆矩阵 `J^{-1}`。
4. 与该逆矩阵配套的 Jacobian 行列式 `det(J)`。

并将结果写入 `.h5` 文件。

## 输入约束

1. 继续复用现有 `info.txt`、`part_*.txt`、`NURBSExtraction1/2`、`IEN`、`CP` 数据。
2. 不改动当前 `matrix_free` 的主计算流程和文本分区格式。
3. 积分点配置默认沿用 matrix-free 当前实现，即 `x` 与 `y` 方向分别使用 `p+1`、`q+1` 个 Gauss 点。

## 输出约束

1. 每个分区输出一个 `.h5` 文件。
2. HDF5 至少包含：
   - 分区网格元信息。
   - 积分点与权重。
   - `inv_jacobian`，形状为 `[nElem, nQP, 2, 2]`。
   - `det_jacobian`，形状为 `[nElem, nQP]`。
3. 数据布局要固定、可复现，并在文档中明确说明单元与积分点的展平顺序。

## 实现约束

1. 尽量复用 `ElementMF` 现有 Jacobian 公式，不要重新发明一套几何映射逻辑。
2. 新增独立预处理入口，避免把 HDF5 写盘耦合进 matvec 主路径。
3. 构建系统需要显式接入 HDF5。

## 验收标准

1. 新程序可以从现有 partition 文本文件生成 `.h5`。
2. 构建通过。
3. HDF5 文件中能直接读取每个单元、每个积分点的 `J^{-1}` 与 `det(J)`。
