# MFSF

这个目录是从 `matrix_free/` 里的 `ElementMFSF` 思路抽出来的独立版本。

## 设计目标

- 保留 `matrix_free` 的核心思路：
  - 每个方向先给 1D basis 和 1D derivative
  - 通过张量积组装单元上的多维 basis
  - 通过控制点恢复几何映射与 Jacobian
- 用模板把维度泛化成 `Dim`
- 当前支持 `Dim = 2` 和 `Dim = 3`

## 核心类

- `MFSFElement<Dim>`
- `BSplineBasis`
- `BernsteinBasis`
- `RefElement`
- `QuadraturePoint`
- `NURBSExtractionGenerator`
- `MFSFIENGenerator<Dim>`
- `MFSFIDGenerator<Dim>`
- `MFSFControlPointGenerator`
- `MFSFPrecompute<Dim>`

它处理单个积分点上的：

- tensor-product basis
- 参数域导数
- 物理坐标
- Jacobian
- inverse Jacobian
- determinant
- 物理空间导数

## 新增生成器

这次把 `matrix_free/` 里的两类基础数据结构也带进来了，并做了泛化：

- `MFSFIENGenerator<Dim>`
  - 生成单元到局部基函数的 connectivity
  - 支持直接从 `num_funcs + degrees` 生成
  - 也支持从 `BSplineBasis` 数组生成

- `MFSFIDGenerator<Dim>`
  - 生成自由度编号
  - 边界自由度记为 `-1`
  - 内部自由度按连续编号 `0,1,2,...`

- `MFSFControlPointGenerator`
  - 从 `BSplineBasis` 和几何包围盒生成张量积控制点
  - 当前支持和 `matrix_free` 一样的 1D 求解，再推广到 2D/3D

- `MFSFPrecompute<Dim>`
  - 串起 basis / extraction / ID / IEN / control points / quadrature
  - 预计算每个单元每个积分点上的：
    - `inv_jacobian`
    - `det_jacobian`

- `MFSFH5Writer`
  - 把 `MFSFPrecomputeData<2/3>` 直接写成 `.h5`
  - 当前输出分组：
    - `/meta`
    - `/topology`
    - `/mesh`
    - `/quadrature`
    - `/geometry`

说明：

- 相比 `matrix_free` 旧版 2D `ID` 编号，这里统一成了“内部自由度连续编号”的风格，三维时更自然。
- `IEN` 采用张量积顺序，`x` 方向是最快变化方向，这和 `matrix_free` 当前二维实现是一致的。

## 输入数据约定

`EvaluateSingleQP(...)` 需要三组输入：

1. `basis_1d`
   - 类型：`std::array<std::vector<double>, Dim>`
   - 第 `a` 个方向长度应为 `degrees[a] + 1`

2. `dbasis_1d`
   - 类型：`std::array<std::vector<double>, Dim>`
   - 和 `basis_1d` 一一对应
   - 默认表示对参数坐标的导数

3. `control_points`
   - 展平方式：`[local_basis_id][coord]`
   - 长度必须是 `GetNumLocalBasis() * Dim`

## 构建

```bash
cmake -S . -B build
cmake --build build
./build/mfsf_demo
```

运行 demo 后会在当前工作目录生成：

- `mfsf_precompute_2d.h5`
- `mfsf_precompute_3d.h5`
