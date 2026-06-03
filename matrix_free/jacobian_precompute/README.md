# Jacobian Precompute

这个目录收纳了本次 Jacobian 预计算需求的提示词和 HDF5 写盘实现。

## HDF5 布局

生成的每个分区文件包含以下数据集：

- `/meta/p`
- `/meta/q`
- `/meta/nlocalelemx`
- `/meta/nlocalelemy`
- `/meta/nqp1`
- `/meta/nqp2`
- `/quadrature/points_x`
- `/quadrature/points_y`
- `/quadrature/weights_x`
- `/quadrature/weights_y`
- `/mesh/elem_size1`
- `/mesh/elem_size2`
- `/geometry/inv_jacobian`
- `/geometry/det_jacobian`

其中：

- `inv_jacobian` 的形状是 `[nElem, nQP, 2, 2]`，按行主序存储。
- `det_jacobian` 的形状是 `[nElem, nQP]`。
- `nElem = nlocalelemx * nlocalelemy`。
- 单元顺序与 `matrix_free` 现有代码一致：`elem = jy * nlocalelemx + ix`。
- 积分点顺序与 `ElementMF::GenerateElement(...)` 一致：`qp = qy * nqp1 + qx`。

## 运行方式

构建后可执行：

```bash
./precompute_jacobian [info.txt] [partition_dir] [output_dir]
```

默认行为：

- `info.txt` 默认为当前目录下的 `info.txt`
- `partition_dir` 默认为 `info.txt` 所在目录
- `output_dir` 默认为 `partition_dir/jacobian_precompute_h5`
