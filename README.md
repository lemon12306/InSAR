# InSAR

车载毫米波 SAR / InSAR MATLAB 仿真代码，用于比较全精度与零阈值单比特量化的成像、干涉相位和高程反演结果。

## 环境

- 本项目在 MATLAB R2024b 下运行验证。
- 需要提供 `hann` 等函数的 Signal Processing Toolbox。
- 在 MATLAB 中将当前文件夹切换到本项目目录后运行脚本。

## 主要脚本

| 文件 | 用途 |
| --- | --- |
| `SAR_Parameters.m` | 雷达、观测几何、高斯山丘及散射点参数 |
| `RawDataGen_rand.m` | 双通道原始 LFM 回波或理想距离压缩回波生成 |
| `test.m` | 多组 SNR 及无噪声条件下的全精度与单比特对比实验 |
| `test_random_amplitude_3d.m` | 固定散射点几何，使用三个随机幅度种子运行无噪声实验，保存三维高程与剖面结果 |
| `inspect_onebit_results.m` | 读取已保存结果，进行正弦拟合检查，打开并保留可旋转的三维高程和剖面图窗 |
| `RDA_Imaging.m` | RDA 成像实验 |
| `BP_Imaging.m` | BP 成像实验脚本 |
| `PhaseFilter_Unwrap.m` | 读取 RDA 结果进行相位处理的实验脚本 |

## 运行

完整 SNR 对比实验：

```matlab
test
```

随机幅度实验和交互式三维查看：

```matlab
test_random_amplitude_3d
inspect_onebit_results
```

`inspect_onebit_results` 需要先生成 `随机幅度_无噪声_高程结果.mat`。查看脚本从该文件读取已反演的高程矩阵，黑线标出中心距离向剖面；脚本结束后不退出 MATLAB，也不关闭图窗。

当前 `SAR_Parameters.m` 的山丘增量高度为 5 m，基底高度为 0.1 m；散射幅度采用 `[0.5, 1.5)` 区间内的随机正实数，不引入随机散射相位。可直接修改参数后重新运行实验。

## 结果与数据

运行生成的 `.mat`、`.fig`、图像、统计 CSV 和日志保留在本地，不作为源码上传。`PhaseFilter_Unwrap.m` 需要的 RDA 数据可由相应成像实验生成。

当前三组随机幅度实验可观察到单比特高程起伏和尖峰，部分局部段接近正弦，但这些结果尚不能证明存在固定周期或固定相位的正弦规律。
