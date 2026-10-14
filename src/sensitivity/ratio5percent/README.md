# K/π 比例测量与纯 pion 排除灵敏度

本目录接收 Michel 模板，计算并绘制两张图：

- **比例测量**：固定事件数下扫描 χ²，以 `min(χ²) + 1` 标记原 notebook 的 1σ 区间。
- **纯 pion 排除**：计算不同 K/π 比例下，达到原分析 χ² = 9 阈值所需的事件数。

上游原始数据处理复用 `src/analysis_pmu_from_MusAirS/`，不复制筛选与模板生成代码。

## 目录职责

```text
ratio5percent/
├── analysis.py                     # 模板接口、混合、两种分析的数值计算
├── run_end_to_end.py                # 原始数据/标准模板 → 两种分析与两张图
├── run_pipeline.py                  # 历史模板 → notebook 的 1σ 图
├── plot_ratio.py                    # 单独重绘比例扫描结果
├── plot_exclusion_sensitivity.py    # 单独重绘纯 pion 排除结果
├── tests/                          # 单元、历史回归、合成 ROOT 集成测试
├── requirements.txt
├── README.md
└── results/                        # 本地生成结果，不提交 Git
```

两种入口有意保留：`run_end_to_end.py` 用于新分析，`run_pipeline.py` 用于历史图复现。
两个绘图脚本也可独立使用，不必为调整样式重新运行分析。

## 安装

在仓库根目录执行：

```bash
python3 -m venv src/sensitivity/ratio5percent/.venv
source src/sensitivity/ratio5percent/.venv/bin/activate
python -m pip install -r src/sensitivity/ratio5percent/requirements.txt
```

只有从原始数据开始时需要命令行 `root`（支持 C++17）。已有模板或只绘图时不需要 ROOT 可执行程序。

## 1. 完整分析入口

### 从原始 MusAirS 数据开始

```bash
python src/sensitivity/ratio5percent/run_end_to_end.py \
  --proton /path/to/proton.root \
  --helium /path/to/helium.root \
  --output-dir /path/to/new_results
```

流程为：原始质子/氦数据 → 极化分布 → μ−/μ+ Michel 模板 → 两种分析 → 两张图。
不依赖 `dev/` 中预存的历史文件，不修改原始输入。

上游默认：单位事例权重、质子/氦来源系数 `9.704844263299673`/`100.0`、
历史去重、50-bin 极化分布、100-bin Michel 分布，每个来源/电荷采样 100M 次，
随机种子 `20241025`。六个模板共需采样 600M 次。

**权重模式必须根据数据含义选择。** 若原始 `Weight` 已包含完整通量和样本归一化：

```bash
python src/sensitivity/ratio5percent/run_end_to_end.py \
  --proton /path/to/weighted_proton.root \
  --helium /path/to/weighted_helium.root \
  --raw-weights --proton-scale 1 --helium-scale 1 \
  --output-dir /path/to/new_weighted_results
```

`--no-legacy-deduplicate` 保留所有通过初筛的候选。原始字段、筛选与权重约定见
[上游 README](../../analysis_pmu_from_MusAirS/README.md)。

试跑可使用 `--template-events 1000000`，但模板噪声增大，结果不应直接用于正式结论。
模板采样数 `--template-events` 与实验事件数 `--events` 是不同参数。

### 从已有标准模板开始

```bash
python src/sensitivity/ratio5percent/run_end_to_end.py \
  --templates /path/to/michelData.root \
  --output-dir /path/to/new_analysis_results
```

常用参数：

| 参数 | 默认值/含义 |
| --- | --- |
| `--events` | 10M，比例测量的实验事件数 |
| `--true-ratio` | 0.05 |
| `--scan-start`、`--scan-end`、`--points` | 0.01、0.10、1000 |
| `--measurement-charge` | `muon`；也可选 `antimuon` |
| `--exclusion-charge` | `antimuon`；也可选 `muon` |
| `--ratio-axis` | `physical`，正确物理横坐标；可选 `notebook` 保留旧索引映射 |
| `--seed`、`--template-events` | 上游采样种子与每模板采样数 |

### 上下游接口

接口是 `michelData.root` 中的独立 `cosθ` TH1 模板，而不是图片或已计算好的灵敏度。

| 名称 | 含义 |
| --- | --- |
| `piMuon`、`KMuon`、`KLongMuon` | μ− 的三个母粒子来源 |
| `piAntiMuon`、`KAntiMuon`、`KLongAntiMuon` | μ+ 的三个母粒子来源 |
| `provenance` | 上游筛选、权重、种子、生成方法等记录 |

下游只使用所选电荷的 pion/kaon 模板，不混入 KLong。
要求 bin 内容有限、非负且积分为正；同一对模板的 bin 边界必须相同，范围为 `[-1,1]`。
下游分别归一化各来源，不从模板积分推断实验事件数。缺少所选电荷时明确报错。

上游 `generate_michel.C` 默认仍只输出三个 μ− 模板；新增最后一个参数
`includeAntiMuon=true` 时额外生成 μ+，保留原 μ− 随机流：

```bash
root -l -b -q 'src/analysis_pmu_from_MusAirS/generate_michel.C("pmu_output/polarizationData.root","templates_output",100000000LL,20241025,false,true)'
```

### 输出和安全

```text
new_results/
├── manifest.json
├── analyze_pmu.log / generate_michel.log    # 原始数据模式
├── upstream/                               # 原始数据模式：ROOT、CSV、溯源
├── measurement/
│   ├── michelDis.root
│   ├── measure_ratio_chi2test.root
│   ├── measure_ratio.pdf / .png
│   └── summary.json
└── exclusion/
    ├── eigen_sensitiviy_output.root
    └── 241013kpiRatioNSensitivity.pdf / .png
```

默认输出为本目录 `results/end_to_end/`。**完整入口拒绝覆盖非空输出目录。**
失败状态和原因写入 manifest，ROOT 诊断见日志；失败后重跑也需使用新目录。
manifest 保存输入路径、参数、ROOT 调用、可用的模板溯源及运行状态，不是内容校验和。
正式归档可额外保存原始文件 SHA-256 和 Git 提交号。

## 2. 历史 notebook 复现入口

```bash
python src/sensitivity/ratio5percent/run_pipeline.py
```

默认读取 `dev/megaMichelDis.root` 的 `piMichelDis`、`kaMichelDis`（μ−），
复现 `dev/2601reviseArticle/ratio5percent.ipynb` 第二张图。
这些模板与历史 `michelDis.root` 的参考形状一致；不是另一组
`2601_michelElectronWithWeight.root` 模板。
ROOT 数据不提交 Git，全新克隆需提供模板：

```bash
python src/sensitivity/ratio5percent/run_pipeline.py \
  --templates /path/to/templates.root \
  --pion-name piMichelDis --kaon-name kaMichelDis \
  --output-dir /path/to/results
```

此入口默认 `--ratio-axis notebook`，保留旧结果的横坐标偏移：

```text
physical(i) = 0.01 + i       * (0.10 - 0.01) / 999
notebook(i) = 0.01 + (i - 1) * (0.10 - 0.01) / 1000
```

旧图最低点横坐标为 `0.04987`，实际比例为 `0.05`。
正式分析推荐 `--ratio-axis physical`；树中始终保存 `physical_ratio`。
默认复现的两个蓝点横坐标为 `0.045415`、`0.054415`。

此入口输出到 `results/`，包括 `michelDis.root`、`measure_ratio_chi2test.root`、
`measure_ratio.pdf/.png` 与 `summary.json`，重复运行会覆盖同名文件。
旧 `measure_ratio_mega.root` 是第一张图的产物，不再被当前流程使用。

## 3. 只绘图

```bash
# 1σ 比例测量图
python src/sensitivity/ratio5percent/plot_ratio.py \
  /path/to/measure_ratio_chi2test.root \
  --output-dir /path/to/plots

# 纯 pion 排除事件数图
python src/sensitivity/ratio5percent/plot_exclusion_sensitivity.py \
  --input /path/to/eigen_sensitiviy_output.root \
  --output-dir /path/to/plots
```

排除图脚本不传 `--input` 时默认读取
`dev/2601reviseArticle/sensitivity_revise/eigen_sensitiviy_output.root`。
它读取 `sensitivityResult/ratio` 和 `targetStatistics`，不重新计算分析。
样式保留原虚线、对数轴及旋转文字，阴影为 `#f6f6ffff`，文字为
`#8a8af5ff`、加粗，均完全不透明。

绘图范围固定：比例测量 `x=[0.04,0.06]`、`y=[0,5]`；
纯 pion 排除 `x=[0.01,0.20]`、`y=[3e4,2e7]`。
排除图沿用原 `ratio (%)` 标签，但横坐标实际为分数 0.01–0.20（即 1%–20%）。
只绘图入口会覆盖同名图像，调整样式无须重新生成模板。

## 统计定义

- 混合分布为 `(p_pi + r*p_K)/(1+r)`，使用 Asimov 分布，不施加 toy 涨落。
- 等总计数且两组 bin 方差均为计数时，ROOT `WW CHI2` 化简为
  `χ² = Σ(n-μ)²/(n+μ)`；不是第一张历史图使用的单样本 Neyman χ²。
  两组都为零的 bin 不贡献 χ²，只有一组为零的 bin 仍参与计算。
- 比例测量以 `min(χ²)+1` 求交，图中保留 notebook 的中点近似，
  Brent 精确交点另存到 JSON。扫描范围没有两个交点时明确报错。
- 排除计算使用线性缩放 `N_target=9/χ²(N=1)`，比例列表沿用历史九点。
  两种模板完全不可区分时明确报错，不输出伪造的有限事件数。
- 保留原图的 1σ/3σ 阈值约定，不验证覆盖率，也不称为单样本 profile likelihood。
  新模板采用直接采样，不保证与历史加权蒙特卡洛模板逐点一致。

## 测试

```bash
python -m unittest discover -s src/sensitivity/ratio5percent -v
python -m unittest discover -s src/analysis_pmu_from_MusAirS -p test_pipeline.py -v
```

测试文件集中在 `tests/`，保留计算、绘图、拒绝覆盖、原始数据端到端、
μ− 随机流兼容等验证。无 ROOT 时跳过原始数据集成测试；
无本地 `dev/` 参考数据时跳过历史回归测试。
