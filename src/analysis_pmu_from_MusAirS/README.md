# 从 MusAirS 原始数据生成 μ 子极化图和 Michel 角分布图

本目录保存可纳入 Git 的入口、处理规则、绘图脚本和测试，不依赖 `dev/` 的代码或中间数据。

```
质子 + 氦 MusAirS ROOT 文件
  → analyze_pmu.C：初筛、历史去重、母粒子/动能筛选、分量合并、归一化
  → polarizationData.root + polarizationData.csv + provenance.txt
  ├─ plot_pmu.py → combinedPmu_revise.pdf
  └─ generate_michel.C：三个独立来源、无权逆 CDF 直接采样角分布
       → michelData.root + michelData.csv + michel_provenance.txt
       → plot_michel.py → michelComponents.pdf
```

这里直接生成最终绘图所需的六个直方图，不生成 `muonData.root`；母粒子动能 `ParEk` 不参与这张图的筛选，因此不需要提取或保存。

## 依赖

- ROOT 6，支持 C++17 / `ROOT::RDataFrame`；本地测试版本为 6.38.02。
- Python 3.9+、Matplotlib 和 SciPy；Python 端只读取 CSV，不需要 PyROOT 或 uproot。已有环境请重新安装 `requirements.txt`，新增的 Michel 插值需要 SciPy。

以下命令均以 **GME 工程根目录**为工作目录：

```bash
python3 -m venv src/analysis_pmu_from_MusAirS/.venv
source src/analysis_pmu_from_MusAirS/.venv/bin/activate
python -m pip install -r src/analysis_pmu_from_MusAirS/requirements.txt
```

## 一条命令完成生成和绘图

在依赖已安装的 Python 环境中：

```bash
bash src/analysis_pmu_from_MusAirS/run_pipeline.sh \
  /Users/sunmingchen/codeDIR/2502_MusAirSValidation/CRmuData/MusAirS/MusAirS_AMS_polarization.root \
  /Users/sunmingchen/codeDIR/2502_MusAirSValidation/CRmuData/MusAirS/musAirSHelium.root \
  pmu_output
```

上面的两个路径是本机找到的同名文件位置，不是程序内置路径。已经完整试跑：新输出的六个直方图与 `dev/polarizationData20.root` 的所有 bin 内容和误差（含 underflow/overflow）完全一致，最大差值均为 0。这证明它们能复现旧绘图数据，但不等于证明它们与历史外接盘文件逐字节相同。换机器时传入自己的文件路径。脚本可从任意工作目录调用；相对输入/输出路径相对于调用者的工作目录。可用 `PYTHON=/path/to/python` 指定解释器。

输出目录默认 `pmu_output`。一键入口现在同时生成两张 PDF：`combinedPmu_revise.pdf`（极化图）和 `michelComponents.pdf`（π、K、K⁰L 三个独立来源的 Michel 极坐标图），以及对应的 ROOT、CSV 和溯源记录。程序在读取大文件前检查所有目标文件是否已存在，默认拒绝覆盖；重跑时请使用新目录。

第四、第五个可选参数分别是**每个独立成分的采样次数**和随机种子（不是组分比例的采样单位）：

```bash
bash src/analysis_pmu_from_MusAirS/run_pipeline.sh /path/proton.root /path/helium.root pmu_high_stats 100000000 20241025
```

默认使用逆 CDF 直接采样，每个成分采样 `100000000` 次（100M），三个来源合计 3 亿次，随机种子 `20241025`。上述示例显式指定默认统计量。试跑时可显式传入 `1000000` 降低耗时。三个成分始终独立处理、采样量相同，不使用 `100:15:5` 混合比例。

## 分别运行各阶段

### C++：从原始 ROOT 文件生成绘图数据

```bash
root -l -b -q 'src/analysis_pmu_from_MusAirS/analyze_pmu.C("/path/proton.root","/path/helium.root","pmu_output")'
```

ROOT 宏函数参数为：

```cpp
analyze_pmu(proton, helium, outputDirectory = "pmu_output",
            bins = 50, legacyDeduplicate = true, overwrite = false)
```

例如，改用当前 `src/getMichel.cc` 的“保留全部 μ 子”行为：

```bash
root -l -b -q 'src/analysis_pmu_from_MusAirS/analyze_pmu.C("/path/proton.root","/path/helium.root","pmu_all_muons",50,false)'
```

输出：

| 文件                    | 内容                                                        |
| ----------------------- | ----------------------------------------------------------- |
| `polarizationData.root` | 六个 TH1D，以及 `provenance` 元数据对象                     |
| `polarizationData.csv`  | `histogram,bin,x,y,yerr`，只保存常规 bin                    |
| `provenance.txt`        | 输入绝对路径、大小、ROOT 版本、筛选和去重规则、权重、bin 数 |

文件名不沿用 `polarizationData20.root`，避免把历史命名中的 `20` 误认为 bin 数。默认实际为 50 bins，范围 `[-1,1)`。

元数据不是内容校验和，也不自动记录 Git 提交。正式归档时请额外保存 `git rev-parse HEAD` 与原始文件 SHA-256，例如 `shasum -a 256 /path/proton.root /path/helium.root`。大文件的校验会额外读取整个文件。

### Python：从最终 CSV 一键绘图

```bash
python src/analysis_pmu_from_MusAirS/plot_pmu.py pmu_output/polarizationData.csv
```

默认在 CSV 所在目录保存 `combinedPmu_revise.pdf`。其他选项：

```bash
python src/analysis_pmu_from_MusAirS/plot_pmu.py pmu_output/polarizationData.csv \
  --output /path/figure.pdf --overwrite --show
```

默认使用无界面的 Agg 后端，适合远程机器。`--show` 才启用交互显示。支持 Matplotlib 可保存的其他输出格式，例如 `.png`。图例优先使用 Times New Roman；未安装时使用 DejaVu Serif。

## Michel 角分布生成和绘图

已有极化数据时，不需要重新读取两份大原始文件：

```bash
root -l -b -q 'src/analysis_pmu_from_MusAirS/generate_michel.C("pmu_output/polarizationData.root","pmu_output",100000000LL,20241025)'
python src/analysis_pmu_from_MusAirS/plot_michel.py pmu_output/michelData.csv
```

宏参数：

```cpp
generate_michel(input = "pmu_output/polarizationData.root",
                outputDirectory = "pmu_output", eventsPerComponent = 100000000LL,
                seed = 20241025, overwrite = false)
```

- 输入三个 μ− 极化直方图，各自按常规 bin 内容计算平均极化；不叠加、不按母粒子比例缩放。
- `piMuonPmu → piMuon`：π 来源，独立生成 `eventsPerComponent` 次。
- `KMuonPmu → KMuon`：K 来源，独立生成 `eventsPerComponent` 次。
- `KLongMuonPmu → KLongMuon`：K⁰L 来源，独立生成 `eventsPerComponent` 次。
- 对应 `dev/plot.ipynb` 中读取 `michelElectron.root` 并绘制 `KMuon`、`piMuon`、`KLongMuon` 的来源对比图，而不是读取 `combineData*.root` 的组分叠加图。
- 原权重满足 `a²+b²+c²=1`。对独立均匀辅助角及极化积分后，边缘密度为 `p(x)=(1+⟨Pmu⟩x/3)/2`，其中 `x=cos(theta)`。这与旧加权模型的边缘角分布期望一致，不是新物理模型。
- 平均极化用 bin 中心和常规 bin 内容计算，与原 `GetRandom()` 的 bin 内均匀采样一致；不使用可能保留原始填充统计的 ROOT `GetMean()`。
- 用 `TRandom3(seed)` 每事例抽一个随机数，通过稳定逆 CDF 直接得到 `cosθ`，仍实际生成指定数量的样本。使用 64 位 bin 计数器累积，避免逐事件 `TH1::Fill`。种子 0 被拒绝；相同输入、参数及 ROOT 版本可复现相同 CSV。
- `michelData.root` 保存三个 100-bin 独立 Michel `cosθ` TH1D 和溯源对象；输入极化仍保存在上游 `polarizationData.root`。CSV 保存**无权计数及 `sqrt(count)` 误差**，不逐曲线重归一化。每条直方图的填充次数及常规 bin 积分严格等于 `eventsPerComponent`。ROOT 角度矩统计使用 bin 中心。
- `michel_provenance.txt` 保存随机种子、每成分采样次数、来源映射、ROOT 版本和上游极化元数据，并明确记录 `mode=independent_components`。
- 生成的是 **Michel 电子边缘角分布**，不包含电子能谱、空间顶点、探测器接受度或分辨率，也不保留辅助变量的联合分布。若增加这些变量或相关筛选，必须重新推导；不能直接沿用此替换。

绘图支持 `--output`、`--overwrite`、`--show`，默认输出 `michelComponents.pdf`：

```bash
python src/analysis_pmu_from_MusAirS/plot_michel.py pmu_output/michelData.csv --output /path/michel.png
```

使用 `arccos` 将 bin 中心变成角度，再镜像到整圆。绘图径向值仍是 **`dN/dcosθ` 的 bin 计数**，并不是经过 Jacobian 变换的 `dN/dθ`。采用周期三次样条闭合曲线；低统计量下负的样条过冲会警告并截为零，原始 CSV 不改变，建议增大采样量。没有把统计误差绘制成误差带。

π、K、K⁰L 曲线分别为红、绿、蓝实线，图例与 μ− 来源一一对应。没有总组合曲线或黑色虚线。保存前统一设置径向范围，使用固定随机种子。

**迁移语义更正**：初版错误地迁移了 `dev/combine.cpp` 的逐项叠加流程。现已改为三个独立来源；原先生成的 `combinedMichel.pdf` 和历史验证报告只能说明旧叠加流程的执行情况，不能作为当前目标图的验证结果。新结果另存为 `michelComponents.pdf`。

**默认方法已替换为直接采样，默认每来源 100M 次**。与原加权算法的边缘角分布期望相同，但不逐事例复现其 RNG 序列，计数与误差含义也不同。固定随机种子及周期插值仍保留。极化数据处理部分不变，仍能精确复现已验证的历史六个直方图。

## 原加权算法对照与兼容入口

原算法保留为 `generate_michel_weighted.C`，仅用于验证和比较，不再由一键入口调用。对照结果必须放在另一个目录：

```bash
root -l -b -q 'src/analysis_pmu_from_MusAirS/generate_michel_weighted.C("pmu_output/polarizationData.root","michel_weighted_output",100000000LL,20241025)'
python src/analysis_pmu_from_MusAirS/compare_michel_methods.py pmu_output/michelData.csv michel_weighted_output/michelData.csv pmu_output/polarizationData.csv --output pmu_output/comparison.pdf
```

原算法输出加权计数及 `Sumw2` 误差；直接采样输出无权计数及 `sqrt(count)`。比较脚本按各自 provenance 中的采样量缩放，不按各曲线实际积分重归一化。上排显示两种方法与解析期望，下排显示以各自 bin 误差标准化的残差。

此前试验入口 `generate_michel_direct.C` 保留为调用主生成器的薄兼容包装，共用同一实现，默认也是每来源 100M，输出目录仍为 `michel_direct_output`。正在运行的旧进程及已有数据不会被替换或覆盖。

## 极化数据：与旧图保持一致的处理规则

依据 `src/getMichel.cc`、2024 年 Git 版本 `2620910` 和 `dev/plot.ipynb`：

1. 两个输入文件分别读取 `G4Run0/ReactionChain(13)` 与 `G4Run0/ReactionChain(-13)`。
2. 初筛：`PDGID == ±13`、`KillProc == "<0|"`、`Zenith < 1 - sqrt(3)/2`。 原始字段 `Zenith` 定义是 **`1 - cos(zenith angle)`**，因此这个条件对应天顶角小于 30°；不是直接比较角度数值。
3. `Helicity` 直接作为 `Pmu`。
4. 默认复现旧 `polarizationEnergyRelation()` 的去重：同一棵输入树中，每个 `(EvtID, ParTrkID)` 只保留**最后一个通过初筛**的候选，再做母粒子和动能筛选。不能先筛动能再去重，否则会改变旧算法结果。
5. 母粒子分类和 μ 子动能条件：

   | 分类 | `abs(ParPDGID)` | 动能条件（原始 `Ek` 单位，通常为 MeV） |
   | ---- | --------------- | -------------------------------------- |
   | π    | 211             | `100 < Ek < 500`                       |
   | K    | 321             | `100 < Ek < 500`                       |
   | K⁰L  | 130             | `0 < Ek < 1e9`                         |

1. 质子、氦直方图分别乘以 `9.704844263299673`、`100.0`，相加后，各自按常规 bin 内容之和归一化。
2. 误差为 ROOT `Sumw2` 的加权计数误差，再乘同一归一化因子；不额外引入归一化分母的统计误差。这与原直方图处理方式一致。
3. 输出六个直方图：`piMuonPmu`、`piAntiMuonPmu`、`KMuonPmu`、`KAntiMuonPmu`、`KLongMuonPmu`、`KLongAntiMuonPmu`。绘图使用三个 μ− 直方图。
4. 图保持原有 5×5 英寸、红/绿/蓝误差点、对数纵轴和 `1e-4..1` 范围。纵轴保留历史标签 `counts`，但实际含义是**归一化 bin 含量**。零 bin 不绘制，也不添加伪计数。

### 边界和局限

- **历史去重与当前实现不同**：2026 年 `findRelationPanM()` 改用按事件处理并保留全部候选。`legacyDeduplicate=false` 可切换到这个计数规则；两个模式可能产生不同结果。
- 本地真实原始文件中存在非连续重复的 `EvtID`。历史模式因此按**整棵树**去重，不能在事件编号变化时提前写入直方图。它只保存候选的绘图字段，内存随唯一 `(EvtID, ParTrkID)` 数量增长；每处理完一棵树释放。全候选模式直接填充直方图，内存占用基本固定。关闭 ROOT 隐式多线程以保证最后候选的顺序语义。
- 合并文件如果不同模拟批次重复使用事件编号，历史代码会把这些编号视为同一事件并覆盖候选。本脚本默认保留这种行为以便复现；这不代表该去重规则在物理上正确。需要全候选计数时使用 `legacyDeduplicate=false`，不能将该模式与旧结果混为一谈。
- `Pmu == +1` 进入 ROOT overflow，`Pmu == -1` 进入第一个 bin，沿用原行为。ROOT 输出保留 underflow/overflow，CSV 和归一化分母不包含它们；出现时打印警告。
- 空直方图保持全零并警告，不除以零。绘图要求三个 μ− 直方图非空且归一化。
- 输入缺失、字段缺失、非有限能量/极化值会报错。该入口针对上述旧 MusAirS 格式，不自动支持带事件权重的新格式；不要用固定 P/He 缩放替代新数据中的逐事件权重。
- 本机上述两份输入已验证能精确复现六个历史直方图的 bin 内容及误差；其他输入或全候选模式不能直接假定结果相同。字体和 Matplotlib 版本仍可能造成 PDF 渲染差异。
- 所有输入处理完成后才写结果；磁盘写入失败可能留下不完整输出目录，应换新目录重跑。

## 测试

```bash
python src/analysis_pmu_from_MusAirS/test_pipeline.py -v
```

测试在临时目录创建小型原始 ROOT 文件，检查六个极化直方图、权重、误差、历史去重（包含非连续重复事件）/全候选模式、输入错误、覆盖保护、CSV 校验，以及 Michel 权重表达式与旧模型的一致性、固定种子复现、周期插值和一键双图流程。共 12 项测试；新增直接采样与原加权算法的边缘分布一致性检验，兼容入口与主入口结果一致性检验；新增直接采样逆 CDF 边界/单调性、精确计数、误差、分布期望、固定种子复现及覆盖保护验证；额外验证修改 K 的输入极化只影响 K 曲线，不影响 π 或 K⁰L，防止再次误迁移为混合组分。缺少 Matplotlib/SciPy 时跳过对应绘图测试；缺少 ROOT 时跳过整个集成测试。不读取或修改 `dev` 的历史数据。
