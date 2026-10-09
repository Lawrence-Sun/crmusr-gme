# 从 MusAirS 原始数据生成 μ 子极化图

本目录保存可纳入 Git 的入口、处理规则、绘图脚本和测试，不依赖 `dev/` 的代码或中间数据。

```
质子 + 氦 MusAirS ROOT 文件
  → analyze_pmu.C：初筛、历史去重、母粒子/动能筛选、分量合并、归一化
  → polarizationData.root + polarizationData.csv + provenance.txt
  → plot_pmu.py
  → combinedPmu_revise.pdf
```

这里直接生成最终绘图所需的六个直方图，不生成 `muonData.root`；母粒子动能 `ParEk` 不参与这张图的筛选，因此不需要提取或保存。

## 依赖

- ROOT 6，支持 C++17 / `ROOT::RDataFrame`；本地测试版本为 6.38.02。
- Python 3.9+ 和 Matplotlib；Python 端只读取 CSV，不需要 PyROOT 或 uproot。

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

输出目录默认 `pmu_output`。程序默认拒绝覆盖已有数据，重跑时请使用新目录。已有 PDF 也拒绝覆盖；完整一键入口不提供强制覆盖选项。

## 分别运行两个阶段

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

## 与旧图保持一致的处理规则

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

测试在临时目录创建小型原始 ROOT 文件，检查六个直方图、权重、误差、历史去重（包含非连续重复事件）/全候选两种模式、输入错误、覆盖保护、CSV 校验和 PDF 绘图。缺少 Matplotlib 时只跳过绘图测试；缺少 ROOT 时跳过整个集成测试。不读取或修改 `dev` 的历史数据。
