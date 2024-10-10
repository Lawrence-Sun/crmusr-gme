# GME

GME代表的是米歇尔电子产生子，在确定缪子极化的情况下生成米歇尔电子的空间分布。本仓库的结构为
- dev：主要用于装载测试使用的代码以及中间文件在git的过程中记得将文件夹中的代码加入.gitignore文件中
- include：主要用于装载头文件
- lib：主要用于装载头文件 <!-- 虽然现在还有什么头文件 -->
- src：主要用于装载源文件

## `getMichel.cc`
此源文件中有两个可以直接调用使用的类，分别是`getPmu`和`GenMichel`。

### `getPmu`
此类是从原始数据中产生$\mathcal{P}_\mu$的类。几个成员函数主要是
- `polarizationEnergyRelation`：主要用于从MusAirS的原始数据中提取出$\mu^\pm$的数据
- `FillTheTTree`：定义直接将对应数据放进`TTree`中的操作
- `findRelationPandM`：综合上述两个函数从原属数据中提取出缪子的事例以及对应的数据

- `getPAndHeHist`：从$\mu^\pm$的数据中提取出proton和helium产生的模拟数据
- `TH1DDrawing`：一个简单的TH1D的绘图规范
- `findCertainEnergyRangePolarization`：找到在固定能量范围内所有的$\mu^\pm$的极化数据

- ==`getPolarizationFromRaw`==：综合所有的功能直接从原始数据中提取固定立体角固定能量区间不同母粒子对应的缪子的极化分布。

### `GenMichel`
此类是根据$\mathcal{P}_\mu$生成米歇尔电子的类。几个成员函数主要是
- `getWeight` 主要是权重的生成函数，生成的依据是旋转后的Michel电子空间分布
- `generator` 均匀生成Michel电子的空间位置并且回去这个空间位置Michel电子的权重大小
- `PmuMichelDis` 通过制定的$K-\pi$ratio生成指定事例数的Michel电子空间分布
- `IsotropyDis` 生成指定事例数的各向同性均匀分布
- `PmuMichelDis2D` 生成指定$K-\pi$ratio下的Michel电子空间分布的二位直方图

具体的使用方法可以参考头文件最后的示例。$$

## How to use?
```c++
#include "the/path/to/getMichel.cc"

// example usage
int main()
{
    getPmu pmu; // build a object to call the class

    pu.findCertainEnergyRangePolarization(......);

    return 0;
}
```