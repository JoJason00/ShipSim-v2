# ShipSim v2

船舶耐波性与波浪中操纵性计算程序。本仓库是 [ShipSim 参考实现](https://github.com/JoJason00/ShipSim-1.0)
的重构版本。

重构方案见 `../ShipSim-Refactor/`:
[总纲](../ShipSim-Refactor/README.md) ·
[基线策略](../ShipSim-Refactor/testing-strategy.md) ·
[目标架构](../ShipSim-Refactor/architecture.md) ·
[配置设计](../ShipSim-Refactor/config-design.md)

## 当前状态

| 步骤 | 状态 |
|---|---|
| Step 0 老代码上 git 并冻结(tag `reference-before-refactor`) | ✅ |
| Step 1 基础设施骨架 | ✅ |
| Step 2 搬入参考代码,依赖与编译选项保持不变 | ✅ 可编译,5 个可执行文件全部生成 |
| Step 3 检查点:结果与参考实现一致 | ⬜ **下一步** |
| Step 4 建立回归基线 | ⬜ |
| Step 5 依赖换 vcpkg(Eigen / jsoncpp → toml++) | ⬜ |
| Step 6 复核 `/fp:fast` | ⬜ |
| Step 7 拆模块、拆上帝类 | ⬜ |

### Step 3 开始前待定的两件事

1. **线程可复现性** —— 启动时报告 `OpenMP=30, Eigen=30`。并行归约的浮点
   累加顺序可能逐次不同。**先在参考实现上把同一算例连跑两遍**,确认结果
   是否逐位一致。若否,基线要么固定线程数,要么放宽容差。这比直接比新旧
   更基础,应当先做。
2. **选哪个算例** —— 建议先用 `wigleyI`(规模最小)。需要知道它大致耗时。

> **纪律:任何时刻 `ctest` 都必须是绿的。** 不允许出现"重构中,暂时跑不了"的状态。

## 环境

CMake ≥ 3.21、Ninja,以及指向 vcpkg 检出目录的 `VCPKG_ROOT`。

## 构建

任意终端下,`build.bat` 会自动找到并初始化 MSVC 环境:

```bat
build.bat            配置 + 编译 + 测试(debug,默认)
build.bat release    优化构建
build.bat asan       编译 + 测试,查内存错误
```

可执行文件在 `build/<preset>/bin/`。

在 *Developer PowerShell for VS* 里也可直接用 preset(等价):

```bash
cmake --preset debug && cmake --build --preset debug && ctest --preset debug
```

| Preset | 用途 |
|---|---|
| `debug` | 日常开发、测试 |
| `release` | 优化构建 |
| `profile` | 性能分析(`scripts/profile.ps1`) |
| `coverage` | 覆盖率(`scripts/coverage.ps1`) |
| `asan` | AddressSanitizer,抓内存错误 |

> ASan 版慢 2~3 倍,只在排查内存 bug 时用。查算例内存问题:
> `build/asan/bin/shipsim_cli.exe <case>`。日常用 `debug` / `release`。

## 测试分三层

```
tests/unit/         单个函数/类,硬编码期望值
tests/validation/   对标试验或解析解,宽容差 2~5%
tests/regression/   对标冻结基线,严容差 1e-10
```

三者用途不同,容差不可混用 —— 详见
[testing-strategy.md](../ShipSim-Refactor/testing-strategy.md)。

## 目录

```
src/               参考实现代码(Phase 3 拆分为各模块)
app/               shipsim_cli
test/              旧驱动程序(core_test 等,非单元测试)
tests/             新测试:unit / validation / regression
eigen-5.0.0/ json/ jsoncpp.cpp   vendored 依赖(Step 5 前保留)
cases/             算例输入(结果产物由 .gitignore 排除)
cmake/ scripts/    覆盖率、profiling 辅助
```

模块划分与依赖规则见 [architecture.md](../ShipSim-Refactor/architecture.md)。
