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
| Step 1 基础设施骨架 | ✅ 本仓库 |
| Step 2 搬入参考代码,依赖与编译选项保持不变 | ⬜ |
| Step 3 检查点:结果与参考实现一致 | ⬜ |
| Step 4 建立回归基线 | ⬜ |
| Step 5 依赖换 vcpkg(Eigen / jsoncpp → toml++) | ⬜ |
| Step 6 复核 `/fp:fast` | ⬜ |
| Step 7 拆模块、拆上帝类 | ⬜ |

> **纪律:任何时刻 `ctest` 都必须是绿的。** 不允许出现"重构中,暂时跑不了"的状态。

## 环境

CMake ≥ 3.21、Ninja,以及指向 vcpkg 检出目录的 `VCPKG_ROOT`。

**Windows:** 从 *x64 Native Tools Command Prompt* 或 *Developer PowerShell for VS*
运行,否则 Ninja 找不到 `cl.exe` 会静默选中 MinGW —— 构建会直接报错拦住。

## 构建

```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
```

| Preset | 构建类型 | 用途 |
|---|---|---|
| `debug` | Debug | 日常开发、测试 |
| `release` | Release | 优化构建 |
| `profile` | RelWithDebInfo | 性能分析(`scripts/profile.ps1`) |
| `coverage` | Debug + clang-cl | 覆盖率(`scripts/coverage.ps1`) |
| `asan` | Debug + AddressSanitizer | 抓内存错误 |

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
src/core/         数学 / 类型 / 常量
third_party/      Eigen、jsoncpp(Step 5 前保持 vendored)
tests/            unit / validation / regression
cmake/            enable_coverage() 等辅助
scripts/          coverage.ps1  profile.ps1
```

模块划分与依赖规则见 [architecture.md](../ShipSim-Refactor/architecture.md)。
