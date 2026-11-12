# 量衡测量处理

将仪器原始观测按校准证书与质量策略处理，生成时间窗统计及可校验的测量发布包。

Python 3.11 标准库与 SQLite 即可运行。数据、时钟和文件均在本地，不需要服务凭据。示例只使用合成记录。

## 运行

在仓库根目录执行：

```sh
python -m metrology demo
python -m examples.context
python -m unittest discover -s tests -v
python -m compileall -q metrology examples
```

示例使用内存或临时目录，不会改动持久台账。`examples/context.json` 记录一组可重复导入的业务输入，字段与时间口径见 `examples/数据说明.md`。`examples.context` 检查数据并调用现有产品入口；它只展示原有业务数据，不计算分析结论。

`tests` 覆盖现有业务生命周期、非法输入、版本和完整性边界，另有示例导入回归。完整测试和上述演示可在一次短终端会话中执行。

## 方法一致性研究

便携法替代台架法不能只看相关系数：两种方法可以高度相关却始终相差一个固定或比例偏差。`metrology.agreement` 以同一试件的成对观测为单位做 Bland–Altman 方法一致性研究，并在提供误差方差比时另做 Deming 回归。

研究在 JSON 动作文档中固定两侧方法（参考/候选仪器与通道）、试件标签、重复测量策略（`first`、`last`、`mean`）、有效时间窗、配对最大时间间隔、预先确定的允许差限和预设量程。配对按试件键与时间就近匹配，绝不按行号硬配；未成对、超出配对时间窗、同一观测重复占用、等距歧义、缺失试件键和缺失/撤销证书都在 `exclusions` 中逐条列明。一个试件的多次重复先按策略汇总为一个独立单位，bootstrap 也按试件整体重采样，重复读数不会被当成多个独立试件。

输出包含读数级与试件级 Bland–Altman 差值（偏差、差值标准差、95% 一致性限）、均值相关偏差（差值对均值的回归斜率）、重复测量下的组内方差与重复性限，以及按试件整体重采样的一致性限区间。差值尺度可选绝对差、预设百分比（需正值参考值）或自然对数比（需两侧均为正值）；零值和负值不适用时给出明确诊断而不是静默产出。提供正的误差方差比时执行 Deming 回归并返回斜率、截距及刀切法区间，不用普通最小二乘冒充两侧均有误差的换算；方差比非法或数据退化时只给可支持的结论。

替代判断把预设量程内的一致性限及其 bootstrap 区间与允许差限比较：限和区间全部落在允许差限内且样本与量程覆盖充分才为 `interchangeable`；样本不足（少于 8 个试件）或量程覆盖不足时降为 `conditionally_interchangeable` 并附原因，超限为 `not_interchangeable`，无法计算为 `inconclusive`。研究固化配对、观测版本和证书指纹；后续校准更正后用重放生成新观测版本，再以新研究 ID 建版本，通过版本比较查看哪些试件的结论改变。

```json
{
  "action": "agreement_run",
  "specification": {
    "id": "agreement-2026-01",
    "reference": {"instrument": "bench", "channel": "p"},
    "candidate": {"instrument": "portable", "channel": "p"},
    "start": "2026-11-01T00:00:00Z",
    "end": "2026-12-01T00:00:00Z",
    "specimen_tag": "specimen",
    "replicate_strategy": "mean",
    "pairing_seconds": 5,
    "acceptable_difference": "1.5",
    "measurement_range": {"low": "5", "high": "100"},
    "scale": "absolute",
    "bootstrap_replicates": 2000,
    "bootstrap_seed": 7,
    "variance_ratio": "1"
  }
}
```

相关动作还有 `agreement_get`、`agreement_list`、`agreement_compare`（两个研究版本）和 `agreement_export`（`json`、`pairs`、`exclusions`、`summary` 四种导出，复用现有原子导出）。`tests/test_agreement.py` 用高相关固定偏差、比例偏差、重复测量、缺失配对、撤销证书、零值/负值尺度诊断和小样本案例核对结果。

## 容器

```sh
docker build -t local/metrology .
docker run --rm --network none local/metrology
docker run --rm --network none local/metrology python -m examples.context
```

镜像构建完成后运行和测试不下载文件。持久化使用各自命令的本地路径参数；具体操作及参数可通过 `python -m metrology --help` 查询。第三方依赖清单为空。
