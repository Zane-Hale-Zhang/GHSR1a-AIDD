# 实验结果记录 · p-AMPK 时程（实验 A）

> 完成日期：2026-09-11
> 记录日期：2026-09-12
> **决策点 Day 7：✅ Go —— 进入 W2 机制主线**

---

## 一、实验条件

| 项目 | 内容 |
|---|---|
| 细胞 | H9C2，未分化，DMEM + 10% FBS |
| 处理 | acyl-ghrelin 时程 |
| 分组 | 0 / 15 / 30 / 45 / 60 / 90 min + AICAR + Compound C，共 8 组 |

---

## 二、结果

| 组别 | p-AMPK 条带 |
|---|---|
| **AICAR（阳性对照）** | **明显** ✅ |
| 0 min（基线） | ⚠️ **结果待确认** |
| 15 min | 明显 |
| 30 min | 明显 |
| 45 min | 明显 |
| **60 min** | **最浓（峰值）** |
| 90 min | 明显（弱于 60 min） |
| Compound C | ⚠️ **结果待确认** |

---

## 三、可以得出的结论

### 1. ✅ WB 体系完全正常

AICAR 出现明显条带，证明**抗体、转膜、上样、磷酸酶抑制剂全部有效**。
这是整套实验的前提，现在已确认。

### 2. ✅ acyl-ghrelin 在 H9C2 上诱导 AMPK Thr172 磷酸化

与 Chen et al. 2019（ghrelin 保护 H9c2 抗 H/R，经 AMPK 通路）结论一致。

### 3. ✅ 时程动力学明确：15 min 起效，60 min 达峰

**后续所有实验的 ghrelin 处理时间统一用 60 min。**

90 min 弱于 60 min 提示信号是双相的（上升后回落），符合 GPCR 信号特征。

### 4. 🔑 最关键的一条：H9C2 表达**功能性** GHSR1a

acyl-ghrelin **只能通过 GHSR1a 起作用**（des-acyl ghrelin 不激活该受体）。
能在 H9C2 上看到时间依赖的 AMPK 磷酸化，说明：

- 受体**存在**
- 受体**有功能**
- 下游信号通路**完整**

**这是功能学证据，比 qPCR 测 mRNA 更有说服力。**
课题立论的根基（"H9C2 表达 GHSR1a"）已经由这个实验间接证实。

---

## 四、还不能得出的结论

| 待证问题 | 缺什么 | 何时补 |
|---|---|---|
| 效应是否**GHSR1a 特异** | GHSR 拮抗剂（[D-Lys³]-GHRP-6）或 GHSR siRNA | Day 11–12 |
| 效应是否**AMPK 依赖** | **Compound C 预处理 + ghrelin 共处理组** | Day 9–10 |
| 效应量是否**统计显著** | 灰度定量 + n≥3 独立重复 | 现在就做 |

⚠️ **"条带浓"是定性描述，不能进论文。**
必须做 densitometry，算 **p-AMPK / total AMPK 比值**，做 3 次独立重复后统计。

---

## 五、两个必须澄清的关键信息

### 1. 0 min 组到底有没有条带？

这决定了论文的表述方式：

| 0 min 情况 | 含义 | 论文写法 |
|---|---|---|
| **无条带或极弱** | ghrelin **诱导**（从无到有） | "ghrelin **induces** AMPK phosphorylation" |
| **有条带，但明显浅于处理组** | ghrelin **增强**（基础磷酸化上调） | "ghrelin **increases / enhances** AMPK phosphorylation" |

两种都说得通，但**不能混用**，否则会被审稿人抓表述不严谨。

### 2. Compound C 组是怎么设计的？

**单独一组 Compound C 不能证明 ghrelin 的效应是 AMPK 依赖的。**

必须这四组齐备：

```
Vehicle
ghrelin (60 min)              ← 效应组
Compound C alone              ← 排除抑制剂自身效应
Compound C 预处理 + ghrelin    ← ★ 关键，阻断组
```

只有最后一组条带明显变浅/消失，才能说
"ghrelin 诱导的 AMPK 磷酸化是 AMPK 通路依赖的"。

如果你只做了 Compound C alone，需要补做共处理组。

---

## 六、下一步（按优先级）

### 立即做

- [ ] **ImageJ / Image Lab 灰度定量**：p-AMPK / total AMPK 比值
- [ ] 确认 **total AMPK 各组是否一致**（应该不变；若变则说明上样有问题）
- [ ] 确认内参（GAPDH/β-actin）是否均匀

### 本轮重复（n=3）

- [ ] 同一时程再重复 **2 次**（共 3 次独立实验）
- [ ] 建议补 **5 min** 和 **120 min** 两个点，把曲线两端封住
- [ ] 60 min 这个点每次必做

### 进入 W2

- [ ] Day 9–10 流式阻断核心组：**ghrelin + Compound C** / **ghrelin + SSO**
- [ ] ghrelin 处理时间统一 **60 min**
- [ ] Day 8 并行启动 NRVM

---

## 七、对论文的直接价值

| 图 | 内容 | 状态 |
|---|---|---|
| 主图 A | ghrelin 时程诱导 p-AMPK（15–90 min，60 min 峰） | ✅ 数据已有，需定量+重复 |
| 主图 B | Compound C 阻断 | 待补共处理组 |
| 主图 C | CD36 转位/摄取（SSO 阻断） | Day 9–10 |

---

## 备注

- 时程实验已过 Day 7 决策点，**结果为 Go**，无需启动应急方案 A
- 60 min 达峰偏慢（AMPK 磷酸化常在 5–30 min 达峰），
  可能提示 ghrelin 经 Ca²⁺–CaMKKβ 间接激活，而非直接快速磷酸化。
  这一猜想可用 STO-609（CaMKK 抑制剂）在 Day 11–12 检验。
