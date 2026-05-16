# agri-ackermann-planner

农机阿克曼底盘路径重规划系统 — 将直角弓字形路径转换为符合阿克曼转向约束的可执行路径。

## 核心功能

### 1. 阿克曼运动学约束路径重规划
- **自行车模型**：δ = arctan(L/R)，最小转弯半径约束（R_min=5m，L=2.5m → δ_max=26.57°）
- **曲率法检测急转弯**：κ = |x'y'' - y'x''| / (x'²+y'²)^(3/2)，κ > 1/R_min 触发圆弧过渡
- **两圆交点法构造圆弧**：沿 heading_entry/exit 方向延伸 R_min 得切点，弦中点法处理 U 型掉头（|Δheading - π| < π/6）
- **幅宽衔接校准**：消除作业带间隙/重叠，保证连续覆盖

### 2. 地形适应性分析（RTK 高程 → 坡度 → 限速）
| 坡度 | 限速 | 标签颜色 |
|------|------|---------|
| < 5° | 6 km/h | 绿色 |
| 5° ~ 15° | 4 km/h | 橙色 |
| > 15° | 2 km/h | 红色 |

### 3. 作业质量评估
- 路径索引格网法覆盖率计算（O(n×k)，k≈126，无 scipy 依赖）
- 漏喷/重喷热力图可视化
- 有效覆盖率 ECR 统计

## 快速开始

```bash
# 安装依赖（仅 numpy + matplotlib，scipy 可选）
pip install numpy matplotlib

# 运行演示（自动生成对比图 + 坡度图 + 覆盖质量图）
python ackermann_path_planner.py
```

## 输出文件

```
demo_results/
├── ackermann_comparison.png   # 原始直角 vs 阿克曼圆弧路径对比
├── terrain_slope_map.png      # 坡度热力图 + 限速标注
├── coverage_quality.png       # 覆盖率热力图（绿=正常/红=漏/蓝=重）
└── optimized_ackermann_path.csv # 含航向角路径（ROS 接入用）
```

## 算法参数（可调）

```python
planner = AckermannPathPlanner(
    min_turning_radius=5.0,   # 最小转弯半径 R_min (m)
    wheelbase=2.5,            # 轴距 L (m)
    max_steering_angle=None,  # 最大转向角（None=自动由R_min计算）
    working_width=3.0         # 作业幅宽 (m)
)
```

## 适用场景

- 农机自动驾驶路径跟踪（前馈 + RTK 反馈校正）
- 不规则田块弓字形路径生成
- 多机协同作业路径规划基础模块
- ROS/Apollo 路径跟踪接口

## 技术栈

Python 3 / NumPy / Matplotlib / SciPy（可选）

## 许可证

MIT
