//! **软体域（绳索）通道** + 2b 边界的**成组状态**。
//!
//! **为什么 `FluidBoundary` 定义在这里**：`World` 是 god 门棘轮下的**记录型**结构（成员数只准减），
//! 加一个域就得腾一个成员位 ⇒ 把 2b 那一族的 3 个散字段收成一个结构（成员 23 → 22，再 +1 给绳索）。
//! 定义随域走、`World` 只留一个字段。
//!
//! **绳索通道**（T1，`docs/SURVEY-SOFT-CLOTH-AND-CONVERSION.md` §8）：与液体域同款——每 tick **一次**
//! （体子步全部完成之后）推进；接触走**统一提供者通道**（`Providers` 实现 `ProviderColliders`，
//! 与窄相同一个 id 空间，见 `world_step.rs` 的 `fluid_pass` 同款段位）。**没有绳索的场景逐位不变**
//! （`ropes` 空 ⇒ `rope_pass` 首行短路）⇒ 默认档判据不受影响。
//!
//! **本片边界**：绳索只读**提供者**（地形）；与刚体的双向耦合（Akinci 边界）、摩擦（切向）、
//! 自碰撞都属后续切片。
use super::*;

/// 2b（Akinci 边界粒子）那一族的成组状态（原为 `World` 的 3 个散字段：开关 / 暂存 / 覆盖集）。
#[derive(Default)]
pub struct FluidBoundary {
    /// 与 `fluids` 同序的 2b 开关：`true` = 该流体每 tick 重建边界粒子并回流反作用。
    pub two_b: Vec<bool>,
    /// 边界粒子生成的暂存 `(体 id, 形状, 位姿)`（复用免每 tick 分配）。
    pub scratch: Vec<(u32, Shape, vxl_phys_fluid::BodyPose)>,
    /// **覆盖集**（与 `bodies` 同序，每 tick 重建）：上次进了边界粒子集的体。
    pub covered: Vec<bool>,
}

impl FluidBoundary {
    /// 该流体是否开了 2b（越界一律 `false`）。**写成方法而不是在调用点展开链**：
    /// 展开式在原地超 `chain_width` 会被 rustfmt 折成 5 行 —— 而 `world_body.rs` 受尺寸棘轮
    /// （只准减），折行会让它"变胖"。
    pub fn is_two_b(&self, fluid: usize) -> bool {
        self.two_b.get(fluid).copied().unwrap_or(false)
    }
}

impl World {
    /// 注册一条绳索（`vxl_phys_soft::Rope`），返回其索引。
    pub fn add_rope(&mut self, rope: vxl_phys_soft::Rope) -> usize {
        self.ropes.push(rope);
        self.ropes.len() - 1
    }

    /// 已注册绳索（判据/渲染读 `pos` / `vel`）。
    pub fn ropes(&self) -> &[vxl_phys_soft::Rope] {
        &self.ropes
    }

    /// 第 `i` 条绳索。
    pub fn rope(&self, i: usize) -> Option<&vxl_phys_soft::Rope> {
        self.ropes.get(i)
    }

    /// **域通道**（每 tick 一次、体子步全部完成之后）：液体域 → 软体域。
    /// 顺序固定 ⇒ 确定性不受影响（两条通道互不读对方状态）。
    pub(crate) fn domain_pass(&mut self) {
        self.fluid_pass();
        self.rope_pass();
    }

    /// **软体域通道**：每条绳索按自身 `substeps` 推进一个 `config.dt`；接触走统一提供者通道
    /// （`0..providers.len()` 全量 id）。**空集 ⇒ 零成本短路**。
    pub(crate) fn rope_pass(&mut self) {
        if self.ropes.is_empty() {
            return;
        }
        let dt = self.config.dt;
        let gravity = self.config.gravity;
        let count = self.providers.len() as u32;
        // 两个字段互不重叠 ⇒ 可同时借（`&mut self.ropes` + `&self.providers`）。
        let providers = &self.providers;
        for rope in &mut self.ropes {
            rope.step(dt, gravity, providers, count);
        }
    }
}
