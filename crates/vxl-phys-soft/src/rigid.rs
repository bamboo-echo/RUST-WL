//! **粒子 ↔ 刚体耦合**（`SPEC.md` §4.6「粒子-刚体距离约束（Akinci 式边界处理）」的最小实现）：
//! 刚体**代理视图** + 穿透查询 + **反作用回填**。
//!
//! **为什么自带查询而不是用窄相**：软体 crate 在依赖图上只挂 `core`（`soft <- 门面`）——
//! 引窄相会多一条边，而这里只需要"**球 vs 形状**"这一种查询（粒子 = 半径 `radius` 的球）。
//! 本片只做 **Sphere / Box / Capsule**（解析、共 ~60 行）；`Cylinder`/`Cone`/`ConvexHull`/`Compound`
//! 返回 `None`（写清边界，属后续切片——与窄相那条"六形状 × 四提供者"是两码事）。
//!
//! **反作用口径**（与 2b 流体反作用同段位）：软体侧累加**冲量**（粒子的动量变化取反），
//! 门面 `÷dt` 后按**力**加到体上（引擎每子步施加一次 ⇒ 一个 tick 的冲量 = `F·dt`，账平）。
use vxl_phys_core::{Mat3, Quat, Shape, Vec3};

/// 刚体代理：门面每 tick 填一次，软体侧只读（含接触所需的最小字段）。
pub struct RigidProxy {
    /// 门面的体索引（反作用回填用）。
    pub body: u32,
    pub shape: Shape,
    pub pos: Vec3,
    pub rot: Quat,
    /// 体心线速度（摩擦用：粒子滑移量取**相对体**的，体在动时绳才不会被"粘"在原地）。
    pub linvel: Vec3,
    /// `0` = 静态（仍参与接触，但不接收反作用）。
    pub inv_mass: f32,
}

/// 反作用（**冲量**口径，每 tick 累加；门面 `÷dt` 后作为力施加）。
#[derive(Clone, Copy)]
pub struct RigidReaction {
    pub body: u32,
    /// 作用在体上的冲量（= 粒子所受冲量取反）。
    pub impulse: Vec3,
    /// 绕**体原点**的角冲量（接触点 × 冲量）。
    pub torque: Vec3,
}

/// 球（心 `p`、半径 `radius`）对 `shape`（位姿 `pos`/`rot`）的**穿透**：
/// 返回 `(外向法线, 穿透深度, 体表接触点)`；不接触 / 形状不支持 ⇒ `None`。
///
/// 法线约定与提供者通道一致：**从体表面指向粒子**（`pos += n·depth` 即推出）。
pub fn shape_penetration(
    shape: &Shape,
    pos: Vec3,
    rot: Quat,
    p: Vec3,
    radius: f32,
) -> Option<(Vec3, f32, Vec3)> {
    match *shape {
        Shape::Sphere { radius: r } => {
            let d = p - pos;
            let dist = d.length();
            if dist < 1e-9 {
                return Some((Vec3::Y, r + radius, pos)); // 同心退化：任取方向
            }
            if dist > r + radius {
                return None;
            }
            let n = d * (1.0 / dist);
            Some((n, r + radius - dist, p - n * radius))
        }
        Shape::Box { half } => {
            let m = Mat3::from_quat(rot);
            // 局部系里做"点 vs 盒"（夹取求最近点），再变换回去。
            let local = m.transpose_mul_vec3(p - pos);
            let c = Vec3::new(
                local.x.clamp(-half.x, half.x),
                local.y.clamp(-half.y, half.y),
                local.z.clamp(-half.z, half.z),
            );
            let d = local - c;
            let dist = d.length();
            if dist > 1e-6 {
                if dist > radius {
                    return None; // 外侧且够远
                }
                let n_local = d * (1.0 / dist);
                return Some((m.mul_vec3(n_local), radius - dist, pos + m.mul_vec3(c)));
            }
            // 点在盒内：法线取**最近面**的外向，深度 = 到该面 + 半径。
            let dx = half.x - local.x.abs();
            let dy = half.y - local.y.abs();
            let dz = half.z - local.z.abs();
            let (deep, axis) = if dx <= dy && dx <= dz {
                (dx, Vec3::new(local.x.signum(), 0.0, 0.0))
            } else if dy <= dz {
                (dy, Vec3::new(0.0, local.y.signum(), 0.0))
            } else {
                (dz, Vec3::new(0.0, 0.0, local.z.signum()))
            };
            let n = m.mul_vec3(axis);
            Some((n, deep + radius, p - n * radius))
        }
        Shape::Capsule {
            half_height,
            radius: r,
        } => {
            let m = Mat3::from_quat(rot);
            let axis = m.mul_vec3(Vec3::Y);
            // 点到轴段的最近点（夹取投影参数）。
            let t = (p - pos).dot(axis).clamp(-half_height, half_height);
            let q_axis = pos + axis * t;
            let d = p - q_axis;
            let dist = d.length();
            let total = r + radius;
            if dist < 1e-9 {
                return Some((Vec3::Y, total, q_axis));
            }
            if dist > total {
                return None;
            }
            let n = d * (1.0 / dist);
            Some((n, total - dist, p - n * radius))
        }
        _ => None, // Cylinder / Cone / ConvexHull / Compound / Provider / HeightField：待补
    }
}
