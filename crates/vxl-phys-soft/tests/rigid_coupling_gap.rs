//! **粒子↔刚体耦合的已知缺口：钉住判据**（T1 第四/五片，2026-09-26）。
//!
//! **现状（实测）**：绳把落下的盒子**接住**了（前 60 tick 的动量确实被吃掉），但**十几~几十 tick 后
//! 会被弹飞**（+4 m/s 量级），之后自由落体 ⇒ 1800 tick 的长窗判据**过不了**。
//!
//! **已排除的成因（8 组二维扫描，别重复试）**：接触**柔度** `α_c ∈ {0, 1e-6, 1e-5, 1e-4, 1e-3}`
//! × 绳**内摩擦** `damping ∈ {1.0, 0.999, 0.995, 0.99}` —— **全部 8 组都弹飞**
//! （稳态窗摆幅 2000–4200 m）⇒ **不是刚度错配、也不是阻尼不足**。
//!
//! **缺口已闭合（2026-09-27，见 `docs/SURVEY-SOFT-CLOTH-AND-CONVERSION.md` §8.4.10）**：
//! 加了**位置口径回填** `Rope::body_dx`（被速度钳位压掉的那一份**只补位置、不补速度**），
//! 本判据的引擎按**门面同口径**吃它 ⇒ **长窗 1800 tick 盒子停在 y = +0.982**（短窗 0.9802）。
//! 机制：下沉 = **位置漏**（体每 tick 先按 `v·dt` 走过 `g·dt² = 2.72 mm`，而接触只取消接近速度）
//! ⇒ 补足 ~3.1 mm/tick 正好抵掉（实测 Σ|dx| = 5.6 m / 1800 tick）。
//!
//! **仍未闭的是门面侧**：门面也喂 `body_dx`，盒子能停住几百 tick，但之后仍会滑穿
//! （实测 `|ω| = 0`、无睡眠 ⇒ 剩下的就是 §8.4.8 那个**离散接触集**：接触集一跳，支撑就断）。
//! 那一刀落地前，门面级判据（`crates/vxl-phys/tests/rope_scene.rs`）还不能翻。
use vxl_phys_core::{interop::NoProviders, Quat, Shape, Vec3};
use vxl_phys_soft::{RigidProxy, Rope};

const DT: f32 = 1.0 / 60.0;
const G: Vec3 = Vec3::new(0.0, -9.81, 0.0);

/// **托住判据**（原"钉住缺口"，2026-09-27 缺口闭合后**翻过来** ⇒ 现在判"**停在绳上**"）：
/// 紧绳 + 0.6 宽盒 + 1800 tick ⇒ 短窗（60 tick）在绳上、**长窗仍停在 y ≈ 0.98**。
#[test]
fn box_on_rope_is_held() {
    let mut r = Rope::line(
        Vec3::new(-0.5, 1.0, 0.0),
        Vec3::new(0.5, 1.0, 0.0),
        33,
        0.02,
    );
    r.damping = 0.999;
    for _ in 0..600 {
        r.step(DT, G, &NoProviders, 0, &[]);
    }
    let shape = Shape::Box {
        half: Vec3::new(0.3, 0.05, 0.3),
    };
    let (mut y, m) = (1.2f32, 1.0f32);
    let mut v = 0.0f32;
    let mut y_at_60 = 0.0f32;
    for t in 0..1800 {
        let proxy = RigidProxy {
            body: 0,
            shape,
            pos: Vec3::new(0.0, y, 0.0),
            rot: Quat::IDENTITY,
            linvel: Vec3::new(0.0, v, 0.0),
            inv_mass: 1.0 / m,
        };
        r.step(DT, G, &NoProviders, 0, std::slice::from_ref(&proxy));
        v += G.y * DT;
        y += v * DT;
        v += r.body_dv.first().map(|d| d.y).unwrap_or(0.0);
        // **位置口径回填**（`Rope::body_dx`，§8.4.10）：门面也这么做 ⇒ 自扮引擎必须同口径才可比。
        // 只回速度 = 体每 tick 按 `v·dt` 走过的 `g·dt²` 一去不回（"缓慢下沉"的真因）。
        y += r.body_dx.first().map(|d| d.y).unwrap_or(0.0);
        if t == 59 {
            y_at_60 = y;
        }
    }
    println!("短窗(60 tick) y={y_at_60:.4} | 长窗(1800 tick) 末 y={y:.4} v={v:+.3}");
    assert!(
        y_at_60 > 0.70,
        "短窗该被**接住**（阈值 0.70；实测 0.98）——短窗红了说明耦合本身坏了"
    );
    assert!(
        y > 0.5,
        "长窗该**停在绳上**（实测末 y={y:.4}、v={v:+.3}）——红了说明位置口径回填被改坏。\
         注：这条原先是钉住缺口用（当年末 y ≈ −2091）⇒ 2026-09-27 缺口闭合后翻成托住"
    );
}
