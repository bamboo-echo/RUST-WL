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

/// **钉住缺口（3 维口径，2026-09-27 重立）**：紧绳 + 0.6 宽盒 + 1800 tick ⇒ 短窗（60 tick）接住
/// （0.9968），**长窗横向逃逸**（y = −2599.7、x = +141.4）。
/// ⚠️ 这一条曾经被翻成"托住"——那是**1 维自扮引擎**的假象（只看 `y`、只吃反作用的 `y` 分量
/// ⇒ 摩擦的**横向**分量被丢掉 ⇒ 盒子不可能横向滑出）。升到 3 维后**门面与自扮引擎结论一致**
/// （都托不住）⇒ 原先"门面 vs 自扮差 300×"整条线索是**仪器维数**造成的（§8.4.16）。
/// 修好后（横向滑出被治住）把断言翻成 `y > 0.5`。
#[test]
fn box_on_rope_is_caught_but_slips_sideways_pinned_gap() {
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
    // **三维自扮引擎**（§8.4.16）：`pos`/`linvel` 是全 `Vec3`、反作用吃**整个向量**
    // （`v += dv`、`x += dx`）。**原来是 1 维的**（`pos = Vec3::new(0.0, y, 0.0)`、只吃 `.y`）
    // ⇒ **摩擦反作用的横向分量被静默丢掉** ⇒ 与 3 维门面对拍时混进"1D vs 3D"这一整类差异。
    let shape = Shape::Box {
        half: Vec3::new(0.3, 0.05, 0.3),
    };
    let m = 1.0f32;
    let (mut pos, mut vel) = (Vec3::new(0.0, 1.2, 0.0), Vec3::ZERO);
    let mut pos_at_60 = Vec3::ZERO;
    for t in 0..1800 {
        let proxy = RigidProxy {
            body: 0,
            shape,
            pos,
            rot: Quat::IDENTITY,
            linvel: vel,
            inv_mass: 1.0 / m,
        };
        r.step(DT, G, &NoProviders, 0, std::slice::from_ref(&proxy));
        vel += G * DT;
        pos += vel * DT;
        if let Some(dv) = r.body_dv.first() {
            vel += *dv; // 速度口径（整向量：含摩擦的横向分量）
        }
        // **位置口径回填**（`Rope::body_dx`，§8.4.10）：门面也这么做 ⇒ 自扮引擎必须同口径才可比。
        // 只回速度 = 体每 tick 按 `v·dt` 走过的 `g·dt²` 一去不回（"缓慢下沉"的真因）。
        if let Some(dx) = r.body_dx.first() {
            pos += *dx;
        }
        if t == 59 {
            pos_at_60 = pos;
        }
    }
    let (y, y_at_60, x) = (pos.y, pos_at_60.y, pos.x);
    let v = vel.y;
    println!(
        "短窗(60 tick) y={y_at_60:.4} | 长窗(1800 tick) 末 y={y:.4} v={v:+.3} | 横向 x={x:+.4}"
    );
    assert!(
        y_at_60 > 0.70,
        "短窗该被**接住**（阈值 0.70；3 维实测 0.9968）——短窗红了说明耦合本身坏了"
    );
    assert!(
        y < -1.0,
        "长窗**仍会逃逸**（3 维实测末 y={y:.4}、v={v:+.3}、横向 x={x:+.1}）⇒ **缺口是开着的**。\
         逃逸方式是**横向滑出**（x 跑到 +141 m）——1 维自扮引擎看不见这条路径，\
         所以它当年托住是**仪器假象**（§8.4.16）；修好后把这条翻成 y > 0.5"
    );
}
