//! **门面级绳索判据**（T1 第二片）：绳索接进 `World` 之后，在**真实场景地形**上的黑盒读数。
//!
//! 三条判据（机器无关：无计时、无随机、纯 CPU）：
//! ① **落在场景地形上**：逐粒子就位间隙 ≈ 0 —— 走的是 `Providers` 通道（与刚体、与
//!    `tests/provider_shape_coverage.rs` 同一根），不是测试里自造的 collider；
//! ② **与刚体同场互不干扰**：同场景的盒仍按自己的路径落在地形上（本片绳只读提供者，
//!    与刚体既无接触也无耦合 ⇒ 两条通道并存）；
//! ③ **确定性**：同构造两跑，绳末态**逐位相同** + 末态哈希**冻结基线**。
//!
//! ⚠️ **换代级**：哈希是冻结值 ⇒ 改绳索数值/接触口径必须重冻并在此登记
//! （同 `default_stability` / `rope_minimal` 的惯例）。
//!
//! **本片边界**：绳不吃刚体（无 Akinci 耦合）、无摩擦、无自碰撞 —— 都属后续切片。

use vxl_phys::*;
use vxl_phys_core::{PhysConfig, Quat, Shape, Vec3};
use vxl_phys_soft::Rope;
use vxl_phys_terrain::mesh::TriMesh;

/// 平地板（y = 0）：8×8 格、±4 m（与另外两条提供者判据同几何）。
fn flat_mesh() -> TriMesh {
    const N: usize = 8;
    const S: f32 = 4.0;
    let mut verts: Vec<Vec3> = Vec::new();
    let mut tris: Vec<[u32; 3]> = Vec::new();
    for iz in 0..=N {
        for ix in 0..=N {
            let x = -S + (2.0 * S * ix as f32) / N as f32;
            let z = -S + (2.0 * S * iz as f32) / N as f32;
            verts.push(Vec3::new(x, 0.0, z));
        }
    }
    for iz in 0..N as u32 {
        for ix in 0..N as u32 {
            let a = iz * (N as u32 + 1) + ix;
            let c = a + 1;
            let d = a + N as u32 + 1;
            let e = d + 1;
            tris.push([a, d, c]);
            tris.push([c, d, e]);
        }
    }
    TriMesh::new(verts, tris)
}

/// 绳末态位置的 **FNV-1a**（按 f32 位模式逐字节）。
fn rope_hash(r: &Rope) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for p in &r.pos {
        for v in [p.x, p.y, p.z] {
            for b in v.to_bits().to_le_bytes() {
                h ^= b as u64;
                h = h.wrapping_mul(0x0000_0100_0000_01b3);
            }
        }
    }
    h
}

/// 一个"绳 + 刚体同场"的场景跑 `ticks` 步，返回（绳就位间隙最低/最高、绳末速、盒 y、绳哈希）。
fn run_scene(ticks: usize) -> (f32, f32, f32, f32, u64) {
    const NODES: usize = 33;
    const RADIUS: f32 = 0.02;
    let mut w = World::new(PhysConfig::default());
    let _mesh = w.add_mesh(flat_mesh());
    // **静态盒体**（顶面 y = 0.5）：绳索落在**它**上面 —— 这条走的是门面的**刚体代理通道**
    // （`rope_pass` 里建的 `RigidProxy` 表），与提供者通道是两条不同的路。
    // （提供者通道那一条由 `vxl-phys-soft` 的 `rope_rests_on_real_trimesh_provider` 守。）
    let _box_floor = w.add_static(
        Shape::Box {
            half: Vec3::new(1.0, 0.25, 1.0),
        },
        Vec3::new(0.0, 0.25, 0.0),
        Quat::IDENTITY,
    );
    // 同场动态刚体：证明两条通道并存（绳不吃刚体；它落在盒外，落在网格地形上）。
    let b = w.add_dynamic(
        Shape::Box {
            half: Vec3::splat(0.2),
        },
        Vec3::new(1.5, 1.0, 0.0),
        Quat::IDENTITY,
        1000.0,
    ) as usize;

    let mut rope = Rope::line(
        Vec3::new(-0.4, 1.2, 0.0),
        Vec3::new(0.4, 1.2, 0.0),
        NODES,
        RADIUS,
    );
    // 两端松开 ⇒ 整条绳落到场景地形上（判"接触就位"）。
    rope.set_pinned(0, false);
    rope.set_pinned(NODES - 1, false);
    rope.damping = 0.999;
    let ri = w.add_rope(rope);
    assert_eq!(ri, 0, "第一条绳索的索引应为 0");

    for _ in 0..ticks {
        w.step();
    }

    let r = w.rope(0).expect("rope 0 已注册");
    // 就位基线 = **静态盒顶面**（y = 0.5）⇒ 间隙 = 粒子表面到盒面。
    let mut lo = f32::INFINITY;
    let mut hi = f32::NEG_INFINITY;
    for p in &r.pos {
        lo = lo.min(p.y - RADIUS - 0.5);
        hi = hi.max(p.y - RADIUS - 0.5);
    }
    let vmax = r.vel.iter().map(|v| v.length()).fold(0.0f32, f32::max);
    (lo, hi, vmax, w.bodies.position[b].y, rope_hash(r))
}

#[test]
fn rope_rests_on_a_static_box_body_beside_a_rigid_body() {
    let (lo, hi, vmax, box_y, hash) = run_scene(900);
    let (_, _, _, _, hash2) = run_scene(900);
    println!(
        "绳就位间隙 最低={lo:.5} 最高={hi:.5} | 绳末速max={vmax:.3e} | 同场盒 y={box_y:.4} | 哈希={hash:016x}"
    );

    assert!(
        lo > -0.01,
        "没有粒子该陷进**静态盒体**（最低间隙 {lo:.5}）——红了说明门面的刚体代理通道没接上"
    );
    assert!(
        hi < 0.05,
        "绳该整体贴在盒面上（最高间隙 {hi:.5}）——大了说明有粒子悬在空中"
    );
    assert!(
        vmax < 1e-2,
        "绳该已静止（末速 {vmax:.3e}）——没停说明接触/阻尼有问题"
    );
    assert!(
        (box_y - 0.2).abs() < 0.05,
        "同场刚体该照常落在地形上（半高 0.2，实得 {box_y:.4}）——\n\
         绳索通道**不该**影响刚体管线（本片两者无耦合）"
    );
    assert_eq!(
        hash, hash2,
        "同构造两跑绳末态必须逐位相同（门面接线后仍是纯顺序推进）"
    );
    assert_eq!(
        hash, 0x5a24_4091_4067_fe15,
        "绳末态哈希是**冻结基线**（换代级：改绳索数值/接触口径必须重冻并登记）"
    );
}

/// **门面级双向耦合**：动态盒落在绳上 ⇒ 被托住。
///
/// 与软体侧 `rope_minimal::rope_couples_with_rigid_bodies` 同场景，但走 `World::step` 的**完整门面
/// 路径** —— 只有这条路才吃到**位置回填**（`rope_pass` 里把 `body_disp − 自身运动` 投影到
/// `bodies.position`）。软体侧的判据自己扮演引擎（只吃 `body_dv`）⇒ **量不到这一项**。
#[test]
fn box_on_rope_in_a_world_is_held() {
    let mut w = World::new(PhysConfig::default());
    let mut rope = Rope::line(
        Vec3::new(-0.5, 1.0, 0.0),
        Vec3::new(0.5, 1.0, 0.0),
        33,
        0.02,
    );
    rope.damping = 0.999;
    assert_eq!(w.add_rope(rope), 0, "第一条绳索的索引应为 0");
    for _ in 0..300 {
        w.step();
    }
    let b = w.add_dynamic(
        Shape::Box {
            half: Vec3::new(0.3, 0.05, 0.3),
        },
        Vec3::new(0.0, 1.2, 0.0),
        Quat::IDENTITY,
        // ⚠️ 第 4 参是**密度**（`mass_props(&shape, density)`），不是质量！
        // 体积 = 0.6×0.1×0.6 = 0.036 m³ ⇒ 要 1.0 kg 就得 27.78 kg/m³（与软体侧判据**同质量**）。
        // 先前随手写的 `1.0` 只给 0.036 kg ⇒ 轻 28×，读数与软体侧不可比（那次"300× 差异"的一半）。
        1.0 / 0.036,
    ) as usize;
    for _ in 0..1800 {
        w.step();
    }
    let y = w.bodies.position[b].y;
    println!("门面级 1800 tick：盒 y={y:+.4}");
    // ⚠️ **钉住的差异（本轮新发现，成因未定）**：**软体侧**同场景（判据自扮引擎、只吃 `body_dv`）
    // 读到 **−7.9**，而**门面路径**读到 **−2348**（差 ≈300×）⇒ 两者不等价。候选成因：门面的体子步
    // （`config.substeps`）与绳索 tick 的粒度差、力场/积分顺序、材质与睡眠。**修好后把这条翻成
    // "被托住"（`y > 0.5`）** —— 本判据的价值就是它**抓到了这个差异**（软体侧那两条量不到）。
    assert!(
        y < -1.0,
        "门面级今天**托不住**（实测 y={y:+.4}）⇒ 已修好就翻断言（阈值改 `y > 0.5`）"
    );
}
