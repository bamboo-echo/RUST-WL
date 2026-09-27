//! Cloth world channel: fixed tick, stable registration order, static provider contact.

use super::*;

#[derive(Default)]
pub(crate) struct SoftDomain {
    pub ropes: Vec<vxl_phys_soft::Rope>,
    pub cloths: Vec<vxl_phys_soft::ClothSheet>,
}

impl World {
    /// Register a cloth sheet. Indices remain stable for the lifetime of this world.
    pub fn add_cloth(&mut self, cloth: vxl_phys_soft::ClothSheet) -> usize {
        self.soft.cloths.push(cloth);
        self.soft.cloths.len() - 1
    }

    pub fn cloths(&self) -> &[vxl_phys_soft::ClothSheet] {
        &self.soft.cloths
    }

    pub fn cloth(&self, index: usize) -> Option<&vxl_phys_soft::ClothSheet> {
        self.soft.cloths.get(index)
    }

    pub(crate) fn cloth_pass(&mut self) {
        if self.soft.cloths.is_empty() {
            return;
        }
        let dt = self.config.dt;
        let gravity = self.config.gravity;
        let count = self.providers.len() as u32;
        for cloth in &mut self.soft.cloths {
            cloth.step(dt, gravity, &self.providers, count);
        }
    }
}
