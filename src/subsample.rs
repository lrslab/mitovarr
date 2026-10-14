//! Bounded, order-independent sampling of whole named BAM fragments.
use crate::{cli::MappingArgs, io::fastq};
use anyhow::{ensure, Result};
use serde_json::{json, Value};
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufWriter, Write},
};

type Key = (u64, Vec<u8>);

/// Seeded fragment-name priority shared by mapping and de novo subsampling.
/// Fixed FNV mixing plus avalanche: selection priority, not an integrity hash.
pub(crate) fn priority(name: &[u8], seed: u64) -> u64 {
    let mut x = 0xcbf29ce484222325u64 ^ seed;
    for b in name {
        x = (x ^ u64::from(*b)).wrapping_mul(0x100000001b3);
    }
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
    x ^ (x >> 31)
}
pub(crate) struct Selection {
    entries: BTreeMap<Key, usize>,
    cutoff: Option<Key>,
    limit: usize,
    seed: u64,
    pub seen: usize,
    pub retained: usize,
}

impl Selection {
    pub fn new(limit: usize, seed: u64) -> Self {
        Self {
            entries: BTreeMap::new(),
            cutoff: None,
            limit,
            seed,
            seen: 0,
            retained: 0,
        }
    }
    fn key(&self, name: &[u8]) -> Key {
        (priority(name, self.seed), name.to_vec())
    }
    pub fn observe(&mut self, name: &[u8]) {
        self.seen += 1;
        let key = self.key(name);
        if self.cutoff.as_ref().is_some_and(|c| &key >= c) {
            return;
        }
        *self.entries.entry(key).or_default() += 1;
        self.retained += 1;
        while self.retained > self.limit {
            let (key, count) = self.entries.pop_last().expect("nonempty selection");
            self.retained -= count;
            self.cutoff = Some(key);
        }
    }
    pub fn contains(&self, name: &[u8]) -> bool {
        self.entries.contains_key(&self.key(name))
    }
    pub fn report(&self) -> Value {
        json!({"enabled":true,"applied":self.seen>self.retained,"method":"whole-query-name priority prefix under exact target-BAM record cap", "seed":self.seed,"input_target_records":self.seen,"retained_target_records":self.retained,"retained_templates":self.entries.len(),"max_records":self.limit,"evidence_scope":"subsampled target-BAM fragments and their matching complete FASTQ pairs","limitation":"coverage and variant fractions describe retained evidence; subsampling can lose rare alleles or weak structural support"})
    }
    pub fn paired_reads(&self, args: &MappingArgs) -> Result<(MappingArgs, Value)> {
        let mut sampled = args.clone();
        let mut info = self.report();
        if args.local_assembly {
            let r1 = args.out.join("subsampled.R1.fastq");
            let r2 = args.out.join("subsampled.R2.fastq");
            let mut a = BufWriter::new(File::create(&r1)?);
            let mut b = BufWriter::new(File::create(&r2)?);
            let mut retained_pairs = 0usize;
            let scan = fastq::scan_pairs(
                args.r1.as_ref().unwrap(),
                args.r2.as_ref().unwrap(),
                |_, one, two| {
                    if self.contains(one.name(1)?.as_bytes()) {
                        one.write(&mut a)?;
                        two.write(&mut b)?;
                        retained_pairs += 1;
                    }
                    Ok(())
                },
            )?;
            a.flush()?;
            b.flush()?;
            ensure!(
                retained_pairs > 0,
                "subsampling retained no matching complete FASTQ pairs"
            );
            info["input_fastq_pairs"] = json!(scan.pairs);
            info["retained_fastq_pairs"] = json!(retained_pairs);
            sampled.r1 = Some(r1);
            sampled.r2 = Some(r2);
        }
        Ok((sampled, info))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn keeps_entire_fragments_under_cap_independent_of_record_order() {
        let names = ["a", "b", "c", "d", "e", "a", "b", "b", "c", "d", "e", "e"];
        let mut a = Selection::new(5, 42);
        let mut b = Selection::new(5, 42);
        for n in names {
            a.observe(n.as_bytes());
        }
        for n in names.iter().rev() {
            b.observe(n.as_bytes());
        }
        assert_eq!(a.entries, b.entries);
        assert!(a.retained <= 5 && a.retained > 0);
        for (key, count) in &a.entries {
            assert_eq!(
                *count,
                names.iter().filter(|n| n.as_bytes() == key.1).count()
            );
        }
    }
    #[test]
    fn below_cap_keeps_every_record_and_seed_changes_selection() {
        let mut all = Selection::new(100, 42);
        let mut a = Selection::new(10, 42);
        let mut b = Selection::new(10, 43);
        for i in 0..50 {
            for _ in 0..2 {
                let n = format!("read{i}");
                all.observe(n.as_bytes());
                a.observe(n.as_bytes());
                b.observe(n.as_bytes());
            }
        }
        assert_eq!(all.retained, 100);
        assert_eq!(a.retained, 10);
        assert_ne!(a.entries, b.entries);
    }
}
