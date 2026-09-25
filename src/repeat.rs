use std::collections::BTreeSet;

use serde::Serialize;

use crate::pileup::{dominant_base, histogram_median, PileupPosition};

#[derive(Debug, Clone, Copy)]
pub struct RepeatSettings {
    pub enabled: bool,
    pub window: usize,
    pub ambiguity_min_triggers: usize,
    pub repeat_min_mapq: u8,
}

#[derive(Debug, Clone, Serialize)]
pub struct WindowMetrics {
    pub start: usize,
    pub end: usize,
    pub median_mapq: f64,
    pub sec_supp_frac: f64,
    pub depth_ratio: f64,
    pub mean_entropy: f64,
    pub triggers: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct AmbiguityInterval {
    pub start: usize,
    pub end: usize,
    pub len: usize,
    pub triggers: String,
    pub median_mapq: f64,
    pub sec_supp_frac: f64,
    pub depth_ratio: f64,
    pub mean_entropy: f64,
    pub action: String,
}

pub fn detect_ambiguity(
    pileup: &[PileupPosition],
    settings: RepeatSettings,
) -> Vec<AmbiguityInterval> {
    if !settings.enabled || pileup.is_empty() {
        return Vec::new();
    }

    let window = settings.window.max(1);
    let step = (window / 2).max(1);

    let baseline_depth = median_depth(pileup);

    let mut flagged_windows = Vec::new();

    let mut start = 0usize;
    while start < pileup.len() {
        let end = (start + window).min(pileup.len()) - 1;
        let span = &pileup[start..=end];

        let mut hist = [0u32; 61];
        let mut sec_supp_sum = 0u64;
        let mut total_obs_sum = 0u64;
        let mut depth_sum = 0f64;
        let mut entropy_sum = 0f64;
        let mut branch_sites = 0usize;

        for pos in span {
            for (idx, count) in pos.mapq_hist.iter().copied().enumerate() {
                hist[idx] = hist[idx].saturating_add(count);
            }
            sec_supp_sum = sec_supp_sum.saturating_add(u64::from(pos.sec_supp_obs));
            total_obs_sum = total_obs_sum.saturating_add(u64::from(pos.total_obs));
            depth_sum += f64::from(pos.depth);
            entropy_sum += pos.entropy();

            if is_graph_branch_site(pos) {
                branch_sites = branch_sites.saturating_add(1);
            }
        }

        let span_len = span.len() as f64;
        let mean_depth = if span_len > 0.0 {
            depth_sum / span_len
        } else {
            0.0
        };
        let depth_ratio = if baseline_depth > 0.0 {
            mean_depth / baseline_depth
        } else {
            0.0
        };

        let sec_supp_frac = if total_obs_sum > 0 {
            sec_supp_sum as f64 / total_obs_sum as f64
        } else {
            0.0
        };

        let mean_entropy = if span_len > 0.0 {
            entropy_sum / span_len
        } else {
            0.0
        };

        let median_mapq = histogram_median(&hist);

        let mut triggers = Vec::new();

        if median_mapq < f64::from(settings.repeat_min_mapq) {
            triggers.push("low_mapq".to_string());
        }
        if sec_supp_frac > 0.15 {
            triggers.push("sec_supp_high".to_string());
        }
        if !(0.55..=1.8).contains(&depth_ratio) {
            triggers.push("depth_anomaly".to_string());
        }
        if mean_entropy > 1.0 {
            triggers.push("high_entropy".to_string());
        }

        let min_branch_sites = (window / 50).max(1);
        let branch_frac = if span_len > 0.0 {
            branch_sites as f64 / span_len
        } else {
            0.0
        };
        if branch_sites >= min_branch_sites && branch_frac >= 0.05 {
            triggers.push("graph_branch".to_string());
        }

        if triggers.len() >= settings.ambiguity_min_triggers {
            flagged_windows.push(WindowMetrics {
                start,
                end,
                median_mapq,
                sec_supp_frac,
                depth_ratio,
                mean_entropy,
                triggers,
            });
        }

        if end == pileup.len() - 1 {
            break;
        }
        start = start.saturating_add(step);
    }

    merge_flagged_windows(&flagged_windows, window)
}

fn is_graph_branch_site(pos: &PileupPosition) -> bool {
    if pos.depth < 10 {
        return false;
    }

    let total_weight = pos.total_base_weight();
    let (major_idx, major_weight) = dominant_base(&pos.base_weights);
    let minor_weight = second_best_weight(&pos.base_weights, major_idx);

    let major_frac = if total_weight > 0.0 {
        major_weight / total_weight
    } else {
        0.0
    };
    let minor_frac = if total_weight > 0.0 {
        minor_weight / total_weight
    } else {
        0.0
    };

    let base_branch = major_frac <= 0.80 && minor_frac >= 0.20;
    if base_branch {
        return true;
    }

    let mut best = 0.0f64;
    let mut second = 0.0f64;
    for stats in pos.insertions.values().chain(pos.deletions.values()) {
        let w = stats.weight;
        if w > best {
            second = best;
            best = w;
        } else if w > second {
            second = w;
        }
    }

    best > 0.0 && second > 0.0 && second / best >= 0.80
}

fn second_best_weight(weights: &[f64; 5], major_idx: usize) -> f64 {
    weights
        .iter()
        .copied()
        .enumerate()
        .filter(|(idx, _)| *idx != major_idx)
        .map(|(_, w)| w)
        .fold(0.0, f64::max)
}

pub fn mask_from_intervals(reference_len: usize, intervals: &[AmbiguityInterval]) -> Vec<bool> {
    let mut mask = vec![false; reference_len];
    for interval in intervals {
        if interval.start >= reference_len {
            continue;
        }
        let end = interval.end.min(reference_len.saturating_sub(1));
        for bit in mask.iter_mut().take(end + 1).skip(interval.start) {
            *bit = true;
        }
    }
    mask
}

fn merge_flagged_windows(flagged: &[WindowMetrics], window: usize) -> Vec<AmbiguityInterval> {
    if flagged.is_empty() {
        return Vec::new();
    }

    let max_gap = window / 2;
    let mut merged = Vec::new();

    let mut cur_start = flagged[0].start;
    let mut cur_end = flagged[0].end;
    let mut triggers = BTreeSet::new();
    let mut metric_count = 0usize;
    let mut median_mapq_sum = 0.0;
    let mut sec_supp_frac_sum = 0.0;
    let mut depth_ratio_sum = 0.0;
    let mut entropy_sum = 0.0;

    let flush_current = |merged: &mut Vec<AmbiguityInterval>,
                         cur_start: usize,
                         cur_end: usize,
                         triggers: &BTreeSet<String>,
                         metric_count: usize,
                         median_mapq_sum: f64,
                         sec_supp_frac_sum: f64,
                         depth_ratio_sum: f64,
                         entropy_sum: f64| {
        let count = metric_count.max(1) as f64;
        merged.push(AmbiguityInterval {
            start: cur_start,
            end: cur_end,
            len: cur_end.saturating_sub(cur_start) + 1,
            triggers: triggers.iter().cloned().collect::<Vec<_>>().join("|"),
            median_mapq: median_mapq_sum / count,
            sec_supp_frac: sec_supp_frac_sum / count,
            depth_ratio: depth_ratio_sum / count,
            mean_entropy: entropy_sum / count,
            action: "strict_filtered".to_string(),
        });
    };

    for win in flagged {
        if win.start <= cur_end.saturating_add(max_gap) {
            cur_end = cur_end.max(win.end);
        } else {
            flush_current(
                &mut merged,
                cur_start,
                cur_end,
                &triggers,
                metric_count,
                median_mapq_sum,
                sec_supp_frac_sum,
                depth_ratio_sum,
                entropy_sum,
            );

            cur_start = win.start;
            cur_end = win.end;
            triggers.clear();
            metric_count = 0;
            median_mapq_sum = 0.0;
            sec_supp_frac_sum = 0.0;
            depth_ratio_sum = 0.0;
            entropy_sum = 0.0;
        }

        for t in &win.triggers {
            triggers.insert(t.clone());
        }

        metric_count = metric_count.saturating_add(1);
        median_mapq_sum += win.median_mapq;
        sec_supp_frac_sum += win.sec_supp_frac;
        depth_ratio_sum += win.depth_ratio;
        entropy_sum += win.mean_entropy;
    }

    flush_current(
        &mut merged,
        cur_start,
        cur_end,
        &triggers,
        metric_count,
        median_mapq_sum,
        sec_supp_frac_sum,
        depth_ratio_sum,
        entropy_sum,
    );

    merged
}

fn median_depth(pileup: &[PileupPosition]) -> f64 {
    let mut depths = pileup.iter().map(|p| p.depth).collect::<Vec<_>>();
    if depths.is_empty() {
        return 0.0;
    }

    depths.sort_unstable();

    let mid = depths.len() / 2;
    if depths.len() % 2 == 1 {
        f64::from(depths[mid])
    } else {
        (f64::from(depths[mid - 1]) + f64::from(depths[mid])) / 2.0
    }
}

#[cfg(test)]
mod tests {
    use crate::pileup::{init_pileup, PileupPosition};

    use super::*;

    #[test]
    fn test_mask_from_intervals() {
        let intervals = vec![AmbiguityInterval {
            start: 2,
            end: 4,
            len: 3,
            triggers: "low_mapq".into(),
            median_mapq: 10.0,
            sec_supp_frac: 0.2,
            depth_ratio: 2.0,
            mean_entropy: 1.2,
            action: "strict_filtered".into(),
        }];

        let mask = mask_from_intervals(10, &intervals);
        assert!(!mask[1]);
        assert!(mask[2]);
        assert!(mask[4]);
        assert!(!mask[5]);
    }

    #[test]
    fn test_detect_ambiguity_triggers() {
        let mut pileup = init_pileup(b"AAAAAAAAAA");
        for pos in &mut pileup {
            pos.mapq_hist[5] = 10;
            pos.total_obs = 10;
            pos.sec_supp_obs = 3;
            pos.depth = 40;
            pos.base_weights = [10.0, 10.0, 5.0, 5.0, 0.0];
        }

        let intervals = detect_ambiguity(
            &pileup,
            RepeatSettings {
                enabled: true,
                window: 4,
                ambiguity_min_triggers: 2,
                repeat_min_mapq: 30,
            },
        );

        assert!(!intervals.is_empty());
    }

    #[test]
    fn test_median_depth() {
        let mut pileup = vec![PileupPosition::new(b'A'); 3];
        pileup[0].depth = 10;
        pileup[1].depth = 30;
        pileup[2].depth = 20;
        assert_eq!(median_depth(&pileup), 20.0);
    }
}
