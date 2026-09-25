# mitovarr development tasks

## Release preparation

- [ ] Choose a license and add its declaration and text before public release.
- [ ] Add the actual repository URL after the hosting location is chosen.
- [ ] Create the initial commit and prepare the release when ready to publish.

## Method development

- [x] Require linked unique-flank evidence for local contractions, with a masked-repeat negative control and a linked positive control.

- [x] Gate CDS terminal-stop interpretation on supported 3′ boundaries, with clipped T/TA and retained stop/gap regression tests.

- [ ] Extend sample-derived anchors and read threading to explore repeat exits more effectively. Preserve unresolved repeat tails and require support across both junctions; do not use evaluation-reference sequence to generate or select paths.
- [ ] Pool sample-derived connections or complete paths across main-k runs for a shared-fragment comparison. Do not select by gene completeness or uncalibrated scores from separate runs.
- [ ] Improve RNA terminal-mismatch handling using existing sample sequence and homology. Retain genuine missing-end cases and distinguish length variants from undetected models.
- [ ] Compare homologous NCRs for insertion, deletion, duplication, inversion and unresolved repeat structure, retaining locus-specific read support.

Use real reads for method evaluation and small portable fixtures for regression
coverage. Keep experiment protocols, result matrices and research logs outside
this repository; document implemented behavior and limitations in the method
guides.
