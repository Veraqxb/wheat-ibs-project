# Decision Rules

## Primary Principle

Use random-site DNA VCF as the main identity evidence.

Use exon / transcriptome VCF as supporting evidence.

## Pairwise IBS Variables

For each sample, derive:

- `ibs_expected`
- `ibs_best`
- `second_best`
- `margin = ibs_best - second_best`
- `status`

## Geno Threshold Selection

Candidate thresholds:

- `0.2`
- `0.3`
- `0.4`
- `0.5`

Selection rule:

1. Exclude thresholds with too few retained SNPs
2. Prefer the threshold with the highest number of correct expected matches
3. If tied, prefer the threshold with the larger mean margin
4. If still tied, prefer the stricter threshold

## Main DNA Judgement

Recommended logic:

- `MATCH`
  Expected pair is the best match and margin is clear
- `SUSPICIOUS`
  Expected pair is still best, but margin is small or IBS is lower than expected
- `MISMATCH`
  Expected pair is not the best match

## Transcriptome Judgement

Use `Z23` as the first reference.

For each transcriptome sample:

1. Check whether it matches the expected `Z23` sample
2. If not, check whether it matches the corresponding `B25` sample
3. Check the same-row `2group` reference context before accepting a B25 rescue
4. If the expected `B25` reference is trusted or review-only, use B25 as rescue evidence
5. If the expected `B25` reference is `True_mismatch` or `No_data`, do not silently rescue; write `REVIEW_B25_CONTEXT_RISK`
6. If it matches neither Z23 nor usable B25/cluster neighbors, mark it as a removal candidate

## 2group Reference Context

The `2group` DNA reference has two evidence layers:

- `cluster_table`
  Strict within-reference high-similarity neighborhoods built with `IBS >= 0.99`.
- `reference_context`
  Row-wise bidirectional Z23/B25 status from B25-to-Z23 and Z23-to-B25 scans.

`reference_context` classes:

- `Exact_match`
  Both B25-to-Z23 and Z23-to-B25 scans return the expected map pair.
- `B25_only_exact`
  B25-to-Z23 returns the expected Z23, but the reverse direction is shifted.
- `Z23_only_exact`
  Z23-to-B25 returns the expected B25, but the reverse direction is shifted.
- `Expected_high_but_best_shifted`
  Expected pair has high IBS, but at least one direction prefers another high-IBS neighbor.
- `low_ibs_match`
  No high-confidence swapped identity is found, but expected IBS is below the strict DNA threshold.
- `True_mismatch`
  Expected pair is low while another candidate passes the DNA threshold.
- `No_data`
  Required pairwise IBS evidence is unavailable.

Reference confidence used by RNA rescue:

- `trusted`: `Exact_match`
- `review`: `B25_only_exact`, `Z23_only_exact`, `Expected_high_but_best_shifted`, `low_ibs_match`
- `exclude`: `True_mismatch`, `No_data`

## Final Labels

- `KEEP_Z23`
  Transcriptome sample matches expected `Z23`
- `KEEP_B25`
  Transcriptome sample does not match expected `Z23`, but matches corresponding usable `B25`
- `REVIEW`
  Match is ambiguous, usually because margin is too small
- `REVIEW_B25_CONTEXT_RISK`
  RNA matches B25, but the 2group reference context says the corresponding B25 is mismatched or missing
- `REMOVE`
  Sample matches neither expected `Z23` nor corresponding `B25`

## Removal Guidance

Recommended removal candidates:

- not matched to expected `Z23`
- not matched to corresponding `B25`
- weak best-vs-second-best separation
- low-confidence identity pattern across both datasets

## Retention Guidance

Recommended retention candidates:

- matched to expected `Z23`
- or not matched to `Z23` but clearly matched to corresponding `B25`
