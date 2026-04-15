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
3. If it matches `B25`, keep it for later resequencing
4. If it matches neither, mark it as a removal candidate

## Final Labels

- `KEEP_Z23`
  Transcriptome sample matches expected `Z23`
- `KEEP_B25`
  Transcriptome sample does not match expected `Z23`, but matches corresponding `B25`
- `REVIEW`
  Match is ambiguous, usually because margin is too small
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

