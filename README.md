# Wheat IBS Project Record

This folder keeps a local record of the wheat IBS identification project, including:

- workflow notes
- decision logic for sample filtering and retention
- code snapshots used during the project
- config template examples

## Structure

- `docs/project_log.md`
  Project summary, workflow, threshold logic, and sample decision rules.
- `docs/server_notes.md`
  Recommended server-side directory layout and execution notes.
- `scripts/`
  Code snapshots for the pipeline and report scripts.
- `configs/`
  Config template snapshot.

## Current Project Scope

The project focuses on:

- VCF filtering from chromosome-level VCF files
- IBS calculation with PLINK
- within-group repeatability checks
- between-group identity checks
- Z23-based decision logic for transcriptome samples

## Main Datasets

- Random-site resequencing VCFs for `Z23` and `B25`
- Exon / genotype / transcriptome VCFs for `Z23`, `B25`, `TC`, `SC`, `FC`

## Practical Analysis Principle

- Use random-site DNA VCF as the primary identity evidence
- Use exon / RNA-based VCF as supporting evidence
- Use `best match + second best + margin` rather than a single IBS value
- Use `Z23` as the first identity baseline for transcriptome sample decisions

