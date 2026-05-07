# BCMA scFv combination-design workflow

## Input files

- `relaxed_teclistamab_bcma.pdb`: chain A is the scFv template and chain B is BCMA.
- `bcma_experiment.csv`: partial wet-lab results used as design priors; higher `property_value` means stronger binding.
- `design_bcma_combinations.py`: reproducible generator for the 20 constrained combination designs.

## Residue-region assumptions

- VH/heavy region: chain A residues 1-121.
- VL/light region: chain A residues 137-244.
- Candidate mutations are limited to the 17 user-provided sites.
- Heavy-chain designs use exactly 6 mutations, split as 2 HCDR1 + 2 HCDR2 + 2 HCDR3.
- Light-chain designs use 5 mutations, with at least 1 mutation from each LCDR bucket.

## Output files

- `bcma_combination_designs.csv`: the 20 requested combination-mutation forms plus validation, rationale, and blank columns reserved for SWISS-MODEL/PISA outputs.
- `bcma_mutant_sequences.fasta`: chain A amino-acid sequences for the same 20 mutants.
- `bcma_design_run.log`: generation log for the CSV/FASTA result files, stored in this `summary` folder together with the results.
- `bcma_validation.log`: validation log confirming the generated result files exist in this `summary` folder and the 20 designs pass the requested mutation-count constraints.

## Online SWISS-MODEL and EMBL-PISA status

The online modeling/evaluation fields are intentionally marked as pending in `bcma_combination_designs.csv` because this execution environment cannot complete custom-template uploads through the interactive SWISS-MODEL and EMBL-PISA web sessions on the user's behalf.  The generated FASTA, design table, and log files are stored directly in the `summary` folder and are ready for those online runs.

Recommended manual/online sequence for each design:

1. Submit the corresponding sequence from `bcma_mutant_sequences.fasta` to SWISS-MODEL.
2. Use `relaxed_teclistamab_bcma.pdb` as the custom template, preserving chain A-chain B complex context where the web interface permits.
3. Download each modeled mutant PDB.
4. Upload each modeled mutant PDB to EMBL-PISA.
5. Record the chain A-chain B `interface area` and `ΔiG` values back into the `interface_area_A2` and `delta_iG_kcal_per_mol` columns in `bcma_combination_designs.csv`.

## Design-prior notes

- The experimentally best provided row was used as a strong prior for combining `S31D/P63A/K66R/V104T` with light-chain mutations such as `G159V/D188M/H232T`.
- `K166Q` was preferred over `K166R` where possible because the available `K166R`-containing historical combinations were frequently detrimental in the supplied experiment table.
- Adjacent sites were avoided when feasible; rows with unavoidable or intentional near-adjacent choices are flagged in `near_adjacent_pairs_pos_delta_le_2`.
