# Validation and limits

## How to read the results

ligint applies **geometric criteria** to the atomic coordinates of a model and lists the contacts that satisfy them. Keep in mind:

- A listed interaction means that the geometry of the model fits the criterion, not that the interaction contributes to binding, and
  not that it exists in solution. A missing interaction can be a real absence, a limit of the model (low resolution, a poorly defined
  ligand: look at the B-factor note) or a consequence of an assumption about the chemistry.
- Hydrogens are normally absent from crystal structures, so donor and acceptor roles and charges are **inferred**. When a ligand is
  chemically unusual, or when the protonation state matters, provide an SDF or MOL2 file ([templates](interactions.md#sdf-and-mol2)).
- Water bridges depend on the waters that the structure contains, and many of them are not well defined.
- Unfavorable contacts (clashes, repulsions) are reported to help you judge a pose or a model; they are not errors of the program.
- The thresholds are conventions. Changing them changes the results; the **PLIP-like** preset is one such change, chosen to get closer to
  another widely used program.

## How the results were checked

ligint has been checked in several ways, none of which is a validation against experiment:

- **Automatic tests** on synthetic fragments with known geometry and known expected outcome, on real complexes whose interactions were
  checked by hand, and on malformed input.
- **Comparison with PLIP** on 150 complexes, reported in detail in [Comparison with PLIP](benchmark.md). The agreement is moderate
  and the disagreements are mostly conventions of the two programs, not errors of one of them.
- **Equivalence tests** between the formats: the same structure read from PDB and from mmCIF gives the same results.

<!-- include: README.md "Verification" -->

<!-- include: README.md "Known limits and differences from PLIP" -->

<!-- include: README.md "Test status" -->
