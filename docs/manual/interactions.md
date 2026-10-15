# Interactions and assumptions

This page says what ligint recognizes and with which criteria. The criteria are geometric: they are applied to the coordinates of the
model, together with a few chemical assumptions that the next sections list. Every threshold can be changed (see the
[command-line reference](cli.md#thresholds)).

<!-- include: README.md "What it recognizes (default thresholds)" -->

## Interaction type names

The `type` column of the CSV, the **Type** column of the tables and the library constants use these names:
`Hydrogen bond`, `Salt bridge`, `Electrostatic repulsion`, `Hydrophobic contact`, `Apolar contact (backbone)`, `pi-stacking (parallel)`,
`pi-stacking (T-shaped)`, `Cation-pi`, `Halogen bond`, `Metal coordination`, `Water bridge`, `Steric clash` and `Possible covalent bond`.
Steric clashes and electrostatic repulsions are **unfavorable** contacts, listed with `-` in the favorable column. A possible covalent bond
(distance below 2.0 A) is only a **flag** that the contact is too close to be a non-covalent one (it is listed with `+`). None of the three
enters the [fingerprint](fingerprint.md).
