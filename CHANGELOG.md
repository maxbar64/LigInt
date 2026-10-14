# Changelog

## 0.1.0

Prima versione organizzata come progetto CMake.

- Libreria `ligint_core` con API pubblica (`analyzeFile`, `analyzeStream`, `printText`, `printCsv`, `loadSdf`, `loadMol2`);
  il risultato include le coordinate delle interazioni, gli atomi del ligando e dei residui a contatto per la grafica.
- Programma a riga di comando `ligint`: output identico a quello della versione a file singolo (verificato su 1075
  esecuzioni di 119 complessi in quattro modalità), più `--version`.
- GUI Qt `ligint-gui` con vista 3D, tabelle, filtri per tipo, esportazione CSV e immagine, template SDF/MOL2.
- Interazioni: legami a idrogeno (anche lunghi con geometria possibile), ponti salini (anche per centroidi), idrofobici,
  π-stacking, cation-π, alogeni, metalli, ponti d'acqua, clash; ligandi di più residui uniti; DNA/RNA come partner.
- Test: C++ (`core_unit`), Python (66 test), test di fumo della GUI offscreen, confronto e benchmark contro PLIP.
