// ligint - interazioni non covalenti ligando-proteina da un file PDB
//
// API pubblica della libreria. Il flusso tipico:
//
//     ligint::Options opt;            // soglie, ligando da analizzare, template SDF/MOL2...
//     opt.ligName = "ATP";
//     ligint::Result r = ligint::analyzeFile("1ATP.pdb", opt);
//     for (const auto& lig : r.ligands)
//         for (const auto& h : lig.hits) { /* h.type, h.ligandAtom, h.proteinAtom, h.distance ... */ }
//
// Senza ligName si ottiene solo l'elenco dei ligandi presenti (Result::listingOnly).
// La libreria non stampa nulla: avvisi e informazioni sono restituiti in Result.
#pragma once

#include <array>
#include <climits>
#include <iosfwd>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef LIGINT_VERSION
#define LIGINT_VERSION "0.1.0"
#endif

namespace ligint {

using V3 = std::array<double, 3>;

// ---------------------------------------------------------------------------------------------
// Soglie
// ---------------------------------------------------------------------------------------------
struct Params {
    static constexpr double kMaxGridRadius = 6.0;   // massimo raggio di ricerca supportato (Angstrom)

    double hbMin = 2.2, hbMax = 3.5, hbHA = 2.7, hbDHA = 120.0;          // legami a idrogeno
    double hydMin = 3.3, hydMax = 4.0, saltMax = 4.0, repMax = 4.5;       // idrofobici, ponti salini, repulsioni
    double halMax = 3.5, halAngle = 140.0, metalMax = 3.5;                // alogeni, metalli (tetto globale)
    double piMax = 5.5, piOffset = 3.5, piAngle = 30.0, cpMax = 6.0, clash = 0.5;
    double wbMin = 75.0, wbMax = 140.0;                                   // angolo ligando-O(acqua)-proteina
    double hbLong = 0.0, hbLongTol = 40.0;                                // H-bond lunghi (0 = spento)
    double saltCentroid = 0.0;                                            // ponti salini per centroidi (0 = spento)
    double hbAnte = 90.0;                                                 // angolo antecedente-atomo...partner

    // Imposta una soglia per nome (come nella riga di comando, senza "--"); false se il nome e' sconosciuto.
    bool set(const std::string& name, double value);
    // Preset per il confronto con PLIP: hbLong 4.1, saltCentroid 5.5, halMax 4.0.
    void applyPlipLike();
    // Limita a kMaxGridRadius le soglie usate dalla ricerca a griglia; restituisce gli avvisi.
    std::vector<std::string> clampToGrid();
};

// ---------------------------------------------------------------------------------------------
// Template chimici (SDF / MOL2)
// ---------------------------------------------------------------------------------------------
struct TAtom { std::string el, type; V3 p{0, 0, 0}; int fc = 0; };
struct TBond { int a, b, order; };          // indici 0-based; ordine 4 = aromatico
struct Tmpl { std::string name; std::vector<TAtom> at; std::vector<TBond> bd; };

// Leggono tutte le molecole del file. Ritorna false se il file non si apre; i problemi su singole
// molecole (V3000, file malformato) sono aggiunti a `warnings` e la molecola e' saltata.
bool loadSdf(const std::string& path, std::vector<Tmpl>& out, std::vector<std::string>& warnings);
bool loadMol2(const std::string& path, std::vector<Tmpl>& out, std::vector<std::string>& warnings);

// ---------------------------------------------------------------------------------------------
// Input
// ---------------------------------------------------------------------------------------------
struct Options {
    Params params;
    std::string ligName;            // nome del residuo HETATM (es. "ATP"); vuoto = solo elenco dei ligandi
    std::string chain;              // opzionale: catena del ligando
    int resi = INT_MIN;             // opzionale: numero di residuo (INT_MIN = qualsiasi)
    bool useWater = true;           // considera le acque (ponti d'acqua)
    bool metalDonors = false;       // O/N coordinati a un metallo restano possibili donatori
    bool backboneHyd = false;       // contatti apolari con CA/C della catena principale
    bool mergeLinked = true;        // unisce i residui HETATM legati covalentemente in un solo ligando
    std::vector<Tmpl> templates;    // annotazione chimica dei ligandi (abbinamento per coordinate)
};

// ---------------------------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------------------------
struct AtomRecord {
    std::string name, res, chain, element;
    char icode = ' ';
    int resi = 0;
    V3 pos{0, 0, 0};
    bool hetero = false;
};
struct BondRecord { int a, b, order; };     // indici nel vettore di atomi; order 0 = ignoto, 4 = aromatico

struct HitRecord {
    std::string type;               // es. "Legame a idrogeno", "Ponte salino", "pi-stacking (parallelo)"
    std::string ligandAtom;         // etichetta "RES nnnC:ATOMO" (o dell'anello/gruppo)
    std::string proteinAtom;
    std::string residue;            // residuo della proteina, es. "ASP 189A"
    std::string note;
    double distance = 0;            // Angstrom
    bool favorable = true;          // false per clash e repulsioni
    V3 ligandPos{0, 0, 0};          // punto del ligando (atomo, centroide dell'anello o del gruppo carico)
    V3 proteinPos{0, 0, 0};
    V3 waterPos{0, 0, 0};           // solo ponti d'acqua
    bool hasWater = false;
};

struct LigandResult {
    std::string id;                 // etichetta del ligando, es. "BZC 201A"
    std::string resName, chain;     // primo residuo
    int resi = 0;
    std::vector<std::string> residues;      // residui uniti nel ligando composto (se piu' di uno)
    std::string annotation;                 // come e' stata ottenuta la chimica (template, H del PDB, euristica)
    std::vector<std::string> info;          // righe informative (composto, annotazione, B-factor)
    std::vector<AtomRecord> atoms;          // atomi pesanti del ligando
    std::vector<BondRecord> bonds;
    std::size_t nRings = 0;
    std::vector<HitRecord> hits;            // ordinati per tipo, atomo, distanza
    std::vector<AtomRecord> contactAtoms;   // tutti gli atomi dei residui di proteina che interagiscono
    std::vector<BondRecord> contactBonds;
    std::vector<V3> waters;                 // acque coinvolte nei ponti d'acqua
    std::string metalSphere;                // testo: sfera di coordinazione dei metalli toccati
};

struct LigandId {                           // un residuo HETATM presente nel file
    std::string res, chain;
    int resi = 0;
    char icode = ' ';
    int atoms = 0;
};

struct Result {
    bool listingOnly = false;               // nessun ligando richiesto o trovato: solo elenco
    std::map<std::string, int> hetCounts;   // nome -> numero di residui (solo se listingOnly)
    std::vector<LigandId> hetResidues;      // idem, residuo per residuo
    std::string notFound;                   // nome richiesto ma assente (solo se listingOnly)
    std::vector<std::string> warnings;      // avvisi (coordinate non valide, soglie limitate...)
    std::vector<std::string> info;          // informazioni globali (H della proteina, template caricati)
    std::vector<LigandResult> ligands;
};

// Analizza il PDB. Lancia std::runtime_error se il file non si apre.
Result analyzeFile(const std::string& path, const Options& opt);
Result analyzeStream(std::istream& in, const Options& opt);

// ---------------------------------------------------------------------------------------------
// Rapporti testuali (identici all'output del programma da riga di comando)
// ---------------------------------------------------------------------------------------------
void printText(std::ostream& out, const Result& r);                       // tabella leggibile
void printCsv(std::ostream& out, std::ostream& info, const Result& r);    // CSV su `out`, informazioni su `info`

}  // namespace ligint
