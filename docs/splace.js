// ========================================================================
// State
// ========================================================================
const state = {
    records: [],       // Array of { accession, organism, features, sequence, source }
    selectedGenes: new Set(),
    selectedFeatureTypes: new Set(["CDS"]),
    detectedDataType: "mt", // auto-detected: "mt" or "cp"
    pendingRemoveIndex: null,
    pendingEditIndex: null,
    headerTemplate: ["species", "accession"],
    hiddenColumns: new Set(["source", "accession", "kingdom", "phylum"]),
};

const MT_DEFAULT_GENES = ["COI", "COII", "COIII", "CYTB", "ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "ATP6", "ATP8"];
const CP_DEFAULT_GENES = ["rbcL", "matK", "ndhF", "atpB", "psaA", "psbA", "psbB", "psbC", "psbD", "psbE", "psbF", "psbH", "psbI", "psbJ", "psbK", "psbL", "psbM", "psbN", "psbT"];

const EXAMPLE_ACCESSIONS = "NC_061537, NC_023122, DQ080041";

const FEATURE_TYPE_INFO = {
    CDS: "Protein-coding sequences (genes translated into proteins)",
    rRNA: "Ribosomal RNA genes (12S, 16S, etc.)",
    tRNA: "Transfer RNA genes",
};

const FEATURE_TYPE_LABELS = {
    CDS: "PCGs",
    rRNA: "rRNAs",
    tRNA: "tRNAs",
};

// ========================================================================
// Gene Name Standardization Dictionaries
// ========================================================================
const MitochondrialGenes = {
    "12S": ["small subunit ribosomal RNA", "s-rRNA", "12S ribosormal RNA", "small ribosomal RNA subunit RNA", "12SrRNA", "12 ribosomal RNA", "rrnS", "12S ribosomal RNA subunit", "12S", "small ribosomal RNA", "small subunit ribosormal RNA", "12 rRNA", "12 S ribosomal RNA", "12S small subunit ribosomal RNA", "trnS", "Product small subunit ribosomal RNA", "12S-rRNA", "rRNA-12S", "12S ribosonal RNA", "12Srrn", "12S ribosome RNA", "12S ribsomal RNA", "12S rRNA", "12S ribosomal RNA", "12S ribosomal ribonucleic acid"],
    "16S": ["large subunit ribosomal RNA", "l-rRNA", "16S bibosomal RNA", "large ribosomal RNA subunit RNA", "rrnL", "16S ribosomal RNA subunit", "16S rivbosomal RNA", "l6S ribosomal RNA", "16S", "16S ribosamal RNA", "large ribosomal RNA", "16S rRNA", "16 S ribosomal RNA", "16S large subunit ribosomal RNA", "l-RNA", "16S-rRNA", "16Srrn", "16S ribosommal RNA", "16S ribosomal ribonucleic acid", "16S ribosomal RNA", "16S ribosome RNA", "16S recombinant RNA", "16S ribosomal RNA gene", "16S recombinant ribonucleic acid"],
    "ATP6": ["ATPase F0 subunit 6", "ATP synthase F0 subunit 6", "ATP synthase subunit 6", "ATPase 6", "ATP6", "ATP synthase FO subunit 6", "ATP synthase protein 6", "ATPase subunits 6", "ATP subunit 6", "ATPase subunit-6", "ATP synthase subunit F0 6", "adenosine triphosphatase subunit 6", "ATP synthase subunit-6", "ATP sythase F0 subunit 6", "ATPase 6 protein", "ATP synthase 6", "ATP synthetase F0 Subunit 6", "ATPase subuint 6", "ATPase sununit 6", "ATPase6 protein", "ATP Synthase Membrane Subunit 6", "TP synthase F0 subunit 6", "ATP sybthase F0 subunit 6", "ATP6ase", "ATP synthase A chain protein 6", "ATP synthetase subunit 6", "F0/F1 ATP synthase subunit 6", "disrupted ATP synthase F0 subunit 6", "ATPsynthase F0 subunit 6", "F1 ATPase subunit 6", "ATP sythase subunit 6", "adenine triphosphatase subunit 6", "F0-ATP synthase subunit 6", "F0-ATP synthase subunit6", "F0-ATPase subunit6", "ATP 6synthase 6", "adenosine triphosphate subunit 6", "ATPase subunit 6", "ATPase6", "adenosine triphosphate synthase-6"],
    "ATP8": ["ATPase F0 subunit 8", "ATP synthase F0 subunit 8", "ATP synthase subunit 8", "ATPase 8", "ATP8", "ATP synthase FO submit 8", "ATPase8", "ATP synthase protein 8", "ATPase subunits 8", "ATP subunit 8", "ATPase subunit-8", "ATP synthase subunit F0 8", "adenosine triphosphatase subunit 8", "ATP synthase subunit-8", "ATP sythase F0 subunit 8", "ATP synthase FO subunit 8", "ATPase 8 protein", "adenosine triphoshatase subunit 8", "ATP synthase 8", "ATPase subunit8", "ATP synthetase F0 Subunit 8", "adenosine triphosphate subunit 8", "ATPase8 protein", "ATP Synthase Membrane Subunit 8", "TP synthase F0 subunit 8", "product ATP synthase F0 subunit 8", "ATP sybthase F0 subunit 8", "ATP8ase", "ATP synthetase subunit 8", "F0/F1 ATP synthase subunit 8", "ATPsynthase F0 subunit 8", "F1 ATPase subunit 8", "ATP sythase subunit 8", "adenine triphosphatase subunit 8", "ATP syntahse F0 subunit 8", "F0-ATP synthase subunit 8", "F0-ATP synthase subunit8", "ATPase subunit 8", "F0-ATPase subunit8", "ATP synthase F0 subunit"],
    "COI": ["COX1", "cytochrome c oxidase subunit 1", "cytochrome oxidase gene", "coxidase subunit I", "COX", "c oxidase subunits I", "cytochrome-c-oxidase I", "Cytochrome Oxidase subunit I region", "c-oxidase subunit I", "c oxi- dase I", "c oxidase subunit-I", "cytochrome oxidase I region", "cytochrome oxydase I", "c oxydase I", "cytochrome-oxidase 1", "C Oxidase Gene Subunit-I", "C Oxidation I", "c oxidase I subunit", "cytochrome c oxidase I", "Cythocrome Oxidase I", "cytochrome oxidase I subunit", "Citochrome Oxidase I", "cytochromec oxidase I", "c oxidase submit I", "c oxidase unit I", "c oxidate subunit I", "cytochrome I", "cytochome oxidase subunit I", "cytochrome-c oxidase, subunit I", "cytochrome b oxidase subunit I", "cytochrome subunit I", "cytochrome-oxidase I", "COX-1", "cytochromoxidase I", "cytochrome oxidase 1", "cytochrome oxidase subunit 1", "cytochrome oxidase I", "C Oxidase type I", "cytochrome oxidase subunit I", "cytochrome oxidase I locus", "c oxidase subunit I sequences", "coxidase I", "c oxidase subunit I locus", "Cytochrome Oxidase unit I", "cytochrome oxidase subunits I", "cytochrome oxidase subunit I mtDNA", "cytochrome C oxidase subunit I", "Markers-Cytochrome Oxidase Subunit I", "C Oxidase, Subunit I", "chytochrome c oxidase subunit I", "I", "cytochrome oxidase subunit-1", "Cytochrome oxydase subunit 1", "cytochrome c oxidase subunit idase subunit I", "cytochrome c-oxidase subunit I", "cytochrome c oxidase subunits I", "cytchrome c oxidase subunit I", "subunit 1 of the cytochrome c oxidase", "cytochorome oxidase subunit I", "COI", "cytochrome c oxydase subunit 1", "cytochrome oxidase1", "COI protein", "cyt oxidase subunit 1", "cytochrome oxidase c subunit 1", "cytochrome oxidase c subunit I", "cytochrome oxydase subunit I", "cytochrome c oxidase polypeptide I", "cytochrome coxidase subunit I", "cytochrome c-oxidase subunit 1", "cytochrome c oxidase polypeptide 1", "CO I", "product: cytochrome c oxidase subunit I", "cytochome c oxidase subunit 1", "Cytochrome c oxidase subunit1", "cytochrome coxidase subunit 1", "cytochrome c oxidase subunit I"],
    "COII": ["COX-2", "c oxidase subunit II gene", "cytochrome oxidase II gene", "cytochrome oxidase subunit II gene", "cytochrome oxidase subunit II region", "coxidase subunit II", "c oxidase II gene", "cytochrome oxidase-II", "cytochrome c oxidase subunit 2", "chytochrome c oxidase subunit II", "II", "cytochrome c oxidase subunit II", "cytochrome oxidase subunit II", "cytochrome oxidase subunit 2", "cytohrome oxidase subunit II", "cytochrome oxidase subunit-2", "Cytochrome oxydase subunit 2", "cytochrome c oxidase subunit idase subunit II", "cytochrome c-oxidase subunit II", "cytochrome c oxidase subunits II", "cytochrome c oxidase II", "cytochrome oxidase II", "cytchrome c oxidase subunit II", "subunit 2 of the cytochrome c oxidase", "COX2", "cytochorom oxidase subunit II", "COII", "cytochrome c oxydase subunit 2", "CO2", "cytochrome oxidase subunit2", "COII protein", "cyt oxidase subunit 2", "cytochrome oxidase c subunit 2", "cytochrome oxidase c subunit II", "cytochrome oxydase subunit II", "cytochrome c oxidase polypeptide II", "cytochrome coxidase subunit II", "cytochome oxidase II", "cytochrome c-oxidase subunit 2", "cytochrome c oxidase polypeptide 2", "CO II", "cytochome c oxidase subunit 2", "Cytochrome c oxidase subunit2", "cytochrome coxidase subunit 2", "cytochome oxidase subunit 2", "C-terminal domain of Cytochrome c Oxidase subunit II"],
    "COIII": ["cytochrome oxidase subunit III gene", "c oxidase subunit III gene", "COX3", "COX 3", "COX-3", "c oxidase mitochondrial subunit III", "cytochrome c oxidase subunit 3", "cytochrome c oxidase subunit III", "chytochrome c oxidase subunit III", "cytochrome oxidase subunit III", "cytochrome oxidase subunit 3", "cytohrome oxidase subunit III", "cyctochrome c oxidase subunit III", "cutochrome oxidase subunit 3", "cytochrome oxidase subunit-3", "Cytochrome oxydase subunit 3", "cytochrome c oxidase subunit idase subunit III", "cytochrome c-oxidase subunit III", "cytochrome c oxidase subunits III", "cytochrome c oxidase III", "cytochrome oxidase III", "cytchrome c oxidase subunit III", "subunit 3 of the cytochrome c oxidase", "cytochorome oxidase subunit III", "COIII", "cytochrome c oxydase subunit 3", "CO3", "cytochrome oxidase subunit3", "cytochrome coxidase subunit III", "COIII protein", "cyt oxidase subunit 3", "cytochrome oxidase c subunit 3", "cytochrome oxidase c subunit III", "cytochrome oxydase subunit III", "cytochrome co oxidase subunit III", "cytochrome c oxidase subnunit III", "cytochrome oxidase sununit 3", "Cytochrome c oxidase polypeptide III", "cytochrome C oxidase asubunit 3", "cytochrome c oxidase sununit III", "cytochrome c-oxidase subunit 3", "cytochrome c oxidase polypeptide 3", "CO III", "cytochome c oxidase subunit 3", "cytochrome c oxidase subunit3", "cytochrome coxidase subunit 3"],
    "CYTB": ["Cytochrome b apoenzyme", "apoenzyme", "cytohrome b", "cytochome b", "cytochorome b", "cytchrome b", "cob", "cytochrome b protein", "Cythocrome b", "Cytb protein", "cytchorome b", "CYTB", "Cytochrome-b", "cytbochrome b", "apocytochrome b", "cytochromeb", "ctyb", "Cyt b", "apocytochome b", "cytochrome b", "cytochrome bc1"],
    "ND1": ["NADH dehydrogenase subunit 1", "NAD dehydrogenase subunit 1", "NADH dehydrogenase subunit-1", "NADH denydrogenase subunit 1", "NADH dehydrogenase 1", "NADH dehydrogenase subunits 1", "NADH dehydogenase subunit 1", "NADH dehydrogenase subunit1", "NADH dehydrogenase subunit #1", "Subunit 1 of the NADH ubiquinone oxidoreductase complex", "NADHsubunit 1", "NADH dehydrogenase subnit 1", "ND1", "NADH dehydrogenase subunit I", "NADH1", "NADH1 protein", "NADH subunit 1", "NADH-ubiquinone oxidoreductase chain 1", "NADH dehydrogenase, subunit 1", "NADH-ubiquinone oxidoreductase subunit I", "NADH 1", "NADH dehydrogenase subumit 1", "NADH-ubiquinone oxidoreductase subunit 1", "NADH ubiquinone oxidoreductase subunit 1", "nicotinamide adenine dinucleotide dehydrogenase subunit 1", "truncated NADH dehydrogenase subunit 1", "NADH-1", "NADH dehdrogenase subunit 1", "NaD1", "NADH dehydrogynase subunit 1"],
    "ND2": ["NADH dehydrogenase subunit 2", "NADH dehydrogenase subunit-2", "#NADH dehydrogenase subunit 2", "NADH denydrogenase subunit 2", "NADH dehydrogenase 2", "NADH dehydrogenase subunits 2", "NADH dehydrogenase subunit #2", "subunit 2 of the NADH ubiquinone oxidoreductase complex", "NADH deydrogenase subunit 2", "NADH dehydrogenase subnit 2", "ND2", "NADH dehydrogenase subunit II", "NADH2", "NADH2 protein", "NADH subunit 2", "NADH-ubiquinone oxidoreductase chain 2", "NADH dehdrogenase subnuit 2", "NADH-ubiquinone oxidoreductase subunit II", "NADH 2", "NADH dehydrogenase subumit 2", "NADH-ubiquinone oxidoreductase subunit 2", "NADH ubiquinone oxidoreductase subunit 2", "NADH dehydrognase subunit II", "nicotinamide adenine dinucleotide dehydrogenase subunit 2", "NADH dehydrogenase subunit2"],
    "ND3": ["NADH dehydrogenase subunit 3", "NAD dehydrogenase subunit 3", "NADH dehydrogenase subunit-3", "NADH denydrogenase subunit 3", "NADH dehydrogenase 3", "NADH dehydrogenase subunits 3", "subunit 3 of the NADH ubiquinone oxidoreductase complex", "NADH deydrogenase subunit 3", "NADH dehydrogenase subnit 3", "ND3", "NADH3", "NADH dehydrogenase subunit III", "NADH3 protein", "NADH subunit 3", "NADH-ubiquinone oxidoreductase chain 3", "ND3 NADH dehydrogenase subunit 3", "NADH dehydrogenasesubunit 3", "NADH dehydrogenase, subunit 3", "NADH-ubiquinone oxidoreductase subunit III", "truncated NADH dehydrogenase subunit 3", "NADH 3", "NADH dehydrogenase subumit 3", "NADH-ubiquinone oxidoreductase subunit 3", "NADH ubiquinone oxidoreductase subunit 3", "NADH dehydrogenase subunit3"],
    "ND4": ["NADH dehydrogenase subunit 4", "NAD dehydrogenase subunit 4", "NADH hehydrogenase subunit 4", "NADH dehrogenase subunit 4", "NADH dehydrogenase subunit-4", "NADH denydrogenase subunit 4", "NADH dehydrogenase 4", "NADH dehydrogenase subunits 4", "NADH dehydrosenase subunit 4", "subunit 4 of the NADH ubiquinone oxidoreductase complex", "NADH deydrogenase subunit 4", "NADH dehydrogenase subnit 4", "ND4", "NADH dehyodrogenase subunit 4", "NADH4", "NADH dehydrogenase subunit4", "NADH dehydrogenase subunit IV", "NADH4 protein", "NADH subunit 4", "NADH dehydrogenase sunbunit 4", "NADH-ubiquinone oxidoreductase chain 4", "NADH-ubiquinone oxidoreductase subunit IV", "NADH 4", "NADH dehydrogenase subumit 4", "NADH-ubiquinone oxidoreductase subunit 4", "NADH ubiquinone oxidoreductase subunit 4", "nicotinamide adenine dinucleotide dehydrogenase subunit 4", "truncated NADH dehydrogenase subunit 4"],
    "ND4L": ["NADH dehydrogenase subunit 4L", "NADH dehydrogenase subunit-4L", "NADH denydrogenase subunit 4L", "NADH dehydrogenase 4L", "ND4L", "NADH dehydrogenase subunits 4L", "subunit ND4L of the NADH ubiquinone oxidoreductase complex", "NADH4L protein", "NADH dehydrogenase subnit 4L", "NADH4L", "NADH dehydrogenase subunit 4 L", "NADH dehydrogenase subunit IV L", "NADH subunit 4L", "NADH-ubiquinone oxidoreductase chain 4L", "NADH-ubiquinone oxidoreductase subunit 4L", "NADH dehydrogenase, subunit 4L (complex I)", "NADH 4L", "NADH dehydrogenase subumit 4L", "HADH dehydrogenase 4L", "NADH dehydrogenase subujnit 4L", "NADH ubiquinone oxidoreductase subunit 4L", "nicotinamide adenine dinucleotide dehydrogenase subunit 4L", "NADH dehydrogenase subunit4L"],
    "ND5": ["NADH dehydrogenase subunit 5", "NAD dehydrogenase subunit 5", "NADH dehrogenase subunit 5", "NADH dehydrogenase subunit-5", "NADH hehydrogenase subunit 5", "NADH denydrogenase subunit 5", "NADH dehydrogenase 5", "ND5", "NADH dehydrogenase subunits 5", "NADH hydrogenase subunit 5", "subunit 5 of the NADH ubiquinone oxidoreductase complex", "NADH deydrogenase subunit 5", "NADH dehydrogenase subnit 5", "NADH dehydrodenase subunit 5", "NADH5", "HADA dehydrogenase subunit 5", "NADH dehydrogenase subunit V", "NADH5 protein", "NADH subunit 5", "NADH dehydrogenase subunit 5-0", "NADH-ubiquinone oxidoreductase chain 5", "NADH-ubiquinone oxidoreductase subunit V", "NADH dehydrogenase, subunit 5", "NADH dehydrogenase, subunit 5 complex I", "NADH 5", "NADH dehydrogenase subumit 5", "NADH-ubiquinone oxidoreductase subunit 5", "NADH ubiquinone oxidoreductase subunit 5", "nicotinamide adenine dinucleotide dehydrogenase subunit 5", "truncated NADH dehydrogenase subunit 5", "NADH dehydrogenase subunit5", "NADH dehydroghenase subunit 5"],
    "ND6": ["NADH dehydrogenase subunit 6", "NAD dehydrogenase subunit 6", "NADH dehydrogenase subunit-6", "NADH denydrogenase subunit 6", "NADH dehydrogenase 6", "ND6", "NADH dehydrogenase subunits 6", "subunit 6 of the NADH ubiquinone oxidoreductase complex", "NADH deydrogenase subunit 6", "NADH dehydrogenase subnit 6", "NADH6", "NADH dehydrogenase subunit VI", "NADH6 protein", "NADH subunit 6", "truncated NADH dehydrogenase subunit 6", "NADH dehygrogenase subunit 6", "NADH dehydrogenease subunit 6", "NADH-ubiquinone oxidoreductase chain 6", "NADH dehydrogenase, subunit 6", "NADH dsehydrogenase subunit 6", "NADH-ubiquinone oxidoreductase subunit VI", "NADH 6", "NADH dehydrogenase subumit 6", "NADH-ubiquinone oxidoreductase subunit 6", "NADH ubiquinone oxidoreductase subunit 6", "nicotinamide adenine dinucleotide dehydrogenase subunit 6", "NADH dehydrogenase subunit6"],
    "CR": ["Control region", "non-coding region", "putative control region", "control region 1", "control region ii", "control region i", "control region 2", "noncoding region", "pseudo control region", "cr", "control region (d-loop)", "d-loop control region", "similar to control region", "non coding region", "conrol region", "region: control region", "a+t-rich region", "putative control region 2", "d-loop region (= control reagion)", "pseudo contorl region", "a+t rich region", "c-rich region", "largest non-coding region", "d-loop", "d loop", "control region cr", "the control region", "control region c-rich sequence", "control region coretas sequence", "putative d-loop/control region", "d-loop containing region", "a+t rich", "at-rich region"],
};

const ChloroplastGenes = {
    "accD": ["acetyl-CoA carboxylase beta subunit", "acetyl-CoA carboxylase carboxyl transferase subunit beta", "beta-carboxyltransferase", "AccD", "acetyl-CoA carboxylase carboxyltransferase beta subunit", "accD protein", "Acetyl-CoA carboxylase, biotin carboxylase", "carboxyl transferase subunit beta", "acetyl-CoAcarboxylase beta subunit", "accD gene product", "Accda; acetyl-CoA carboxylase beta subunit"],
    "atpA": ["ATP synthase CF1 alpha subunit", "ATP synthase alpha subunit", "CF1 alpha subunit", "ATP synthase CF1 alpha chain", "AtpA", "atpA protein", "ATP synthase subunit alpha", "F1-ATPase alpha subunit", "ATP synthase alpha chain", "CF1-alpha", "Atpa; ATP synthase CF1 alpha subunit"],
    "atpB": ["ATP synthase CF1 beta subunit", "ATP synthase beta subunit", "CF1 beta subunit", "ATP synthase CF1 beta chain", "AtpB", "atpB protein", "ATP synthase subunit beta", "F1-ATPase beta subunit", "ATP synthase beta chain", "CF1-beta", "Atpb; ATP synthase CF1 beta subunit"],
    "atpE": ["ATP synthase CF1 epsilon subunit", "ATP synthase epsilon subunit", "CF1 epsilon subunit", "AtpE", "atpE protein", "ATP synthase subunit epsilon", "Atpe; ATP synthase CF1 epsilon subunit"],
    "atpF": ["ATP synthase CF0 subunit I", "ATP synthase CF0 I subunit", "AtpF", "atpF protein", "ATP synthase subunit b", "CF0 subunit I", "Atpf; ATP synthase CF0 subunit I"],
    "atpH": ["ATP synthase CF0 subunit III", "ATP synthase CF0 III subunit", "AtpH", "atpH protein", "CF0 subunit III", "ATP synthase subunit c", "Atph; ATP synthase CF0 subunit III"],
    "atpI": ["ATP synthase CF0 subunit IV", "ATP synthase CF0 IV subunit", "AtpI", "atpI protein", "CF0 subunit IV", "ATP synthase subunit a", "Atpi; ATP synthase CF0 subunit IV"],
    "ccsA": ["cytochrome c biogenesis protein", "CcsA", "cytochrome c heme attachment protein", "ccsA protein", "cytochrome c assembly protein"],
    "cemA": ["envelope membrane protein", "CemA", "cemA protein", "inner envelope membrane protein", "chloroplast envelope membrane protein"],
    "chlB": ["light-independent protochlorophyllide reductase subunit B", "ChlB", "protochlorophyllide reductase subunit B", "chlB protein"],
    "chlL": ["light-independent protochlorophyllide reductase subunit L", "ChlL", "protochlorophyllide reductase subunit L", "chlL protein"],
    "chlN": ["light-independent protochlorophyllide reductase subunit N", "ChlN", "protochlorophyllide reductase subunit N", "chlN protein"],
    "clpP": ["ATP-dependent Clp protease proteolytic subunit", "ClpP", "clpP protein", "Clp protease", "ATP-dependent Clp protease", "proteolytic subunit of Clp protease", "clpP1 protein"],
    "clpP1": ["ATP-dependent Clp protease proteolytic subunit 1", "ClpP1"],
    "infA": ["translation initiation factor 1", "IF1", "InfA", "infA protein", "translational initiation factor 1", "Infa; translation initiation factor 1", "initiation factor 1"],
    "matK": ["maturase K", "MatK", "maturase", "matK protein", "intron maturase", "maturaseK", "mat K", "maturase k protein", "maturaseK protein", "maturase enzyme", "MatK; maturase K", "Maturase K protein", "maturase type II intron", "maturase-K", "maturase protein K", "maturase-like protein", "maturase K-like protein", "intron-encoded maturase", "Group II intron maturase", "maturase protein"],
    "ndhA": ["NADH dehydrogenase subunit A", "NADH-plastoquinone oxidoreductase subunit A", "NdhA", "NAD(P)H dehydrogenase subunit A", "NADH dehydogenase subunit A", "ndhA protein", "NADH dehydrogenase ND1", "Ndha; NADH dehydrogenase subunit A"],
    "ndhB": ["NADH dehydrogenase subunit B", "NADH-plastoquinone oxidoreductase subunit B", "NdhB", "NAD(P)H dehydrogenase subunit B", "ndhB protein", "NADH dehydrogenase ND2", "Ndhb; NADH dehydrogenase subunit B"],
    "ndhC": ["NADH dehydrogenase subunit C", "NADH-plastoquinone oxidoreductase subunit C", "NdhC", "NAD(P)H dehydrogenase subunit C", "ndhC protein", "NADH dehydrogenase ND3", "Ndhc; NADH dehydrogenase subunit C"],
    "ndhD": ["NADH dehydrogenase subunit D", "NADH-plastoquinone oxidoreductase subunit D", "NdhD", "NAD(P)H dehydrogenase subunit D", "ndhD protein", "NADH dehydrogenase ND4", "Ndhd; NADH dehydrogenase subunit D"],
    "ndhE": ["NADH dehydrogenase subunit E", "NADH-plastoquinone oxidoreductase subunit E", "NdhE", "NAD(P)H dehydrogenase subunit E", "ndhE protein", "NADH dehydrogenase ND4L", "Ndhe; NADH dehydrogenase subunit E"],
    "ndhF": ["NADH dehydrogenase subunit F", "NADH-plastoquinone oxidoreductase subunit F", "NdhF", "NAD(P)H dehydrogenase subunit F", "ndhF protein", "NADH dehydrogenase ND5", "Ndhf; NADH dehydrogenase subunit F", "NADH dehydrogenase subunit 5"],
    "ndhG": ["NADH dehydrogenase subunit G", "NADH-plastoquinone oxidoreductase subunit G", "NdhG", "NAD(P)H dehydrogenase subunit G", "ndhG protein", "NADH dehydrogenase ND6", "Ndhg; NADH dehydrogenase subunit G"],
    "ndhH": ["NADH dehydrogenase subunit H", "NADH-plastoquinone oxidoreductase subunit H", "NdhH", "NAD(P)H dehydrogenase subunit H", "ndhH protein", "Ndhh; NADH dehydrogenase subunit H"],
    "ndhI": ["NADH dehydrogenase subunit I", "NADH-plastoquinone oxidoreductase subunit I", "NdhI", "NAD(P)H dehydrogenase subunit I", "ndhI protein", "Ndhi; NADH dehydrogenase subunit I"],
    "ndhJ": ["NADH dehydrogenase subunit J", "NADH-plastoquinone oxidoreductase subunit J", "NdhJ", "NAD(P)H dehydrogenase subunit J", "ndhJ protein", "Ndhj; NADH dehydrogenase subunit J"],
    "ndhK": ["NADH dehydrogenase subunit K", "NADH-plastoquinone oxidoreductase subunit K", "NdhK", "NAD(P)H dehydrogenase subunit K", "ndhK protein", "Ndhk; NADH dehydrogenase subunit K"],
    "petA": ["cytochrome f", "apocytochrome f", "PetA", "cytochrome f apoprotein", "petA protein", "Peta; cytochrome f"],
    "petB": ["cytochrome b6", "apocytochrome b6", "PetB", "cytochrome b6 apoprotein", "petB protein", "Petb; cytochrome b6"],
    "petD": ["cytochrome b6/f complex subunit IV", "PetD", "subunit IV of cytochrome b6/f complex", "petD protein", "Petd; cytochrome b6/f complex subunit IV"],
    "petG": ["cytochrome b6/f complex subunit V", "PetG", "petG protein", "Petg; cytochrome b6/f complex subunit V"],
    "petL": ["cytochrome b6/f complex subunit VI", "PetL", "petL protein", "Petl; cytochrome b6/f complex subunit VI"],
    "petN": ["cytochrome b6/f complex subunit VIII", "PetN", "petN protein", "Petn; cytochrome b6/f complex subunit VIII"],
    "psaA": ["photosystem I P700 chlorophyll a apoprotein A1", "PsaA", "PSI-A", "photosystem I subunit A", "psaA protein", "P700 apoprotein A1", "photosystem I reaction center subunit II", "Psaa; photosystem I P700 chlorophyll a apoprotein A1"],
    "psaB": ["photosystem I P700 chlorophyll a apoprotein A2", "PsaB", "PSI-B", "photosystem I subunit B", "psaB protein", "P700 apoprotein A2", "photosystem I reaction center subunit I", "Psab; photosystem I P700 chlorophyll a apoprotein A2"],
    "psaC": ["photosystem I subunit VII", "PsaC", "PSI-C", "psaC protein", "photosystem I iron-sulfur center", "Psac; photosystem I subunit VII"],
    "psaI": ["photosystem I subunit VIII", "PsaI", "PSI-I", "psaI protein", "Psai; photosystem I subunit VIII"],
    "psaJ": ["photosystem I subunit IX", "PsaJ", "PSI-J", "psaJ protein", "Psaj; photosystem I subunit IX"],
    "psbA": ["photosystem II protein D1", "PsbA", "PSII-D1", "D1 protein", "photosystem II D1 protein", "psbA protein", "photosystem II reaction center protein D1", "Psba; photosystem II protein D1", "photosystem II 32 kDa protein"],
    "psbB": ["photosystem II CP47 chlorophyll apoprotein", "PsbB", "PSII-B", "CP47", "photosystem II 47 kDa protein", "psbB protein", "photosystem II CP47 protein", "Psbb; photosystem II CP47 chlorophyll apoprotein"],
    "psbC": ["photosystem II CP43 chlorophyll apoprotein", "PsbC", "PSII-C", "CP43", "photosystem II 43 kDa protein", "psbC protein", "photosystem II CP43 protein", "Psbc; photosystem II CP43 chlorophyll apoprotein"],
    "psbD": ["photosystem II protein D2", "PsbD", "PSII-D2", "D2 protein", "photosystem II D2 protein", "psbD protein", "photosystem II reaction center protein D2", "Psbd; photosystem II protein D2"],
    "psbE": ["photosystem II cytochrome b559 alpha subunit", "PsbE", "cytochrome b559 alpha subunit", "psbE protein", "Psbe; photosystem II cytochrome b559 alpha subunit"],
    "psbF": ["photosystem II cytochrome b559 beta subunit", "PsbF", "cytochrome b559 beta subunit", "psbF protein", "Psbf; photosystem II cytochrome b559 beta subunit"],
    "psbH": ["photosystem II phosphoprotein", "PsbH", "photosystem II 10 kDa phosphoprotein", "psbH protein", "Psbh; photosystem II phosphoprotein"],
    "psbI": ["photosystem II protein I", "PsbI", "photosystem II reaction center subunit I", "psbI protein", "Psbi; photosystem II protein I"],
    "psbJ": ["photosystem II protein J", "PsbJ", "psbJ protein", "Psbj; photosystem II protein J"],
    "psbK": ["photosystem II protein K", "PsbK", "psbK protein", "Psbk; photosystem II protein K"],
    "psbL": ["photosystem II protein L", "PsbL", "psbL protein", "Psbl; photosystem II protein L"],
    "psbM": ["photosystem II protein M", "PsbM", "psbM protein", "Psbm; photosystem II protein M"],
    "psbN": ["photosystem II protein N", "PsbN", "psbN protein", "Psbn; photosystem II protein N"],
    "psbT": ["photosystem II protein T", "PsbT", "psbT protein", "Psbt; photosystem II protein T"],
    "psbZ": ["photosystem II protein Z", "PsbZ", "psbZ protein", "Psbz; photosystem II protein Z"],
    "rbcL": ["ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit", "ribulose 1,5-bisphosphate carboxylase/oxygenase large subunit", "ribulose bisophosphate carboxylase", "ribulose bisphosphate carboxylase large subunit", "Rubisco", "ribulose-1,5-bisphosphate carboxylase/oxygenase", "large subunit of riblose-1,5-bisphosphate carboxylase/oxygenase", "ribulose bisphosphate carboxylase large chain", "ribulose-1,5-bisphosphate carboxylase", "ribulose-bisphosphate carboxylase large chain", "ribulose-bisphosphate carboxylase large subunit", "RuBisCO large chain", "ribulose bisphosphate carboxylase/oxygenase", "ribulose-bisphosphate carboxylase/oxygenase", "ribulose-bisphosphate carboxylase", "ribulose 1,5-bisphosphate carboxylase/oxygenase", "ribulose 1,5-bisphosphate carboxylase large subunit", "ribulose 1,5-bisphosphate carboxylase", "Rubisco large subunit", "large subunit of RuBisCO", "rbcL protein", "RuBisCO,Ribulose Bisphosphate Carboxylase/Oxygenase large subunit", "ribulose 1,5-bisphosphate carboxylase oxygenase large subunit", "ribulose-1,5-bisphosphate carboxylase oxygenase large subunit", "ribulose-1,5-bisphosphate carboxylase oxygenase", "ribulose bisphosphate carboxylase large subunit precursor", "Ribulose-1,5-Bisphosphate Carboxylase-Oxygenase"],
    "rpl2": ["ribosomal protein L2", "50S ribosomal protein L2", "rpl2 protein", "large subunit ribosomal protein L2"],
    "rpl14": ["ribosomal protein L14", "50S ribosomal protein L14", "rpl14 protein", "large subunit ribosomal protein L14"],
    "rpl16": ["ribosomal protein L16", "50S ribosomal protein L16", "rpl16 protein", "large subunit ribosomal protein L16"],
    "rpl20": ["ribosomal protein L20", "50S ribosomal protein L20", "rpl20 protein", "large subunit ribosomal protein L20"],
    "rpl22": ["ribosomal protein L22", "50S ribosomal protein L22", "rpl22 protein"],
    "rpl23": ["ribosomal protein L23", "50S ribosomal protein L23", "rpl23 protein", "large subunit ribosomal protein L23"],
    "rpl32": ["ribosomal protein L32", "50S ribosomal protein L32", "rpl32 protein"],
    "rpl33": ["ribosomal protein L33", "50S ribosomal protein L33", "rpl33 protein", "large subunit ribosomal protein L33"],
    "rpl36": ["ribosomal protein L36", "50S ribosomal protein L36", "rpl36 protein"],
    "rpoA": ["RNA polymerase alpha subunit", "RpoA", "DNA-directed RNA polymerase alpha subunit", "RNA polymerase subunit alpha", "rpoA protein"],
    "rpoB": ["RNA polymerase beta subunit", "RpoB", "DNA-directed RNA polymerase beta subunit", "RNA polymerase subunit beta", "rpoB protein"],
    "rpoC1": ["RNA polymerase beta' subunit", "RpoC1", "DNA-directed RNA polymerase subunit gamma", "RNA polymerase subunit C1", "rpoC1 protein"],
    "rpoC2": ["RNA polymerase beta'' subunit", "RpoC2", "DNA-directed RNA polymerase beta'' subunit", "RNA polymerase subunit beta''", "rpoC2 protein"],
    "rps2": ["ribosomal protein S2", "30S ribosomal protein S2", "rps2 protein", "small subunit ribosomal protein S2"],
    "rps3": ["ribosomal protein S3", "30S ribosomal protein S3", "rps3 protein"],
    "rps4": ["ribosomal protein S4", "30S ribosomal protein S4", "rps4 protein", "small subunit ribosomal protein S4"],
    "rps7": ["ribosomal protein S7", "30S ribosomal protein S7", "rps7 protein", "small subunit ribosomal protein S7"],
    "rps8": ["ribosomal protein S8", "30S ribosomal protein S8", "rps8 protein"],
    "rps11": ["ribosomal protein S11", "30S ribosomal protein S11", "rps11 protein", "small subunit ribosomal protein S11"],
    "rps12": ["ribosomal protein S12", "30S ribosomal protein S12", "rps12 protein", "small subunit ribosomal protein S12"],
    "rps14": ["ribosomal protein S14", "30S ribosomal protein S14", "rps14 protein"],
    "rps15": ["ribosomal protein S15", "30S ribosomal protein S15", "rps15 protein", "small subunit ribosomal protein S15"],
    "rps16": ["ribosomal protein S16", "30S ribosomal protein S16", "rps16"],
    "rps18": ["ribosomal protein S18", "30S ribosomal protein S18", "rps18 protein"],
    "rps19": ["ribosomal protein S19", "30S ribosomal protein S19", "rps19 protein"],
    "rrn16S": ["16S ribosomal RNA", "small subunit ribosomal RNA", "rrn16S", "ribosomal RNA"],
    "rrn23S": ["23S ribosomal RNA", "large subunit ribosomal RNA", "ribosomal 23S RNA"],
    "rrn4.5S": ["4.5S ribosomal RNA", "ribosomal 4.5S RNA"],
    "rrn5S": ["5S ribosomal RNA", "ribosomal 5S RNA"],
    "ycf1": ["Ycf1 protein", "ycf1"],
    "ycf2": ["Ycf2 protein", "ycf2", "Ycf2"],
};

// ========================================================================
// Gene Name Standardization
// ========================================================================
function standardizeGeneName(rawName, product, dataType) {
    const dict = dataType === "cp" ? ChloroplastGenes : MitochondrialGenes;
    const altDict = dataType === "cp" ? MitochondrialGenes : ChloroplastGenes;

    // First try matching product against synonym lists (case-insensitive)
    if (product) {
        const productLower = product.toLowerCase().trim();
        for (const [key, synonyms] of Object.entries(dict)) {
            for (const syn of synonyms) {
                if (syn.toLowerCase() === productLower) {
                    return key;
                }
            }
        }
        // Try alt dictionary
        for (const [key, synonyms] of Object.entries(altDict)) {
            for (const syn of synonyms) {
                if (syn.toLowerCase() === productLower) {
                    return key;
                }
            }
        }
    }

    // Then try rawName (gene qualifier) against dictionary keys directly
    if (rawName) {
        for (const key of Object.keys(dict)) {
            if (key.toUpperCase() === rawName.toUpperCase()) return key;
        }
        for (const key of Object.keys(altDict)) {
            if (key.toUpperCase() === rawName.toUpperCase()) return key;
        }
        // Try rawName against synonyms too
        const rawLower = rawName.toLowerCase().trim();
        for (const [key, synonyms] of Object.entries(dict)) {
            for (const syn of synonyms) {
                if (syn.toLowerCase() === rawLower) {
                    return key;
                }
            }
        }
        for (const [key, synonyms] of Object.entries(altDict)) {
            for (const syn of synonyms) {
                if (syn.toLowerCase() === rawLower) {
                    return key;
                }
            }
        }
    }

    // Fallback: return raw name
    return rawName || product || null;
}

// ========================================================================
// tRNA Name Standardization (Leu1/Leu2, Ser1/Ser2)
// ========================================================================
function standardizeTrnaName(feat) {
    const product = feat.qualifiers.product || "";
    const gene = feat.qualifiers.gene || "";
    const note = feat.qualifiers.note || "";
    const anticodon = feat.qualifiers.anticodon || "";
    const raw = product || gene;

    if (!raw) return null;

    const rawLower = raw.toLowerCase();
    const geneLower = gene.toLowerCase();
    const combined = (anticodon + " " + product + " " + note + " " + gene).toLowerCase();

    // Detect Leucine tRNAs (product "tRNA-Leu" or gene "trnL*")
    if (rawLower.includes("leu") || geneLower.startsWith("trnl")) {
        // Check gene qualifier for direct numbering (trnL1, trnL2, trnL(UAA), trnL(UAG))
        if (geneLower === "trnl1" || geneLower.includes("trnl(uaa)") || geneLower.includes("trnl(taa)")) {
            return "tRNA-Leu1";
        }
        if (geneLower === "trnl2" || geneLower.includes("trnl(uag)") || geneLower.includes("trnl(tag)")) {
            return "tRNA-Leu2";
        }
        // Check anticodon and other qualifiers
        // Leu1 = trnL(UAA) = anticodon TAA
        if (combined.includes("taa") || combined.includes("uaa")) {
            return "tRNA-Leu1";
        }
        // Leu2 = trnL(UAG) = anticodon TAG, recognizes CUN codons
        if (combined.includes("tag") || combined.includes("uag") || combined.includes("cun")) {
            return "tRNA-Leu2";
        }
        return "tRNA-Leu1"; // Default to Leu1 when anticodon cannot be determined
    }

    // Detect Serine tRNAs (product "tRNA-Ser" or gene "trnS*")
    if (rawLower.includes("ser") || geneLower.startsWith("trns")) {
        // Check gene qualifier for direct numbering (trnS1, trnS2, trnS(UGA), trnS(GCU))
        if (geneLower === "trns1" || geneLower.includes("trns(uga)") || geneLower.includes("trns(tga)")) {
            return "tRNA-Ser1";
        }
        if (geneLower === "trns2" || geneLower.includes("trns(gcu)") || geneLower.includes("trns(gct)")) {
            return "tRNA-Ser2";
        }
        // Check anticodon and other qualifiers
        // Ser1 = trnS(UGA) = anticodon TGA (recognizes UCN)
        if (combined.includes("tga") || combined.includes("uga") || combined.includes("ucn")) {
            return "tRNA-Ser1";
        }
        // Ser2 = trnS(GCU) = anticodon GCT (recognizes AGN/AGY)
        if (combined.includes("gct") || combined.includes("gcu") || combined.includes("agn") || combined.includes("agy")) {
            return "tRNA-Ser2";
        }
        return "tRNA-Ser1"; // Default to Ser1 when anticodon cannot be determined
    }

    // For other tRNAs, return a clean standardized name
    const aaMatch = raw.match(/tRNA-(\w+)/i);
    if (aaMatch) {
        return `tRNA-${aaMatch[1].charAt(0).toUpperCase()}${aaMatch[1].slice(1).toLowerCase()}`;
    }

    return raw;
}

// ========================================================================
// Auto Data Type Detection
// ========================================================================
function detectDataType() {
    let mtScore = 0;
    let cpScore = 0;

    for (const record of state.records) {
        for (const feat of record.features) {
            const rawGene = feat.qualifiers.gene || "";
            const rawProduct = feat.qualifiers.product || "";
            const combined = (rawGene + " " + rawProduct).toLowerCase();

            // Check for mt-specific markers
            for (const key of Object.keys(MitochondrialGenes)) {
                if (key.toLowerCase() === rawGene.toLowerCase() ||
                    key.toLowerCase() === rawProduct.toLowerCase()) {
                    mtScore++;
                    break;
                }
            }

            // Check for cp-specific markers
            for (const key of Object.keys(ChloroplastGenes)) {
                if (key.toLowerCase() === rawGene.toLowerCase() ||
                    key.toLowerCase() === rawProduct.toLowerCase()) {
                    cpScore++;
                    break;
                }
            }
        }
    }

    state.detectedDataType = cpScore > mtScore ? "cp" : "mt";
    return state.detectedDataType;
}

// ========================================================================
// GenBank Parser
// ========================================================================
function parseGenBank(text) {
    const result = { accession: "", organism: "", taxonomy: "", features: [], sequence: "", taxonomyRanks: { kingdom: "", phylum: "", class: "", order: "", family: "" } };

    const locusMatch = text.match(/^LOCUS\s+(\S+)/m);
    if (locusMatch) result.accession = locusMatch[1];

    const versionMatch = text.match(/^VERSION\s+(\S+)/m);
    if (versionMatch) result.accession = versionMatch[1];

    const orgMatch = text.match(/^\s+ORGANISM\s+(.+)/m);
    if (orgMatch) result.organism = orgMatch[1].trim();

    // Parse taxonomy lineage (lines after ORGANISM until next section)
    const orgIdx = text.indexOf("ORGANISM");
    if (orgIdx !== -1) {
        const afterOrg = text.substring(orgIdx);
        const orgLines = afterOrg.split("\n");
        const taxLines = [];
        for (let i = 1; i < orgLines.length; i++) {
            const line = orgLines[i];
            // Taxonomy lines are indented and contain semicolons or end with period
            if (/^\s{10,}/.test(line) && (line.includes(";") || line.trim().endsWith("."))) {
                taxLines.push(line.trim().replace(/\.$/, ""));
            } else {
                break;
            }
        }
        result.taxonomy = taxLines.join(" ");
    }

    const featStart = text.indexOf("FEATURES");
    const originStart = text.indexOf("ORIGIN");
    if (featStart !== -1 && originStart !== -1) {
        const featText = text.substring(featStart, originStart);
        result.features = parseFeatures(featText);
    }

    if (originStart !== -1) {
        const originText = text.substring(originStart);
        const seqLines = originText.split("\n").slice(1);
        const seqParts = [];
        for (const line of seqLines) {
            if (line.startsWith("//")) break;
            seqParts.push(line.replace(/[\s\d]/g, ""));
        }
        result.sequence = seqParts.join("").toLowerCase();
    }

    return result;
}

function extractTaxonomy(record) {
    const organism = record.editedOrganism || record.organism || "";
    const parts = organism.split(/\s+/);
    const genus = parts[0] || "";
    const species = parts.slice(1).join(" ") || "";

    const ranks = record.taxonomyRanks || {};

    const kingdom = ranks.kingdom || "";
    const phylum = ranks.phylum || "";
    const klass = ranks.class || "";
    const order = ranks.order || "";

    // Family: prefer taxonomyRanks, then editedFamily, then lineage parse
    let family = ranks.family || record.editedFamily || "";
    if (!family && record.taxonomy) {
        const taxParts = record.taxonomy.split(/;\s*/);
        for (const part of taxParts) {
            const trimmed = part.trim();
            if (trimmed.match(/idae$|aceae$|ales$/i)) {
                if (trimmed.match(/idae$|aceae$/i)) {
                    family = trimmed;
                } else if (!family) {
                    family = trimmed;
                }
            }
        }
    }

    return { kingdom, phylum, class: klass, order, family, genus, species };
}

function countFeatureTypes(record) {
    let pcgs = 0, rrnas = 0, trnas = 0;
    for (const feat of record.features) {
        if (feat.type === "CDS") pcgs++;
        else if (feat.type === "rRNA") rrnas++;
        else if (feat.type === "tRNA") trnas++;
    }
    return { pcgs, rrnas, trnas };
}

function parseFeatures(featText) {
    const features = [];
    const lines = featText.split("\n");
    let current = null;
    let lastQualKey = null;

    for (let i = 1; i < lines.length; i++) {
        const line = lines[i];
        if (!line || line.trim() === "") continue;

        const featureMatch = line.match(/^     (\S+)\s+([\w<>().,:]+)/);
        if (featureMatch) {
            if (current) features.push(current);
            current = {
                type: featureMatch[1],
                locationStr: featureMatch[2],
                qualifiers: {},
            };
            lastQualKey = null;

            let j = i + 1;
            while (j < lines.length && lines[j].match(/^\s{21}[^/]/) && !lines[j].match(/^\s{21}\//)) {
                current.locationStr += lines[j].trim();
                j++;
            }
            continue;
        }

        if (current) {
            const qualMatch = line.match(/^\s{21}\/(\w+)="?([^"]*)"?$/);
            const qualMatchStart = line.match(/^\s{21}\/(\w+)="([^"]*)$/);
            if (qualMatch) {
                current.qualifiers[qualMatch[1]] = qualMatch[2];
                lastQualKey = qualMatch[1];
            } else if (qualMatchStart) {
                current.qualifiers[qualMatchStart[1]] = qualMatchStart[2];
                lastQualKey = qualMatchStart[1];
            } else if (lastQualKey && line.match(/^\s{21}[^/]/)) {
                const val = line.trim().replace(/"$/, "");
                current.qualifiers[lastQualKey] = (current.qualifiers[lastQualKey] || "") + " " + val;
            }
        }
    }
    if (current) features.push(current);
    return features;
}

function parseLocation(locStr, seqLen) {
    const segments = [];
    const isComplement = locStr.startsWith("complement(");
    let inner = locStr;
    if (isComplement) inner = inner.replace(/^complement\(/, "").replace(/\)$/, "");

    const isJoin = inner.startsWith("join(");
    if (isJoin) inner = inner.replace(/^join\(/, "").replace(/\)$/, "");

    const parts = inner.split(",");
    for (const part of parts) {
        const rangeMatch = part.trim().match(/<?\s*(\d+)\s*\.\.\s*>?\s*(\d+)/);
        if (rangeMatch) {
            segments.push({
                start: parseInt(rangeMatch[1]) - 1,
                end: parseInt(rangeMatch[2]),
                complement: isComplement,
            });
        }
    }
    return segments;
}

function extractSequence(fullSeq, locationStr) {
    const segments = parseLocation(locationStr, fullSeq.length);
    let seq = "";
    for (const seg of segments) {
        seq += fullSeq.substring(seg.start, seg.end);
    }
    if (segments.length > 0 && segments[0].complement) {
        seq = reverseComplement(seq);
    }
    return seq.toUpperCase();
}

function reverseComplement(seq) {
    const comp = {
        a: "t", t: "a", c: "g", g: "c", n: "n",
        A: "T", T: "A", C: "G", G: "C", N: "N"
    };
    return seq.split("").reverse().map(c => comp[c] || c).join("");
}

// ========================================================================
// NCBI Fetcher
// ========================================================================
const NCBI_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";
const FETCH_DELAY = 350;

async function fetchAccession(accession) {
    const cached = await cacheGet(accession);
    if (cached) {
        log(`Cache hit: ${accession}`);
        return cached;
    }

    const url = `${NCBI_BASE}?db=nucleotide&id=${encodeURIComponent(accession)}&rettype=gb&retmode=text`;
    log(`Fetching ${accession} from NCBI...`);

    const resp = await fetch(url);
    if (!resp.ok) throw new Error(`NCBI returned ${resp.status} for ${accession}`);

    const text = await resp.text();
    if (text.includes("Error") && text.length < 500) throw new Error(`NCBI error for ${accession}: ${text.trim()}`);

    const parsed = parseGenBank(text);
    parsed.source = "ncbi";

    await cachePut(accession, parsed, text);
    return parsed;
}

async function fetchMultiple(accessions) {
    const results = [];
    const statusEl = document.getElementById("fetchStatus");
    statusEl.classList.remove("hidden");

    showProgress("Fetching from NCBI", `0 of ${accessions.length}`);

    for (let i = 0; i < accessions.length; i++) {
        const acc = accessions[i].trim();
        if (!acc) continue;
        statusEl.textContent = `Fetching ${i + 1}/${accessions.length}: ${acc}...`;
        updateProgress(i + 1, accessions.length, acc);
        try {
            const record = await fetchAccession(acc);
            results.push(record);
            await new Promise(r => setTimeout(r, FETCH_DELAY));
        } catch (e) {
            log(`Error fetching ${acc}: ${e.message}`);
        }
    }
    statusEl.textContent = `Fetched ${results.length}/${accessions.length} records`;
    setTimeout(() => statusEl.classList.add("hidden"), 3000);

    hideProgress(`Loaded ${results.length} record${results.length !== 1 ? 's' : ''}`);

    return results;
}

// ========================================================================
// IndexedDB Cache
// ========================================================================
const DB_NAME = "splace_cache";
const DB_VERSION = 1;
const STORE_NAME = "genbank_records";
const CACHE_TTL = 7 * 24 * 60 * 60 * 1000;

function openDB() {
    return new Promise((resolve, reject) => {
        const req = indexedDB.open(DB_NAME, DB_VERSION);
        req.onupgradeneeded = () => {
            const db = req.result;
            if (!db.objectStoreNames.contains(STORE_NAME)) {
                db.createObjectStore(STORE_NAME, { keyPath: "accession" });
            }
        };
        req.onsuccess = () => resolve(req.result);
        req.onerror = () => reject(req.error);
    });
}

async function cacheGet(accession) {
    try {
        const db = await openDB();
        return new Promise((resolve) => {
            const tx = db.transaction(STORE_NAME, "readonly");
            const store = tx.objectStore(STORE_NAME);
            const req = store.get(accession);
            req.onsuccess = () => {
                const record = req.result;
                if (record && (Date.now() - record.timestamp < CACHE_TTL)) {
                    resolve(record.data);
                } else {
                    resolve(null);
                }
            };
            req.onerror = () => resolve(null);
        });
    } catch { return null; }
}

async function cachePut(accession, data, rawText) {
    try {
        const db = await openDB();
        const tx = db.transaction(STORE_NAME, "readwrite");
        tx.objectStore(STORE_NAME).put({
            accession,
            data,
            rawText,
            timestamp: Date.now(),
        });
    } catch (e) {
        log(`Cache write error: ${e.message}`);
    }
}

// ========================================================================
// UI Rendering
// ========================================================================
function renderRecords() {
    const section = document.getElementById("recordsSection");
    const tbody = document.getElementById("recordsTableBody");
    const count = document.getElementById("recordCount");

    if (state.records.length === 0) {
        section.classList.add("hidden");
        document.getElementById("featureTypesSection").classList.add("hidden");
        document.getElementById("genesSection").classList.add("hidden");
        document.getElementById("downloadSection").classList.add("hidden");
        return;
    }

    // Auto-detect data type
    detectDataType();

    section.classList.remove("hidden");
    count.textContent = `(${state.records.length})`;

    // Sort records alphabetically by organism name
    const sortedIndices = state.records
        .map((r, i) => ({ record: r, originalIndex: i }))
        .sort((a, b) => (a.record.organism || "").localeCompare(b.record.organism || ""));

    tbody.innerHTML = sortedIndices.map(({ record: r, originalIndex: i }) => {
        const tax = extractTaxonomy(r);
        const counts = countFeatureTypes(r);
        const sourceClass = r.source === 'ncbi' ? 'bg-splace-blue-50 text-splace-blue-600' : 'bg-gray-100 text-gray-600';
        const sourceLabel = r.source === 'ncbi' ? 'NCBI' : 'File';
        const abbrevSpecies = tax.genus ? `${tax.genus.charAt(0)}.${tax.species ? ' ' + tax.species : ''}` : '';

        return `
            <tr class="group">
                <td class="px-3 py-2 font-mono text-xs text-gray-600" data-col="accession">${r.accession}</td>
                <td class="px-3 py-2 text-xs text-gray-500" data-col="kingdom">${tax.kingdom || "—"}</td>
                <td class="px-3 py-2 text-xs text-gray-500" data-col="phylum">${tax.phylum || "—"}</td>
                <td class="px-3 py-2 text-xs text-gray-500" data-col="class">${tax.class || "—"}</td>
                <td class="px-3 py-2 text-xs text-gray-500" data-col="order">${tax.order || "—"}</td>
                <td class="px-3 py-2 text-xs text-gray-500" data-col="family">${tax.family || "—"}</td>
                <td class="px-3 py-2 text-xs text-gray-500" data-col="genus"><em>${tax.genus}</em></td>
                <td class="px-3 py-2 text-xs text-gray-500" data-col="species"><em>${abbrevSpecies}</em></td>
                <td class="px-3 py-2 text-center text-xs text-gray-600" data-col="pcgs">${counts.pcgs}</td>
                <td class="px-3 py-2 text-center text-xs text-gray-600" data-col="rrnas">${counts.rrnas}</td>
                <td class="px-3 py-2 text-center text-xs text-gray-600" data-col="trnas">${counts.trnas}</td>
                <td class="px-3 py-2 text-center" data-col="source"><span class="text-xs px-1.5 py-0.5 rounded ${sourceClass}">${sourceLabel}</span></td>
                <td class="px-3 py-2 text-right whitespace-nowrap">
                    <button onclick="editRecord(${i})" class="text-gray-400 hover:text-splace-blue-600 transition-colors text-sm leading-none mr-2" title="Edit"><i class="fa-solid fa-pen-to-square"></i></button>
                    <button onclick="removeRecord(${i})" class="text-gray-400 hover:text-red-500 transition-colors text-lg leading-none" title="Remove"><i class="fa-solid fa-xmark"></i></button>
                </td>
            </tr>
        `;
    }).join("");

    renderGeneSelection();
    applyColumnVisibility();
}

function applyColumnVisibility() {
    const table = document.querySelector(".records-table");
    if (!table) return;
    const cells = table.querySelectorAll("[data-col]");
    cells.forEach(cell => {
        const col = cell.dataset.col;
        cell.style.display = state.hiddenColumns.has(col) ? "none" : "";
    });
    // Sync eye icons
    document.querySelectorAll(".col-eye-toggle[data-col-toggle]").forEach(btn => {
        const col = btn.dataset.colToggle;
        const icon = btn.querySelector("i");
        if (state.hiddenColumns.has(col)) {
            icon.className = "fa-solid fa-eye-slash";
            btn.classList.add("text-gray-300");
            btn.classList.remove("text-gray-500");
        } else {
            icon.className = "fa-solid fa-eye";
            btn.classList.remove("text-gray-300");
            btn.classList.add("text-gray-500");
        }
    });
}

function toggleColumn(colKey, visible) {
    if (visible) {
        state.hiddenColumns.delete(colKey);
    } else {
        state.hiddenColumns.add(colKey);
    }
    applyColumnVisibility();
}

function renderGeneSelection() {
    if (state.records.length === 0) return;

    document.getElementById("featureTypesSection").classList.remove("hidden");
    document.getElementById("genesSection").classList.remove("hidden");
    document.getElementById("downloadSection").classList.remove("hidden");

    const dataType = state.detectedDataType;

    // Collect all feature types and genes
    const featureTypes = new Map();
    const genesByType = new Map();

    for (const record of state.records) {
        const assignedTrnas = new Set();
        for (const feat of record.features) {
            featureTypes.set(feat.type, (featureTypes.get(feat.type) || 0) + 1);

            let geneName = null;

            if (feat.type === "tRNA") {
                geneName = standardizeTrnaName(feat);
                // Dedup: if same tRNA name already seen in this record, use alternate
                if (geneName && assignedTrnas.has(geneName)) {
                    if (geneName === "tRNA-Ser1") geneName = "tRNA-Ser2";
                    else if (geneName === "tRNA-Leu1") geneName = "tRNA-Leu2";
                }
                if (geneName) assignedTrnas.add(geneName);
            } else {
                const rawGene = feat.qualifiers.gene || null;
                const rawProduct = feat.qualifiers.product || null;
                if (!rawGene && !rawProduct) continue;
                geneName = standardizeGeneName(rawGene, rawProduct, dataType);
            }

            if (!geneName) continue;

            if (!genesByType.has(feat.type)) genesByType.set(feat.type, new Map());
            const geneMap = genesByType.get(feat.type);
            const rawProduct = feat.qualifiers.product || geneName;

            // Compute sequence length for this feature
            let seqLen = 0;
            try {
                const seq = extractSequence(record.sequence, feat.locationStr);
                seqLen = seq.length;
            } catch (e) { /* ignore */ }

            const existing = geneMap.get(geneName) || { product: rawProduct, count: 0, seqLengths: [] };
            existing.count++;
            if (seqLen > 0) existing.seqLengths.push(seqLen);
            geneMap.set(geneName, existing);
        }
    }

    // Render feature type switches
    const ftList = document.getElementById("featureTypesList");
    const relevantTypes = ["CDS", "rRNA", "tRNA"];
    const availableTypes = relevantTypes.filter(t => featureTypes.has(t));

    ftList.innerHTML = availableTypes.map(type => {
        const count = featureTypes.get(type) || 0;
        const active = state.selectedFeatureTypes.has(type);
        const info = FEATURE_TYPE_INFO[type] || `Feature type: ${type}`;
        const displayLabel = FEATURE_TYPE_LABELS[type] || type;
        return `
            <label class="feature-switch inline-flex items-center gap-2 cursor-pointer select-none" data-tippy-content="${info}">
                <span class="switch-track ${active ? 'active' : ''}" onclick="toggleFeatureType('${type}', !${active}); event.preventDefault();">
                    <span class="switch-knob"></span>
                </span>
                <span class="text-sm font-medium ${active ? 'text-gray-800' : 'text-gray-500'}">${displayLabel}</span>
                <span class="text-xs text-gray-400">(${count})</span>
            </label>
        `;
    }).join("");

    renderGenes(genesByType);

    // Init Tippy on feature switches
    if (typeof tippy !== "undefined") {
        tippy('[data-tippy-content]', {
            theme: 'splace',
            placement: 'top',
            arrow: true,
        });
    }
}

function renderGenes(genesByType) {
    if (!genesByType) {
        genesByType = new Map();
        const dataType = state.detectedDataType;
        for (const record of state.records) {
            const assignedTrnas = new Set();
            for (const feat of record.features) {
                let geneName = null;

                if (feat.type === "tRNA") {
                    geneName = standardizeTrnaName(feat);
                    if (geneName && assignedTrnas.has(geneName)) {
                        if (geneName === "tRNA-Ser1") geneName = "tRNA-Ser2";
                        else if (geneName === "tRNA-Leu1") geneName = "tRNA-Leu2";
                    }
                    if (geneName) assignedTrnas.add(geneName);
                } else {
                    const rawGene = feat.qualifiers.gene || null;
                    const rawProduct = feat.qualifiers.product || null;
                    if (!rawGene && !rawProduct) continue;
                    geneName = standardizeGeneName(rawGene, rawProduct, dataType);
                }

                if (!geneName) continue;
                if (!genesByType.has(feat.type)) genesByType.set(feat.type, new Map());
                const geneMap = genesByType.get(feat.type);
                const rawProduct = feat.qualifiers.product || geneName;

                let seqLen = 0;
                try {
                    const seq = extractSequence(record.sequence, feat.locationStr);
                    seqLen = seq.length;
                } catch (e) { /* ignore */ }

                const existing = geneMap.get(geneName) || { product: rawProduct, count: 0, seqLengths: [] };
                existing.count++;
                if (seqLen > 0) existing.seqLengths.push(seqLen);
                geneMap.set(geneName, existing);
            }
        }
    }

    const genesList = document.getElementById("genesList");
    const visibleGenes = new Map();

    for (const [type, geneMap] of genesByType) {
        if (!state.selectedFeatureTypes.has(type)) continue;
        for (const [name, info] of geneMap) {
            const existing = visibleGenes.get(name) || { product: info.product, count: 0, types: [], seqLengths: [] };
            existing.count += info.count;
            existing.types.push(type);
            existing.seqLengths.push(...(info.seqLengths || []));
            visibleGenes.set(name, existing);
        }
    }

    const sortedGenes = [...visibleGenes.entries()].sort((a, b) => a[0].localeCompare(b[0]));

    genesList.innerHTML = sortedGenes.map(([name, info]) => {
        const active = state.selectedGenes.has(name);
        const lengths = info.seqLengths;
        let avgLen = 0, minLen = 0, maxLen = 0;
        if (lengths.length > 0) {
            avgLen = Math.round(lengths.reduce((a, b) => a + b, 0) / lengths.length);
            minLen = Math.min(...lengths);
            maxLen = Math.max(...lengths);
        }

        const tooltipHtml = `<strong>${name}</strong><br>` +
            `Feature: ${info.types.join(", ")}<br>` +
            `Records: ${info.count}/${state.records.length}<br>` +
            (lengths.length > 0
                ? `Avg: ${avgLen.toLocaleString()} bp<br>Min: ${minLen.toLocaleString()} bp | Max: ${maxLen.toLocaleString()} bp`
                : "No sequence data");

        const escapedName = name.replace(/'/g, "\\'");

        return `
            <div class="gene-chip ${active ? 'active' : ''}" data-gene="${name}"
                onclick="toggleGene('${escapedName}', !state.selectedGenes.has('${escapedName}'))">
                <span class="switch-track-sm ${active ? 'active' : ''}">
                    <span class="switch-knob-sm"></span>
                </span>
                <span class="gene-name">${name}</span>
                <span class="gene-count">${info.count}/${state.records.length}</span>
                ${lengths.length > 0 ? `<span class="gene-len">~${avgLen} bp</span>` : ''}
            </div>
        `;
    }).join("");

    if (sortedGenes.length === 0) {
        genesList.innerHTML = '<p class="text-sm text-gray-400 italic">No genes found for selected feature types</p>';
    }

    // Init Tippy tooltips for gene chips
    if (typeof tippy !== "undefined") {
        // Destroy old instances
        document.querySelectorAll('.gene-chip').forEach(el => {
            if (el._tippy) el._tippy.destroy();
        });

        sortedGenes.forEach(([name, info]) => {
            const lengths = info.seqLengths;
            let avgLen = 0, minLen = 0, maxLen = 0;
            if (lengths.length > 0) {
                avgLen = Math.round(lengths.reduce((a, b) => a + b, 0) / lengths.length);
                minLen = Math.min(...lengths);
                maxLen = Math.max(...lengths);
            }

            const tooltipHtml = `<strong>${name}</strong><br>` +
                `Feature: ${info.types.join(", ")}<br>` +
                `Records: ${info.count}/${state.records.length}<br>` +
                (lengths.length > 0
                    ? `Avg: ${avgLen.toLocaleString()} bp<br>Min: ${minLen.toLocaleString()} bp | Max: ${maxLen.toLocaleString()} bp`
                    : "No sequence data");

            const el = document.querySelector(`.gene-chip[data-gene="${name}"]`);
            if (el) {
                tippy(el, {
                    content: tooltipHtml,
                    allowHTML: true,
                    theme: 'splace',
                    placement: 'top',
                    arrow: true,
                });
            }
        });
    }
}

// ========================================================================
// UI Event Handlers
// ========================================================================
function toggleFeatureType(type, checked) {
    if (checked) state.selectedFeatureTypes.add(type);
    else state.selectedFeatureTypes.delete(type);
    renderGeneSelection();
}

function toggleGene(name, checked) {
    if (checked) state.selectedGenes.add(name);
    else state.selectedGenes.delete(name);
    renderGenes();
}

window.removeRecord = function (index) {
    state.pendingRemoveIndex = index;
    const r = state.records[index];
    const tax = extractTaxonomy(r);
    document.getElementById("removeModalText").innerHTML =
        `Remove <em>${tax.genus} ${tax.species}</em> (${r.accession})? This action cannot be undone.`;
    document.getElementById("removeModal").classList.remove("hidden");
};

window.toggleFeatureType = toggleFeatureType;
window.toggleGene = toggleGene;

function selectAllGenes() {
    document.querySelectorAll(".gene-chip").forEach(el => {
        const name = el.dataset.gene;
        if (name) state.selectedGenes.add(name);
    });
    renderGenes();
}

function selectNoneGenes() {
    state.selectedGenes.clear();
    renderGenes();
}

function selectDefaultGenes() {
    state.selectedGenes.clear();
    // Ensure CDS feature type is enabled (default genes are PCGs)
    if (!state.selectedFeatureTypes.has("CDS")) {
        state.selectedFeatureTypes.add("CDS");
    }
    const defaults = state.detectedDataType === "mt" ? MT_DEFAULT_GENES : CP_DEFAULT_GENES;
    defaults.forEach(g => state.selectedGenes.add(g));
    renderGeneSelection();
}

function fillExampleAccessions() {
    document.getElementById("accessionInput").value = EXAMPLE_ACCESSIONS;
}

// ========================================================================
// Progress Modal
// ========================================================================
function showProgress(title, subtitle) {
    document.getElementById("progressTitle").textContent = title;
    document.getElementById("progressSubtitle").textContent = subtitle || "";
    document.getElementById("progressBar").style.width = "0%";
    document.getElementById("progressDetail").textContent = "";
    document.getElementById("progressIcon").className = "fa-solid fa-spinner fa-spin text-splace-blue-600 text-xl";
    document.getElementById("progressModal").classList.remove("hidden");
}

function updateProgress(current, total, detail) {
    const pct = Math.round((current / total) * 100);
    document.getElementById("progressBar").style.width = pct + "%";
    document.getElementById("progressSubtitle").textContent = `${current} of ${total}`;
    if (detail) document.getElementById("progressDetail").innerHTML = detail;
}

function hideProgress(successMsg) {
    document.getElementById("progressBar").style.width = "100%";
    document.getElementById("progressIcon").className = "fa-solid fa-circle-check text-green-600 text-xl";
    if (successMsg) document.getElementById("progressTitle").textContent = successMsg;
    document.getElementById("progressSubtitle").textContent = "";
    setTimeout(() => {
        document.getElementById("progressModal").classList.add("hidden");
    }, 1200);
}

// ========================================================================
// GBIF Taxonomy Lookup (dataFishing - Rabelo et al. 2025)
// ========================================================================
async function fetchTaxonomyFromGBIF(speciesName) {
    const url = `https://api.gbif.org/v1/species?name=${encodeURIComponent(speciesName)}`;
    const resp = await fetch(url, {
        method: "GET",
        headers: { "Accept": "application/json" }
    });
    if (!resp.ok) throw new Error(`GBIF returned ${resp.status}`);
    const data = await resp.json();
    if (data && data.results && data.results.length > 0) {
        const accepted = data.results.find(r => r.taxonomicStatus === "ACCEPTED") || data.results[0];
        return {
            kingdom: accepted.kingdom || "",
            phylum: accepted.phylum || "",
            class: accepted.class || "",
            order: accepted.order || "",
            family: accepted.family || "",
        };
    }
    return null;
}

async function fetchAllTaxonomy() {
    if (state.records.length === 0) return;

    const btn = document.getElementById("fetchTaxonomyBtn");
    btn.disabled = true;

    showProgress("Fetching Taxonomy from GBIF", `0 of ${state.records.length}`);
    document.getElementById("progressDetail").textContent = "dataFishing (Rabelo et al. 2025)";

    let successCount = 0;
    let errorCount = 0;
    let completed = 0;
    const BATCH_SIZE = 5;

    for (let i = 0; i < state.records.length; i += BATCH_SIZE) {
        const batch = state.records.slice(i, i + BATCH_SIZE);
        const promises = batch.map(async (record) => {
            const organism = record.editedOrganism || record.organism || "";

            if (!organism) {
                log(`Skipping ${record.accession}: no organism name`);
                errorCount++;
                completed++;
                updateProgress(completed, state.records.length, record.accession);
                return;
            }

            try {
                const taxonomy = await fetchTaxonomyFromGBIF(organism);
                if (taxonomy) {
                    record.taxonomyRanks = {
                        kingdom: taxonomy.kingdom,
                        phylum: taxonomy.phylum,
                        class: taxonomy.class,
                        order: taxonomy.order,
                        family: taxonomy.family,
                    };
                    if (taxonomy.family && !record.editedFamily) {
                        record.editedFamily = taxonomy.family;
                    }
                    successCount++;
                    log(`Taxonomy: ${record.accession} → ${taxonomy.kingdom}, ${taxonomy.phylum}, ${taxonomy.class}, ${taxonomy.order}, ${taxonomy.family}`);
                } else {
                    log(`No GBIF results for ${organism}`);
                    errorCount++;
                }
            } catch (e) {
                log(`Error fetching taxonomy for ${organism}: ${e.message}`);
                errorCount++;
            }

            completed++;
            updateProgress(completed, state.records.length, `<em>${organism}</em>`);
        });

        await Promise.all(promises);
        if (i + BATCH_SIZE < state.records.length) {
            await new Promise(r => setTimeout(r, 200));
        }
    }

    hideProgress(`Taxonomy fetched: ${successCount} OK${errorCount > 0 ? `, ${errorCount} failed` : ""}`);
    btn.disabled = false;
    renderRecords();
}

// ========================================================================
// Edit Record
// ========================================================================
function validateSpeciesName(name) {
    if (!name || !name.trim()) return "Species name is required";
    if (!/^[A-Za-z\s\-]+$/.test(name.trim())) return "Only letters, spaces, and hyphens are allowed";
    const parts = name.trim().split(/\s+/);
    if (parts.length < 2) return "Species name must include genus and epithet (e.g., Rhinella marina)";
    if (parts[0][0] !== parts[0][0].toUpperCase()) return "Genus must start with an uppercase letter";
    return null;
}

function validateFamilyName(name) {
    if (!name || !name.trim()) return null; // family is optional
    if (!name.trim().match(/idae$|aceae$/i)) return "Family name must end with -idae or -aceae";
    return null;
}

function formatFileSize(bytes) {
    if (bytes < 1024) return bytes + " B";
    if (bytes < 1024 * 1024) return (bytes / 1024).toFixed(1) + " KB";
    return (bytes / (1024 * 1024)).toFixed(1) + " MB";
}

window.editRecord = function (index) {
    state.pendingEditIndex = index;
    const r = state.records[index];
    const tax = extractTaxonomy(r);
    const counts = countFeatureTypes(r);

    // Build info section
    let info = '<div class="space-y-1">';
    info += `<div><span class="font-medium">Accession:</span> ${r.accession}</div>`;
    info += `<div><span class="font-medium">Source:</span> ${r.source === "ncbi" ? "NCBI" : "File"}</div>`;
    if (r.fileName) info += `<div><span class="font-medium">File:</span> ${r.fileName}</div>`;
    if (r.fileSize) info += `<div><span class="font-medium">Size:</span> ${formatFileSize(r.fileSize)}</div>`;
    if (r.fileDate) info += `<div><span class="font-medium">Date:</span> ${r.fileDate.toLocaleDateString()}</div>`;
    info += `<div><span class="font-medium">Sequence:</span> ${r.sequence ? r.sequence.length.toLocaleString() + " bp" : "N/A"}</div>`;
    info += `<div><span class="font-medium">Features:</span> ${counts.pcgs} PCGs, ${counts.rrnas} rRNAs, ${counts.trnas} tRNAs</div>`;
    info += "</div>";
    document.getElementById("editModalInfo").innerHTML = info;

    // Populate editable fields
    const currentOrganism = r.editedOrganism || r.organism || "";
    document.getElementById("editSpeciesInput").value = currentOrganism;
    document.getElementById("editGenusDisplay").value = currentOrganism.split(/\s+/)[0] || "";
    document.getElementById("editKingdomInput").value = tax.kingdom || "";
    document.getElementById("editPhylumInput").value = tax.phylum || "";
    document.getElementById("editClassInput").value = tax.class || "";
    document.getElementById("editOrderInput").value = tax.order || "";
    document.getElementById("editFamilyInput").value = tax.family || "";

    // Clear errors/status
    document.getElementById("editSpeciesError").classList.add("hidden");
    document.getElementById("editFamilyError").classList.add("hidden");
    document.getElementById("editGbifStatus").classList.add("hidden");

    // Set GBIF button state based on current species
    const speciesErr = validateSpeciesName(currentOrganism);
    document.getElementById("editGbifBtn").disabled = !!speciesErr;

    document.getElementById("editModal").classList.remove("hidden");
};

// ========================================================================
// Header Builder
// ========================================================================
const HEADER_FIELDS = [
    { key: "species", label: "Species", getter: (r, tax) => (r.editedOrganism || r.organism || "Unknown").replace(/\s+/g, "_") },
    { key: "accession", label: "Accession", getter: (r, tax) => r.accession || "NA" },
    { key: "kingdom", label: "Kingdom", getter: (r, tax) => tax.kingdom || "NA" },
    { key: "phylum", label: "Phylum", getter: (r, tax) => tax.phylum || "NA" },
    { key: "class", label: "Class", getter: (r, tax) => tax.class || "NA" },
    { key: "order", label: "Order", getter: (r, tax) => tax.order || "NA" },
    { key: "family", label: "Family", getter: (r, tax) => tax.family || "NA" },
    { key: "genus", label: "Genus", getter: (r, tax) => tax.genus || "NA" },
];

function getHeaderSeparator() {
    const checked = document.querySelector('input[name="headerSeparator"]:checked');
    return checked ? checked.value : "_";
}

function buildHeader(record, tax) {
    const sep = getHeaderSeparator();
    return state.headerTemplate
        .map(key => {
            const field = HEADER_FIELDS.find(f => f.key === key);
            return field ? field.getter(record, tax) : key;
        })
        .join(sep);
}

function checkHeaderDuplicates() {
    const headers = new Map();
    for (const record of state.records) {
        const tax = extractTaxonomy(record);
        const header = buildHeader(record, tax);
        const existing = headers.get(header) || [];
        existing.push(record.accession);
        headers.set(header, existing);
    }
    const duplicates = [];
    for (const [header, accessions] of headers) {
        if (accessions.length > 1) {
            duplicates.push({ header, accessions });
        }
    }
    return duplicates;
}

function renderHeaderBuilder() {
    const availableContainer = document.getElementById("headerAvailableFields");
    const selectedContainer = document.getElementById("headerSelectedFields");
    const previewEl = document.getElementById("headerPreview");
    const warningEl = document.getElementById("headerDuplicateWarning");
    const warningTextEl = document.getElementById("headerDuplicateText");
    const confirmBtn = document.getElementById("headerModalConfirm");

    // Available fields = all fields NOT in template
    const available = HEADER_FIELDS.filter(f => !state.headerTemplate.includes(f.key));
    availableContainer.innerHTML = available.map(f =>
        `<span class="header-chip available" onclick="addHeaderField('${f.key}')">${f.label}</span>`
    ).join("");

    // Selected fields = fields in template order
    selectedContainer.innerHTML = state.headerTemplate.map((key, idx) => {
        const field = HEADER_FIELDS.find(f => f.key === key);
        const label = field ? field.label : key;
        return `<span class="header-chip selected" draggable="true" data-idx="${idx}" onclick="removeHeaderField(${idx})">${label} <i class="fa-solid fa-xmark text-xs opacity-70"></i></span>`;
    }).join("");

    // Live preview from first record
    const sep = getHeaderSeparator();
    if (state.records.length > 0 && state.headerTemplate.length > 0) {
        const r = state.records[0];
        const tax = extractTaxonomy(r);
        previewEl.textContent = ">" + buildHeader(r, tax);
    } else if (state.headerTemplate.length === 0) {
        previewEl.textContent = "(no fields selected)";
    } else {
        previewEl.innerHTML = "&gt;Species" + sep + "Accession";
    }

    // Duplicate check
    const duplicates = checkHeaderDuplicates();
    if (duplicates.length > 0 && state.headerTemplate.length > 0) {
        const totalDupes = duplicates.reduce((sum, d) => sum + d.accessions.length, 0);
        warningTextEl.textContent = `${totalDupes} records produce ${duplicates.length} duplicated header(s). Add more fields (e.g., Accession) to make headers unique.`;
        warningEl.classList.remove("hidden");
        confirmBtn.disabled = true;
    } else if (state.headerTemplate.length === 0) {
        warningEl.classList.add("hidden");
        confirmBtn.disabled = true;
    } else {
        warningEl.classList.add("hidden");
        confirmBtn.disabled = false;
    }

    // Init drag-drop for selected chips
    initHeaderDragDrop();
}

function addHeaderField(key) {
    if (!state.headerTemplate.includes(key)) {
        state.headerTemplate.push(key);
        renderHeaderBuilder();
    }
}

function removeHeaderField(idx) {
    state.headerTemplate.splice(idx, 1);
    renderHeaderBuilder();
}

window.addHeaderField = addHeaderField;
window.removeHeaderField = removeHeaderField;

function initHeaderDragDrop() {
    const container = document.getElementById("headerSelectedFields");
    const chips = container.querySelectorAll(".header-chip.selected");
    let dragIdx = null;

    chips.forEach(chip => {
        chip.addEventListener("dragstart", (e) => {
            dragIdx = parseInt(chip.dataset.idx);
            chip.classList.add("dragging");
            e.dataTransfer.effectAllowed = "move";
        });

        chip.addEventListener("dragend", () => {
            chip.classList.remove("dragging");
            dragIdx = null;
        });

        chip.addEventListener("dragover", (e) => {
            e.preventDefault();
            e.dataTransfer.dropEffect = "move";
        });

        chip.addEventListener("drop", (e) => {
            e.preventDefault();
            const dropIdx = parseInt(chip.dataset.idx);
            if (dragIdx !== null && dragIdx !== dropIdx) {
                const item = state.headerTemplate.splice(dragIdx, 1)[0];
                state.headerTemplate.splice(dropIdx, 0, item);
                renderHeaderBuilder();
            }
        });
    });
}

function openHeaderBuilder() {
    renderHeaderBuilder();
    document.getElementById("headerModal").classList.remove("hidden");
}

// ========================================================================
// FASTA Generation
// ========================================================================
function generateFastaFiles() {
    const fastaMap = new Map();
    const dataType = state.detectedDataType;

    for (const record of state.records) {
        const bestPerGene = new Map();
        const assignedTrnas = new Set();
        const tax = extractTaxonomy(record);
        const header = ">" + buildHeader(record, tax);

        for (const feat of record.features) {
            if (!state.selectedFeatureTypes.has(feat.type)) continue;

            let geneName = null;
            if (feat.type === "tRNA") {
                geneName = standardizeTrnaName(feat);
                if (geneName && assignedTrnas.has(geneName)) {
                    if (geneName === "tRNA-Ser1") geneName = "tRNA-Ser2";
                    else if (geneName === "tRNA-Leu1") geneName = "tRNA-Leu2";
                }
                if (geneName) assignedTrnas.add(geneName);
            } else {
                const rawGene = feat.qualifiers.gene || null;
                const rawProduct = feat.qualifiers.product || null;
                geneName = standardizeGeneName(rawGene, rawProduct, dataType);
            }

            if (!geneName || !state.selectedGenes.has(geneName)) continue;

            const seq = extractSequence(record.sequence, feat.locationStr);
            if (!seq) continue;

            const existing = bestPerGene.get(geneName);
            if (!existing || seq.length > existing.len) {
                bestPerGene.set(geneName, {
                    header: header,
                    seq: seq,
                    len: seq.length,
                });
            }
        }

        for (const [geneName, { header, seq }] of bestPerGene) {
            const existing = fastaMap.get(geneName) || "";
            fastaMap.set(geneName, existing + header + "\n" + seq + "\n");
        }
    }

    return fastaMap;
}

function downloadIndividualFasta() {
    const files = generateFastaFiles();
    if (files.size === 0) {
        alert("No sequences to download. Select genes and ensure records are loaded.");
        return;
    }
    for (const [name, content] of files) {
        const blob = new Blob([content], { type: "text/plain" });
        saveAs(blob, `${name}.fasta`);
    }
    log(`Downloaded ${files.size} individual FASTA files`);
}

async function downloadZip() {
    const files = generateFastaFiles();
    if (files.size === 0) {
        alert("No sequences to download. Select genes and ensure records are loaded.");
        return;
    }
    const zip = new JSZip();
    for (const [name, content] of files) {
        zip.file(`${name}.fasta`, content);
    }
    const blob = await zip.generateAsync({ type: "blob" });
    saveAs(blob, `splace_${state.detectedDataType || "genes"}_genes.zip`);
    log(`Downloaded ZIP with ${files.size} FASTA files`);
}

// ========================================================================
// Logging
// ========================================================================
function log(msg) {
    const section = document.getElementById("logSection");
    const output = document.getElementById("logOutput");
    section.classList.remove("hidden");
    const time = new Date().toLocaleTimeString();
    output.textContent += `[${time}] ${msg}\n`;
    output.scrollTop = output.scrollHeight;
}

// ========================================================================
// File Handling
// ========================================================================
const GENBANK_EXTENSIONS = /\.(gb|gbk|genbank)$/i;

function handleFiles(fileList) {
    const validFiles = [];
    for (const file of fileList) {
        if (!GENBANK_EXTENSIONS.test(file.name)) {
            log(`Skipping non-GenBank file: ${file.name}`);
            continue;
        }
        validFiles.push(file);
    }

    if (validFiles.length === 0) return;

    showProgress("Loading files", `0 of ${validFiles.length}`);
    let loaded = 0;

    for (const file of validFiles) {
        const reader = new FileReader();
        reader.onload = (e) => {
            const text = e.target.result;
            try {
                const record = parseGenBank(text);
                record.source = "file";
                record.fileName = file.name;
                record.fileSize = file.size;
                record.fileDate = new Date(file.lastModified);
                if (!record.accession) record.accession = file.name.replace(GENBANK_EXTENSIONS, "");
                state.records.push(record);
                log(`Loaded: ${record.accession} - ${record.organism} (${record.features.length} features)`);
            } catch (err) {
                log(`Error parsing ${file.name}: ${err.message}`);
            }

            loaded++;
            updateProgress(loaded, validFiles.length, file.name);

            if (loaded === validFiles.length) {
                renderRecords();
                hideProgress(`Loaded ${validFiles.length} file${validFiles.length !== 1 ? 's' : ''}`);
            }
        };
        reader.readAsText(file);
    }
}

function readEntryAsFile(fileEntry) {
    return new Promise((resolve, reject) => {
        fileEntry.file(resolve, reject);
    });
}

function readDirectoryEntries(dirReader) {
    return new Promise((resolve, reject) => {
        dirReader.readEntries(resolve, reject);
    });
}

async function traverseEntry(entry) {
    const files = [];
    if (entry.isFile) {
        if (GENBANK_EXTENSIONS.test(entry.name)) {
            try {
                const file = await readEntryAsFile(entry);
                files.push(file);
            } catch (e) {
                log(`Error reading file ${entry.name}: ${e.message}`);
            }
        }
    } else if (entry.isDirectory) {
        const reader = entry.createReader();
        let batch;
        do {
            batch = await readDirectoryEntries(reader);
            for (const child of batch) {
                const childFiles = await traverseEntry(child);
                files.push(...childFiles);
            }
        } while (batch.length > 0);
    }
    return files;
}

async function handleDroppedItems(dataTransfer) {
    const items = dataTransfer.items;
    const allFiles = [];

    if (items) {
        const entries = [];
        for (let i = 0; i < items.length; i++) {
            const entry = items[i].webkitGetAsEntry ? items[i].webkitGetAsEntry() : null;
            if (entry) {
                entries.push(entry);
            }
        }

        if (entries.length > 0) {
            showProgress("Scanning folders", "Looking for GenBank files...");
            log("Scanning dropped items for GenBank files...");
            for (const entry of entries) {
                const files = await traverseEntry(entry);
                allFiles.push(...files);
            }
            if (allFiles.length > 0) {
                log(`Found ${allFiles.length} GenBank file(s)`);
                document.getElementById("progressModal").classList.add("hidden");
                handleFiles(allFiles);
            } else {
                log("No GenBank files found in dropped items");
                hideProgress("No GenBank files found");
            }
            return;
        }
    }

    // Fallback: plain file drop
    handleFiles(dataTransfer.files);
}

// ========================================================================
// Event Listeners
// ========================================================================
document.addEventListener("DOMContentLoaded", () => {
    const dropZone = document.getElementById("dropZone");
    const fileInput = document.getElementById("fileInput");

    dropZone.addEventListener("click", () => fileInput.click());
    dropZone.addEventListener("dragover", (e) => { e.preventDefault(); dropZone.classList.add("drop-active"); });
    dropZone.addEventListener("dragleave", () => dropZone.classList.remove("drop-active"));
    dropZone.addEventListener("drop", (e) => {
        e.preventDefault();
        dropZone.classList.remove("drop-active");
        handleDroppedItems(e.dataTransfer);
    });
    fileInput.addEventListener("change", (e) => handleFiles(e.target.files));

    document.getElementById("exampleBtn").addEventListener("click", fillExampleAccessions);

    document.getElementById("fetchBtn").addEventListener("click", async () => {
        const raw = document.getElementById("accessionInput").value;
        const accessions = raw.split(/[\s,;\n]+/).filter(a => a.trim());
        if (accessions.length === 0) return;

        const btn = document.getElementById("fetchBtn");
        btn.disabled = true;
        btn.innerHTML = '<i class="fa-solid fa-spinner fa-spin"></i> Fetching...';
        try {
            const records = await fetchMultiple(accessions);
            state.records.push(...records);
            renderRecords();
        } catch (e) {
            log(`Fetch error: ${e.message}`);
        }
        btn.disabled = false;
        btn.innerHTML = '<i class="fa-solid fa-magnifying-glass"></i> Fetch';
    });

    document.getElementById("clearRecords").addEventListener("click", () => {
        document.getElementById("clearModal").classList.remove("hidden");
    });

    document.getElementById("clearModalCancel").addEventListener("click", () => {
        document.getElementById("clearModal").classList.add("hidden");
    });

    document.getElementById("clearModalOverlay").addEventListener("click", () => {
        document.getElementById("clearModal").classList.add("hidden");
    });

    document.getElementById("clearModalConfirm").addEventListener("click", () => {
        document.getElementById("clearModal").classList.add("hidden");
        state.records = [];
        state.selectedGenes.clear();
        renderRecords();
    });

    document.getElementById("removeModalCancel").addEventListener("click", () => {
        document.getElementById("removeModal").classList.add("hidden");
    });

    document.getElementById("removeModalOverlay").addEventListener("click", () => {
        document.getElementById("removeModal").classList.add("hidden");
    });

    document.getElementById("removeModalConfirm").addEventListener("click", () => {
        document.getElementById("removeModal").classList.add("hidden");
        if (state.pendingRemoveIndex !== null && state.pendingRemoveIndex !== undefined) {
            state.records.splice(state.pendingRemoveIndex, 1);
            state.pendingRemoveIndex = null;
            renderRecords();
        }
    });

    // Edit modal listeners
    document.getElementById("editModalCancel").addEventListener("click", () => {
        document.getElementById("editModal").classList.add("hidden");
    });

    document.getElementById("editModalOverlay").addEventListener("click", () => {
        document.getElementById("editModal").classList.add("hidden");
    });

    document.getElementById("editSpeciesInput").addEventListener("input", function () {
        const parts = this.value.trim().split(/\s+/);
        document.getElementById("editGenusDisplay").value = parts[0] || "";
        // Enable/disable GBIF button based on species validation
        const err = validateSpeciesName(this.value);
        document.getElementById("editGbifBtn").disabled = !!err;
    });

    document.getElementById("editGbifBtn").addEventListener("click", async function () {
        const species = document.getElementById("editSpeciesInput").value.trim();
        if (!species) return;
        const statusEl = document.getElementById("editGbifStatus");
        statusEl.textContent = "Searching GBIF...";
        statusEl.className = "text-xs text-gray-500 mt-1";
        statusEl.classList.remove("hidden");
        this.disabled = true;
        try {
            const taxonomy = await fetchTaxonomyFromGBIF(species);
            if (taxonomy) {
                document.getElementById("editKingdomInput").value = taxonomy.kingdom || "";
                document.getElementById("editPhylumInput").value = taxonomy.phylum || "";
                document.getElementById("editClassInput").value = taxonomy.class || "";
                document.getElementById("editOrderInput").value = taxonomy.order || "";
                document.getElementById("editFamilyInput").value = taxonomy.family || "";
                statusEl.innerHTML = `Found via <em>dataFishing</em> (Rabelo et al. 2025)`;
                statusEl.className = "text-xs text-green-600 mt-1";
                document.getElementById("editFamilyError").classList.add("hidden");
            } else {
                statusEl.textContent = "No results found in GBIF";
                statusEl.className = "text-xs text-orange-500 mt-1";
            }
        } catch (e) {
            statusEl.textContent = `Error: ${e.message}`;
            statusEl.className = "text-xs text-red-500 mt-1";
        }
        this.disabled = false;
    });

    document.getElementById("editModalSave").addEventListener("click", () => {
        const species = document.getElementById("editSpeciesInput").value.trim();
        const family = document.getElementById("editFamilyInput").value.trim();
        const kingdom = document.getElementById("editKingdomInput").value.trim();
        const phylum = document.getElementById("editPhylumInput").value.trim();
        const klass = document.getElementById("editClassInput").value.trim();
        const order = document.getElementById("editOrderInput").value.trim();

        const speciesError = validateSpeciesName(species);
        const familyError = validateFamilyName(family);

        const speciesErrEl = document.getElementById("editSpeciesError");
        const familyErrEl = document.getElementById("editFamilyError");

        if (speciesError) {
            speciesErrEl.textContent = speciesError;
            speciesErrEl.classList.remove("hidden");
        } else {
            speciesErrEl.classList.add("hidden");
        }

        if (familyError) {
            familyErrEl.textContent = familyError;
            familyErrEl.classList.remove("hidden");
        } else {
            familyErrEl.classList.add("hidden");
        }

        if (speciesError || familyError) return;

        const r = state.records[state.pendingEditIndex];
        r.editedOrganism = species;
        r.organism = species;
        if (family) r.editedFamily = family;

        // Save taxonomy ranks
        r.taxonomyRanks = {
            kingdom: kingdom,
            phylum: phylum,
            class: klass,
            order: order,
            family: family || (r.taxonomyRanks && r.taxonomyRanks.family) || "",
        };

        document.getElementById("editModal").classList.add("hidden");
        renderRecords();
    });

    document.getElementById("selectAllGenes").addEventListener("click", selectAllGenes);
    document.getElementById("selectNoneGenes").addEventListener("click", selectNoneGenes);
    document.getElementById("selectDefaultGenes").addEventListener("click", selectDefaultGenes);
    document.getElementById("downloadBtn").addEventListener("click", () => openHeaderBuilder());

    // Fetch Taxonomy button
    document.getElementById("fetchTaxonomyBtn").addEventListener("click", fetchAllTaxonomy);

    // Column visibility eye toggles
    document.querySelectorAll(".col-eye-toggle[data-col-toggle]").forEach(btn => {
        btn.addEventListener("click", () => {
            const col = btn.dataset.colToggle;
            const isHidden = state.hiddenColumns.has(col);
            toggleColumn(col, isHidden);
        });
    });

    // Initialize column visibility checkboxes
    applyColumnVisibility();

    // Header builder modal listeners
    document.getElementById("headerModalCancel").addEventListener("click", () => {
        document.getElementById("headerModal").classList.add("hidden");
    });

    document.getElementById("headerModalOverlay").addEventListener("click", () => {
        document.getElementById("headerModal").classList.add("hidden");
    });

    document.getElementById("headerModalConfirm").addEventListener("click", () => {
        document.getElementById("headerModal").classList.add("hidden");
        const format = document.querySelector('input[name="downloadFormat"]:checked').value;
        if (format === "individual") {
            downloadIndividualFasta();
        } else {
            downloadZip();
        }
    });

    // Re-render header builder when separator changes
    document.querySelectorAll('input[name="headerSeparator"]').forEach(radio => {
        radio.addEventListener("change", () => renderHeaderBuilder());
    });
});
