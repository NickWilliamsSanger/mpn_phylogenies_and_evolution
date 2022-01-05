# mpn_phylogenies_and_evolution

Code and documents for "Life histories of myeloproliferative neoplasm inferred from phylogenies"
    
## Directories

### supplementary_reports

Contains reports supporting the supplementary notes.  The front pages of the QC PDF files  (notes 1-4) contain legends for the figures therein.

* supplementary_note_1_patients_colony_qc_first_pass.pdf : Per patient initial QC assessment plots prior to removal of low quality colonies
* supplementary_note_1_patients_colony_qc_final_pass.pdf : Per patient final QC assessment plots after quality control and removal of low quality colonies
* supplementary_note_3a_tree_qc.pdf : Tree QC plots following removal of low quality colonies
* supplementary_note_3b_per_patient/<PDID>_by_colony.pdf : Tree By Colony Quality Assessment
* supplementary_note_3c_normal_contamination_in_trees.pdf  : Tumour-In-Normal Per Branch Assessment
* supplementary_note_5_per_patient_reports : Per patient reports <PDID>.html showing
    * Timing of driver branches    
    * Timing of chromosomal duplication Events and copy number neutral LOH.
    * Per branch VAF of targeted follow up samples.
* supplementary_note_6_signatures.html : Signature analysis
* supplementary_note_7_telomeres.html : Telomere analysis
* supplementary_note_8_ABCExampleRun.html : Example showing how ABC simulations can be generated using rsimpop ( https://github.com/NickWilliamsSanger/rsimpop )



### example

Reproducible analysis for PD6629 starting from preprocessed files containing read counts of somatic mutations:
* Timing of driver branches
* Timing of chromosomal duplication Events and copy number neutral LOH.
* Per branch VAF of targeted follow up samples.
* Additional tree QC plots.
See example/README.txt for instructions to run the analysis.

### data

Files required by telomere and ABC markdown.

## Tree of contents

```text
.
├── README.md
├── data
│   ├── PDD_TELO.RDS
│   └── example_treedat.RDS
├── example
│   ├── 1805_PD6629aj_ascat_ngs.summary.csv
│   ├── GENES.txt
│   ├── MUTCOUNTBINS.RDS
│   ├── PD6629.RDS
│   ├── PlotVafTreeBig.Rmd
│   ├── README.txt
│   ├── SingleSampleAnalysisExample.Rmd
│   ├── SingleSampleAnalysisExample.html
│   ├── SingleSampleAnalysisExample_RERUN.htm
│   ├── cna.R
│   ├── code.R
│   ├── consequence_col_scheme.txt
│   ├── driver_groups.txt
│   ├── driver_scheme.txt
│   ├── drivers_red.bed
│   ├── plot_tree_annots_extra.R
│   └── stanmodels
│       ├── cntime.stan
│       └── cntimedup.stan
└── supplementary_reports
    ├── supplementary_note_1_patients_colony_qc_final_pass.pdf
    ├── supplementary_note_1_patients_colony_qc_first_pass.pdf
    ├── supplementary_note_3a_tree_qc.pdf
    ├── supplementary_note_3b_per_patient
    │   ├── PD4781_by_colony.pdf
    │   ├── PD5117_by_colony.pdf
    │   ├── PD5147_by_colony.pdf
    │   ├── PD5163_by_colony.pdf
    │   ├── PD5179_by_colony.pdf
    │   ├── PD5182_by_colony.pdf
    │   ├── PD5847_by_colony.pdf
    │   ├── PD6629_by_colony.pdf
    │   ├── PD6634_by_colony.pdf
    │   ├── PD6646_by_colony.pdf
    │   ├── PD7271_by_colony.pdf
    │   └── PD9478_by_colony.pdf
    ├── supplementary_note_3c_normal_contamination_in_trees.pdf
    ├── supplementary_note_5_per_patient_reports
    │   ├── PD4781.html
    │   ├── PD5117.html
    │   ├── PD5147.html
    │   ├── PD5163.html
    │   ├── PD5179.html
    │   ├── PD5182.html
    │   ├── PD5847.html
    │   ├── PD6629.html
    │   ├── PD6634.html
    │   ├── PD6646.html
    │   ├── PD7271.html
    │   └── PD9478.html
    ├── supplementary_note_6_signatures.html
    ├── supplementary_note_7_telomeres.Rmd
    ├── supplementary_note_7_telomeres.html
    ├── supplementary_note_8_ABCExampleRun.Rmd
    └── supplementary_note_8_ABCExampleRun.html
```