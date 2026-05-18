## Code to replicate analyses and figures of the manuscript "Big heads and snapping jaws: Skull shape of the extinct Tasmanian tiger suggests a unique biting style"
Manuscript authors: Vera Weisbecker, Andrew J Pask, Axel H Newton, and Douglass S Rovinsky

Code authors: Douglass S Rovinsky, Vera Weisbecker.

This project is licensed under the terms of the [GNU General Public License, version 3.0 (GPL-3.0)](https://web.archive.org/web/20160412212723/https://opensource.org/licenses/gpl-3.0).

*All scripts are in RMarkdown format (.Rmd), which can be opened in RStudio. There, you can edit and run code chunks as normal or use the Knit button to create HTML versions with both code and output. After cloning this repo, remember to either set your working directory to the correct folder on your computer or open an RStudio project from that folder.*

## [Data](/Data)

Contains all data required to run the analyses, including the [landmark data](/Data/BWT_ABN_lmks.txt) with associated [landmark pairings](/Data/LmkPairs_v2.csv), [development-based partitions](/Data/Sets_Newtonetal2021.csv), and [phylogeny file](/Data/phyloGMM.tre); [body mass estimates for species](/Data/Body_mass.csv); [linear cranial morphometric data](/Data/Thylacine_Canis.csv); and [infraorbital foramen dataset collation](/Data/IOF_Data_Muchlinski.csv).


    
## [Analyses](/Analyses)
**These scripts save several outputs from raw data for downstream use**. These intermediate data are stored as .rda files in the [Data/Processed](/Data/Processed) folder. Figure (and table) outputs are stored in the [Figures](/Figures) These are too large to upload to GitHub so scripts must be run sequentially to generate these outputs. 

* [**01_1_GMM_GPA.Rmd**](/Analyses/01_1_GMM_GPA.Rmd) extracts  3D coordinates from the original Viewbox file and prepares them for analysis. Runs GPA with bilateral symmetry for the whole cranium and all partitions and saves the results in .rda files in the [Data/Processed](/Data/Processed) folder.

* [**01_2_Allometry.Rmd**](Analyses/01_2_Allometry.Rmd) computes analyses of allometry for Supplementary Note 1, including Supplementary Figures 6-10 and Supplementary Tables 2-6.

* [**02_PCA_Functional.Rmd**](Analyses/02_PCA_Functional.Rmd) runs the Principal Components Analyses (PCA) for the whole-cranium and rostral _versus_ neurocranial landmark partitions for outputs in Fig. 1.

* [**03_PCA_Developmental.Rmd**](Analyses/03_PCA_Developmental.Rmd) runs the Principal Components Analyses (PCA) for the developmental partitions and creates outputs for Supplementary Fig. 2

* [**04_Procrustes_Funct_Analyses.Rmd**](Analyses/04_Procrustes_Funct_Analyses.Rmd) computes Procrustes distances and outputs for the whole-cranium and rostral _versus_ neurocranial landmark partitions for Fig. 2a

* [**05_Procrustes_Devo_Analyses.Rmd**](Analyses/05_Procrustes_Devo_Analyses.Rmd) computes Procrustes distances and outputs for the developmental landmark partitions for Fig. 2b and Supplementary Table 1 (combined with distances from previous file).

* [**06_Csize_Analyses.Rmd**](Analyses/06_Csize_Analyses.Rmd) produces body mass _versus_ centroid size plots for Fig. 3

* [**07_Cranial_Linear_Metrics.Rmd**](Analyses/07_Cranial_Linear_Metrics.Rmd) computes regression comparisons between linear measurements _versus_ the geometric mean of all measurements for Fig. 4 and Supp. Table 7.

* [**08_Terminal_Rosette_Analysis.Rmd**](Analyses/08_Terminal_Rosette_Analysis.Rmd) produces the outputs for the "rostral pinch" Fig. 5 and measurement protocol illustration for Supp. Fig. 14

* [**09_IOF.Rmd**](Analyses/09_IOF.Rmd) produces the outputs of Fig. 6 and Supp. Fig. 11 of comparative size of the infraorbital foramen. 

* [**10_Weisbecker_etal_2023_Csize.Rmd**](Analyses/10_Weisbecker_etal_2023_Csize.Rmd) shows that the thylacine's cranial size is an outlier among living marsupials (Supp. Fig. 5 )

*  [**11_Landmark_Config_Vis.Rmd**](Analyses/11_Landmark_Config_Vis.Rmd) creates images of landmark configurations used in Fig. 2 and 3
  
* [**12_Retrodeform.Rmd**](Analyses/12_Retrodeform.Rmd) includes the code used for retrodeformation of fossils as described in Materials and Methods
  
* [**13_phyloMap_Tree_Figures.Rmd**](Analyses/13_phyloMap_Tree_Figures.Rmd) code for the heat maps of distances from the thylacine for Supp. Figs. 3 and 4.

* [**14_Supp_Note_2_Wheatsheaf_Index.Rmd**](Analyses/14_Supp_Note_2_Wheatsheaf_Index.Rmd) Code for replicating the supplementary exploration of the Wheatsheaf index in Supplementary Note 2

### Custom functions
The Benjamini-Hochberg procedure and use of the Evomap package call custom functions, most of which are defined in the [..Data/Functions/Fun.R](/Data/Functions/Fun.R) file. 


