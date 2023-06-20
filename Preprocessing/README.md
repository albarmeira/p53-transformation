# **Single-Cell Multi-Omic Analysis of _TP53_-mutant leukemic transformation (Preprocessing)**

This folder contains the scripts and data involved in the preprocessing of single-cell dataset generated in Rodriguez-Meira, Norfo, Wen, Chedeville *et al.*, 2023. The outputs files from here were used to generate the published figures. 

*   02\_HTMPNAML\_SingleCellaR.revised.github.R: Create SingCellaR object and QC using the raw gene count matrix and cell metadata. This dataset consists of preleukemic (TP53 wildtype) and leukemic (TP53 heterozygous and multi-hit) cells derived from chronic phase (CP) and acute phase (AP) MPN patients, and healthy controls. This dataset was generated and published by Rodriguez-Meira, Norfo, Wen, Chedeville *et al.*, 2023.

*   03\_MF\_MolCell_SingleCellaR.revised.github.R: Create SingCellaR object and QC using raw gene count matrix and cell metadata. This dataset consists of preleukemic (TP53 wildtype) from myelofibrosis (MF) patients (CP MPN) and healthy controls. This dataset was generated and published by Rodriguez-Meira *et al*., 2019.

*   04\_HTMPNAML\_Harmony.AP\_SinglCellaR.revised.github.R: Subset AP MPN patient cells from SingCellaR object generated from 02\_HTMPNAML\_SingleCellaR.revised.github.R, and subsequently perform integration with Harmony, dimension reduction, and louvain clustering.

*   05\_Integration\_MF\_p53MPNAML\_SinglCellaR.revised.github.R: Merge SingCellaR objects generated from 02\_HTMPNAML\_SingleCellaR.revised.github.R and 03\_MF\_MolCell_SingleCellaR.revised.github.R. Cluster assignment into preleukemic, leukemic, and erythroid based on cells from AP MPN patients was additionally performed in this step.

*   06.1_MolCell2019_CD38negative.revised.github.R: Create new raw gene count matrix and cell metadata for only CD38-negative cells from myelofibrosis (MF) patients (CP MPN) and healthy controls from Rodriguez-Meira *et al*., 2019.

*   06.2\_p53MPNAML\_preLSCCD38negative.revised.github.R: Create new raw gene count matrix and cell metadata for only CD38-negative, preleukemic (TP53 wildtype) cells from AP MPN patients from Rodriguez-Meira, Norfo, Wen, Chedeville *et al.*, 2023.

*   06.3\_preleukemic\_NDMFAML\_integration.revised.github.R: Merge raw gene count matrix and cell metadata generated from 06.1_MolCell2019_CD38negative.revised.R and 06.2\_p53MPNAML\_preLSCCD38negative.revised.R

*   06.4\_preleukemic.integration\_SinglCellaR.revised.github.R: Create SingCellaR object using the raw gene count matrix and cell metadata generated from 06.3\_preleukemic\_NDMFAML\_integration.revised.R.
