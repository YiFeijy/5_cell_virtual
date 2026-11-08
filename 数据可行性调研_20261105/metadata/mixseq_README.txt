MIX-seq Data Release 

-----------------------------------------------------------
Contents:
-----------------------------------------------------------
The dataset consists of 10x single-cell RNA-seq data from 30 different experiments, where pools of cancer cell lines were treated with different small-molecule or genetic perturbations. The dataset is described in more detail in [1]. Metadata for each experiment can be found in Supplementary Table 2 (included). 

Each data folder (compressed as a zip file) contains the following files:
	barcodes.tsv: table of identified cell barcodes
	genes.tsv: table of gene information (Ensemble ID and HGNC symbol)
	matrix.mtx: read counts matrix stored as matrix market format
	classifications.csv: table of cell metadata (see below for further information)

In addition to the 30 data folders, the dataset includes:
	- Supplementary Tables.xlsx: 3 Supplementary Tables associated with the MIX-seq paper
		- Table 1 has information about which cancer cell lines were included in which experimental pools 
		- Table 2 has metadata for each experiment
		- Table 3 has metadata for the small-molecules used

	- all_CL_features.rds: Packaged R data file containing a named list of tables, one for each treatment-vs-control comparison used in the manuscript. Each table contains cell line features relevant to that comparison (cell lines that were present in the experimental pool are indicated in the 'in_pool' column). Drug sensitivity ('sens') is computed as 1-AUC_avg, and AUC_avg is given by the average of (quantile-normalized) AUC values from the PRISM [2] and GDSC [3] datasets where available.

*classifications.csv*: Data dictionary for cell information included in each classifications.csv file. See manuscript for more details on SNP-based cell classification
	- barcode: Cell barcode 
	- singlet_ID: CCLE_ID of best-matching reference cell line (using 'singlet' model).
	- num_SNPs: number of SNP sites detected (with > 0 reads)
	- singlet_dev: deviance ratio of best-fitting singlet-model
	- singlet_dev_z: z-score normalized (across reference cell line panel) of best-matching singlet-model deviance
	- singlet_margin: difference between best-matching and second-best-matching singlet-model deviance ratio
	- singlet_z_margin: same as singlet_margin, but using z-score normalized deviance ratios
	- doublet_z_margin: difference between 2nd and 3rd best-fitting singlet models (potentially to be deprecated, not used)
	- tot_reads: total number of SNP reads
	- doublet_dev_imp: difference between doublet and singlet model deviance ratios (higher number is more evidence for doublet)
	- doublet_CL1: first cell line of best-matching doublet pair
	- doublet_CL2: second cell line of best-matching doublet pair
	- percent.mito: percent of reads from mitochondrial genes
	- cell_det_rate: fraction of genes detected with >0 reads
	- cell_quality: classification for the cell. Options:
		- normal (QC-passing singlet. Typically, restrict analysis to these)
		- doublet (doublet, likely mixture of two cells sharing a barcode)
		- empty_droplet (putative empty droplet)
		- low_confidence (cell that's likely a singlet, but fails classification confidence threshold)
		- low-quality (cell failing QC criteria, based on number of SNPs, mitochondrial reads, etc)
	- doublet_GMM_prob: probability assigned to doublet vs singlet models, based on 2-component GMM fit
	- DepMap_ID: DepMap cell line ID



-----------------------------------------------------------
References:
-----------------------------------------------------------
[1] McFarland, J.M., Paolella, B.R., Warren, A. et al. Multiplexed single-cell transcriptional response profiling to define cancer vulnerabilities and therapeutic mechanism of action. Nat Commun 11, 4296 (2020). https://doi.org/10.1038/s41467-020-17440-w

[2] Corsello, S.M., Nagari, R.T., Spangler, R.D., Rossen, J., Kocak, M., Bryan, J.G., Humeidi, R., Peck, D., Wu, X., Tang, A.A., Wang, V.M., Bender, S.A., Lemire, E., Narayan, R., Montgomery, P., Ben-David, U., Chen, Y., Rees, M.G., Lyons, N.J., McFarland, J.M., Wong, B.T., Wang, L., Dumont, N., O’Hearn, P.J., Stefan, E., Doench, J.G., Greulich, H., Meyerson, M., Vazquez, F., Subramanian, A., Roth, J.A., Bittker, J.A., Boehm, J.S., Mader, C.C., Tsherniak, A. and Golub, T.R. 2019. Non-oncology drugs are a source of previously unappreciated anti-cancer activity. BioRxiv.

[3] Iorio, F., Knijnenburg, T.A., Vis, D.J., Bignell, G.R., Menden, M.P., Schubert, M., Aben, N., Gonçalves, E., Barthorpe, S., Lightfoot, H., Cokelaer, T., Greninger, P., van Dyk, E., Chang, H., de Silva, H., Heyn, H., Deng, X., Egan, R.K., Liu, Q., Mironenko, T., Mitropoulos, X., Richardson, L., Wang, J., Zhang, T., Moran, S., Sayols, S., Soleimani, M., Tamborero, D., Lopez-Bigas, N., Ross-Macdonald, P., Esteller, M., Gray, N.S., Haber, D.A., Stratton, M.R., Benes, C.H., Wessels, L.F.A., Saez-Rodriguez, J., McDermott, U. and Garnett, M.J. 2016. A landscape of pharmacogenomic interactions in cancer. Cell 166(3), pp. 740–754.


-----------------------------------------------------------
Version history:
-----------------------------------------------------------
v1: Initial data release
v2: Updated all_CL_features.rds file to include some additional data needed to run certain scripts
v3: Added one dataset that was missing (DMSO_expt10). Also changed the names of files from expt4 to expt10 to be consistent with the naming scheme used in the supplementary tables for the manuscript.
