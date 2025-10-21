# 3D_FCC_Track_scoring
Code used for the analysis presented in Radiat. Env. Biophys. 60, 559-578 (2021) doi 10.1007/s00411-021-00936-4 and 10.48550/arxiv.2105.07159

**Source code**  
*last change*  
Description  

**IC_3D.f**  
*02.04.2021 17:08*  
Determines ionization clusters in particle tracks and saves the cluster positions (and their complexity) to an output file (prefix IC_).  
Uses subroutines in IC_3D_SUBS.f and ROI_3D_Init.f.   

**IC_3D_SUBS.f**  
*02.04.2021 17:08*   
Subroutine CLUSTR called from IC_3D after each track has been read and scores ICs in Wigner Seitz cells and calculates centers of gravity of the clusters 
 
**ROI_3D.f**	 
*17.04.2021 12:37*  
Reads track data (output files of IC_3D or original simulation data) and calculates and outputs   
–	frequency distributions of the number of Wigner-Seitz cells in a (large) spherical region that receive ionization clusters for track at different impact parameters and “true” (infinite radial integral) and conditional (track intersects the spherical target) single event distributions  (output file name 3D_*.dat)  
–	bivariate distributions of Wigner Seitz cells containing single or multiple ionization clusters (Output file name 3B_*.dat)  
–	ratio of bivariate frequency distribution of Wigner Seitz cells containing single or multiple ionization clusters to product of marginal frequencies (Output file name 3C_*.dat)    

**ROI_3D_Init.f** 	 
*31.03.2021 13:01*  
Initializes the Bravais lattice used for scoring  

**ROI_3D_SUBS.f**  
*17.04.2021 12:35*  
Subroutine TARG3D called from ROI_3D after each track has been read and scores ICs in Wigner Seitz cells 

ME_ROI_3C.f	17.04.2021 13:07	Reads data from output files 3B_*.dat produced by ROI_3D and convolutes them with Binomial distributions of a given success probability such as to convert IC to DSB distributions. Outputs:
–	MEA_ multi- and single event frequency distributions after convolution with binomial compared to frequency distributions of ionization clusters
–	MEB_ bivariate multi-event frequency distribution of single and multiple ionization clusters
–	MEC_ ratio of bivariate frequency distribution to product of marginal frequencies
–	MED_ multi- and single event frequency distributions after convolution with binomial only (smaller file size)
–	SEB_ bivariate single-event frequency distribution of single and multiple ionization clusters
–	SEC_ ratio of bivariate frequency distribution to product of marginal frequencies
ME_ROI_3D.f	12.04.2021 17:44	Reads output files 3D_*.dat from ROI_3D, calculates multi event distributions for targets with ionization clusters and produces several output files for further use of the results, named by adding a prefix to the input files:
–	ME_ multi-event distributions
–	SE_ single-event distributions
–	TE_ conditional single event distributions (traversal)
ALL_ all of the three above plus single tracks at certain distances plus tracks passing through annuli around region of interest cross section
ME_ROI_3P.f	15.04.2021 17:40	Uses the bivariate distributions output from ROI_3D (3B_*.dat) and convolutes them with Binomial distributions of a given success probability such as to convert IC to DSB distributions.
This is a faster version of ME_ROI_3C, as it does not process the bivariate distributions.
Output:
–	MEP_ multi- and single event frequency distributions after convolution with binomial compared to frequency distributions of ionization clusters
SE_ROI_3C.f	15.04.2021 12:28	Uses the bivariate distributions output from ROI_3D and convolutes them with Binomial distributions of a given success probability and the convolutes the results 1024 times with itself such as to estimate results for a cell nucleus from results found in 500 nm spheres (from Sonwabile’s data).
–	SEB_ bivariate single-event frequency distribution of single and multiple ionization clusters
–	SEC_ ratio of bivariate frequency distribution to product of marginal frequencies
ROI_2D.f	02.04.2021 17:16	Same as ROI_3D, but for a large cylindrical target (tracks parallel to cylinder axis)
ROI_2D_SUBS.f	02.04.2021 16:12	Subroutine TARG2D called from ROI_2D after each track has been read and scores ICs in Wigner Seitz cells
