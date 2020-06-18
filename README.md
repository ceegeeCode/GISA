# GISA
C code and more for the GISA algorithm. All related to paper (please cite upon use of this code)

Grønbæk C, Hamelryck T, Røgen P. 2020. GISA: using Gauss Integrals to identify rare conformations in protein structures. PeerJ 8:e9159 https://doi.org/10.7717/peerj.9159

The repository contents are:

Branch master:
1. ReadMe pdf-file elaborating on how-to-use.
2. Examples of commands for runs.
3. Compiled C code: .o files. Source C code: .c and .h files.
4. C code: the GISA_main_* and SAonGISA_main_* are wrapping the source GISA_* and SAonGISA_* for command line usage.
5. Python code and Pymol-scripts allowing to fetch out the most exceptional cases of mutual writhe values in the output from the C-code and to generate various plots (distributions of the mutual writhe values and 3d-plots of the stuctures with sub-chain pairs high-lighted).

Branch archive_v2:

Similar files for the first version of GISA (no SAonGISA files).

Branch develop:

Files in development stage.

Branch pythonPrecursor:

Python implementation of a recursion for computing the Gauss Integrals up to order 3 (incl). The code consists in three modules of which git_utils contains the basic computation of writhe of two line segments (by means of Gauss-Bonnet); gitP contains the implementation of the recursion (and much else); git_bf is just an implementation of the Gauss Integrals directly off their definitions, i.e. by brute force. The latter was used only for checking the validity of the recursion and the implementation of it. 
