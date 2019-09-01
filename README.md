# GISA
C code and more for the GISA algorithm. All related to paper

"C.Grønbæk., T.Hamelryck, P.Røgen: "GISA: Using Gauss integrals to identify rare conformations in protein structures", TO APPEAR."

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
