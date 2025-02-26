## CHARMM input scripts used to perform laddered extensions for the B7 TCR simulations for the article published in https://elifesciences.org/reviewed-preprints/104280

Authors: Ana C. Chang-Gonzalez & Wonmuk Hwang ( hwm@tamu.edu ). see LICENSE.

**Primary reference:**

Ana C. Chang-Gonzalez, Aoi Akitsu, Robert J. Mallis, Matthew J. Lang, Ellis L. Reinherz, and Wonmuk Hwang

**Load-based divergence in the dynamic allostery of two TCRs recognizing the same pMHC**

*eLife* (2025).

This repository condtains input scripts for laddered extensions in the above paper. For analysis scripts, see: https://github.com/hwm2746/a6tcr_anal_md

**General instructions:**

[0] Select residues forming contacts between the TCR and pMHC to which a distance restraint were applied. We used residues from contacts which were present during heating and equilibration. These are defined in pair_h0-e0.str and resd.str. 

[1] CPT simulation with distance restraint for the selected residues 

- Inputs:
  - Protein: 1bd2.psf , ./rst/1bd2e0.rst (output from equilibration run)
  - Distance restraints: pair_h0-e0.str , resd.str (generated in step 0)   
- Output: 1bd2d0a.{rst,dcd}

charmm < dyn0a.inp > out_dyn0a.dat

[2] Laddered extensions: incremental increases to the distance between end-terminal MHC and TCR residues. Distance increase is accompanied by a 2-ns CPT simulation. 

- Inputs:
  - ./rst/1bd2d0a.rst , or restart file from previous extension (1bd2id0a_incr4.rst in dyn0a_incr8.inp)
  - Output: 1bd2d0a_incr8.rst. This is the restart file to use in production runs.

Scripts alike dyn0a_incr8.inp were used to increase or decrease extension length by "dx". At a maximum we changed the extension distance by 4 &#x212b; on each end. To get B7<sup>high</sup> which is a 16 &#x212b; increase from the baseline, we first increased the extension by 4 &#x212b; then used the restart file from that CPT simulation for the extension. Interface distance restraints (pair_h0-e0.str and resd.str) are applied.

[3] Production runs were performed with OpenMM. In folders b7low and b7high:

  (a) shift.inp shifts the periodic box so that a corner is at origin as needed for OpenMM. Generates 1bd2id0a.cor in b7low and 1bd2id0a_incr8.cor in b7high. Translation coordinates are used in dyn_ini.py and dyn_cont.py. 
  (b) get_pos.inp prints hold_atom and pull_atom, the coordinates of end C<sub>&#x03B1;</sub> atoms. These are manually added to dyn_ini.py and dyn_cont.py as variables "mhc0_x0" and "tcr_x0."
  (c) dyn_ini.py initializes OpenMM simulation using coordinates from shift.inp, velocities from restart in either step 1 or step 2. Outputs dcd and checkpoint file to continue simulation. 
  (d) dyn_cont.py carries out production runs using checkpoint file from (c).

In dyn_{ini,cont}.py, "mhc_x0" and "tcr_x0" define x-coordinate relative to which end-terminal atoms should be held. 

------------------------------
About "stream include.str" : the file include.str loads charmm parameter files. This needs to be set to your CHARMM installation directory. For example:

set chm_dir /usr/local/charmm/c47a11/
read rtf card name @{chm_dir}toppar/top_all36_prot.rtf
read para card flex name @{chm_dir}toppar/par_all36_prot.prm
stream @{chm_dir}toppar/toppar_water_ions.str ! water & ions
