
# EnzyWizard-Dock

EnzyWizard-Dock is a command-line tool for performing molecular docking
of one or multiple substrates with a cleaned protein structure and
generating a detailed JSON report.
It takes a CIF or PDB protein structure as input and performs docking
using AutoDock Vina. The tool supports both single-substrate docking and
simultaneous multi-substrate docking.
For each substrate, multiple conformations (protomers) can be provided as
separate SDF files. The program automatically enumerates all possible
combinations of conformations of the same substrate and performs docking to identify
the optimal binding result.
Binding pockets are detected using PyVOL with a global docking box as fallback.
The tool outputs structured docking results suitable for downstream
applications such as enzyme-substrate interaction analysis, binding mode
evaluation, and enzyme characterising.


# example usage:

Example command:

enzywizard-dock -i examples/input/cleaned_3GP6.cif -s glucose,fructose -d examples/input/ -o examples/output/


# input parameters:

-i, --input_path
Required.
Path to the input cleaned protein structure file in CIF or PDB format.

-s, --substrate_names
Required.
Input substrate names separated by ','.

Examples:
- glucose
- glucose,fructose
- smiles1,smiles2

This parameter represents a multi-substrate combination for docking.

For input substrate_names "SubstrateA,SubstrateB", the program searches
substrate_dir for matched SDF files:

- SubstrateA.sdf
- SubstrateA_1.sdf
- SubstrateA_2.sdf
- ...

- SubstrateB.sdf
- SubstrateB_1.sdf
- SubstrateB_2.sdf
- ...

Each matched SDF file represents a different conformation.
The program automatically enumerates combinations such as:

- SubstrateA_1 + SubstrateB_1
- SubstrateA_1 + SubstrateB_2
- ...

and performs docking for each combination to identify the optimal result.

Duplicate substrate names are not allowed.

-d, --substrate_dir
Required.
Path to a directory containing input substrate SDF files.

-o, --output_dir
Required.
Directory to save docking outputs and the JSON report.

--max_docking_attempt_num
Optional.
Maximum number of docking attempts.

Default:
  20

--early_stop
Optional.
Whether to stop immediately after the first successful docking result.

Default:
  False

If True, the program returns the first successful result and does not continue searching for a better one.

--exhaustiveness
Optional.
Exhaustiveness of AutoDock Vina search.

Default:
  16

Larger values may improve docking search coverage but increase runtime.

--cpu
Optional.
Number of CPUs used by AutoDock Vina.

Default:
  0

A value of 0 lets Vina decide automatically.

--min_rad
Optional.
Minimum probe radius used in pocket detection.

Default:
  1.8

--max_rad
Optional.
Maximum probe radius used in pocket detection.

Default:
  6.2

--min_volume
Optional.
Minimum pocket volume threshold.

Default:
  50


# output content:

The program outputs the following files into the output directory:

1. A JSON report
   - dock_report_{protein_name}_{substrate_names}.json

   The JSON report contains:

   - "output_type"
     A string identifying the report type:
     "enzywizard_dock"

   - "docked_result"
     A dictionary describing the best docking result.

     It includes:
     - "complex_name"
       Name of the docked protein-substrate complex.

     - "docking_score"
       Docking energy score from AutoDock Vina.

     - "substrate_names"
       Input substrate names.

     - "docking_box_center"
       Center coordinates of the docking box.

     - "docking_box_size"
       Size of the docking box.

     - "docked_substrates"
       A list describing each docked substrate.

       Each entry contains:
       - "substrate_name"
       - "conformation_name"
       - "docked_center_coord"

2. Files:
   - Docked SDF file for each substrate:
     docked_{substrate_name}.sdf

   - Docked protein-substrate complex CIF file:
     docked_{protein_name}_{substrate_names}.cif


# Process:

This command processes the input cleaned protein structure as follows:

1. Load the input structure
   - Read the cleaned CIF or PDB file using Biopython (Bio.PDB).
   - Resolve the protein name from the input filename.

2. Validate input conditions
   - Check that the input file exists.
   - Validate that the structure satisfies the cleaned-structure requirement.

3. Detect pocket regions
   - Use PyVOL to detect pocket regions from the protein structure.

4. Compute global docking box
   - Calculate a bounding box covering the entire protein structure.

5. Parse substrate inputs
   - Split substrate_names by ',' to obtain substrate list.

6. Search substrate files
   - Locate matched SDF files for each substrate in substrate_dir.

7. Enumerate substrate conformations
   - Treat multiple SDF files of the same substrate as alternative conformations.
   - Generate all combinations of substrate conformations.

8. Prepare docking inputs
   - Convert protein structure to receptor PDBQT format.
   - Convert each substrate SDF to ligand PDBQT format.

9. Build docking boxes
   - Use pocket-based boxes.
   - Add one global structure box.

10. Perform docking
   - Iterate over substrate combinations and docking boxes.
   - Perform AutoDock Vina docking for each case.

11. Parse docking results
   - Extract docking poses and energies.
   - Map docked coordinates back to original ligand atoms.

12. Select best result
   - Choose the docking result with lowest energy.
   - Optionally stop early if early_stop=True.
   
13. Save docking outputs
   - Write docked substrate SDF files.
   - Generate protein-substrate complex CIF file.

14. Generate report
   - Save structured JSON report summarizing docking results.


# dependencies:

- AutoDock Vina
- Meeko
- RDKit
- Biopython
- PyVOL
- MSMS


# references:

- Eberhardt et al., AutoDock Vina 1.2.0
  https://doi.org/10.1021/acs.jcim.1c00203

- Trott & Olson, AutoDock Vina
  https://doi.org/10.1002/jcc.21334

- AutoDock Vina:
  https://vina.scripps.edu/

- Meeko:
  https://github.com/forlilab/Meeko

- RDKit:
  https://www.rdkit.org/

- Biopython:
  https://biopython.org/

- PyVOL:
  https://github.com/schlessinger-lab/pyvol

- MSMS:
  https://doi.org/10.1002/jcc.540150805
