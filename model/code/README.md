# How to use Mimics1.5a

> Mimics1.5a is a code for simulating the microwave backscatter from a forest.

## Overview


MIMICS1.5a was written at the University of Michigan, College of Engineering, Electrical Engineering and Computer Science Department, The Radiation Laboratory by Kyle McDonald, with some extensions by Leland Pierce.

Two new ground layer models have been added since the original was written:

1. Purely specular (no backscatter at all)
2. "UMich Empirical" model, based on the article:
Yisok Oh, K. Sarabandi, Fawwaz Ulaby: An Empirical Model and an Inversion Technique for Radar Scattering from Bare Soil Surfaces, IEEE Transactions on Geoscience and Remote Sensing, vol 30., No. 2, 1992. pages 370-381.

The theoretical documentation for MIMICS is given in the code/mimics.pdf file.
This is a pdf document that you can either view on-line or download and print.
It contains the basic theory behind MIMICS along with some example runs.


## Input data

If you look inside the Data directory you will see many files.
Each of these files contains input parameters that are needed in order to run MIMICS.

OK, so there are way too many input files. How does one go about specifying the parameters for a run?
Below is the order in which to edit the files, as well as a quick description of what paraemters are in each:

1. configuration.input -- general overall params for the simulation,
2. sensor.input        -- operating frequencies and incidence angles,
3. environment.input   -- temperatures,
4. ground.input        -- ground params and model,
5. dielectric_ground_table.input -- optionally specify the
                                 ground dielectric constant.
6. dielectric_snow_table.input -- optionally specify the
                                  snow layer dielectric constant.
 
Then optionally edit the remaining vegetation files:

7. trunk_and_gross_canopy.input
8. dielectric_trunk_table.input
9. histogram_trunk_diam.input
10. trunk_hght_f.f  -- relating trunk height (in meters) to diameter (cm),
       This file may be created by using the program make_trnkhgt
       if a simple polynomial of 4th order or less will suffice.
       (you'll need to download make_trnkhgt.f from the src directory and compile and run on your own computer)

11. branch_primary.input
12. dielectric_branch1_table.input
13. branch_secondary.input
14. dielectric_branch2_table.input
       (there are 4 more branch classes, with these same 2 files for each,
        the filenames are obvious)

15. leaf.input
16. dielectric_leaf_table.input
17. needle.input
18. dielectric_needle_table.input
 
19. parameter_nesting.input


At this point there is one file left to fill in:
       parameter_value_table.input
This file will contain the values of parameters that are not to be stepped using a min/max/delta type rule. Hence this file provides for entering several values for a parameter that it will then be stepped through, providing the ability to use non-uniform step sizes. To assist in creating this file, there is a program in the src directory called "make_table.f" that will create this file with empty spaces to fill in with the parameter values. (download, compile and run on your own computer.)

This step must be LAST because it reads the other input files in order to create the right number of parameters in the right order. This is only necessary if any of the parameters are being varied by table entries as specified by the other input files.

## Results

After running MIMICS there will be many files creted in the Results directory.
You should read the theoretical overview (in code/mimics.pdf) in order to understand the terminology used in this section.

Here's an overview of what each contains:
1. forest_kappa_mats.out:  Kappa matrices of individual layers
2. forest_phase_crown.out: Phase differences for each crown-layer mechanism
3. forest_phase_mats.out:  Phase matrices for each full-canopy mechanism
4. forest_phase_terms.out: Phase differences for each full-canopy mechanism
5. forest_phase_trunk.out: Phase differences for each trunk-layer mechanism
6. forest_sigma_cross.out: Cross-polarized Sigma-0's for each full-canopy mechanism
7. forest_sigma_like.out:  Co-polarized Sigma-0's for each full-canopy mechanism
8. forest_trans_mats.out:  Transmissivity matrices for individual layers
9. forest_trans.out:       Transmissivity values for each individual layer
10. mimics2.data:          Experimental
11. polarimetric.out:      Transformation matrices for each full-canopy mechanism

If all you need is the backscattered powers (4 polarizations) these are in forest_sigma_like.out and forest_sigma_cross.out.

Note that the data is printed out in such a way that you need not have access to the input files to understand the output data.
All the data that is not varied is printed at the top, with the varied parameters and the outputs in a tabular form beneath.
This is meant more for human consumption than for another computer code to use.


