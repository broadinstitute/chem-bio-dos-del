# chem-bio-dos-del
Initiated 2021Q4 for code related to the Broad Chemical Biology DNA-encoded library (DEL) analysis and visualization pipeline.

## Content
### /del-git-dev/bruce-hua/tools_utilities/
- lib_enumeration_example_20220515_fixed.py
 
     *Bruce Hua "...put together a python script that consolidates it all in one place, the creation of new library metadata. Should just have to specify the dir_name and lib_id variable." (Also, need to customized the reaction SMARTS strings to match the structures in the new library. For example, **rxn1 = rdChemReactions.ReactionFromSmarts('[CX3:1](=[OX1])NC>>[CX3:1](=[OX1])NC[3H]')** )*
     
- /generate_stats/generate_stats_table.R

     *Bruce Hua put together an R script that generates count generation statistics of samples':*
     
        - total reads,
        
        - total unique reads,
        
        - valid reads,
        
        - valid barcodes,

        - total counts

- /generate_stats/analyses

     *location of generate_stats_table.R outputs (with some example statistics)*
     
- populate_lib_enum_structure_008.py

     *python code to generate library 8 related data, to be appended to lib_enum_struct.csv (which is used during count generation)*

### /del-git-dev/bruce-hua/shiny
- config_file.json
- app.R
 
    *The two files for deploying the Shiny interface of the DEL analysis web application*
### /del-git-dev/bruce-hua/pipe
### /del-git-dev/bruce-hua/fun
- function.R 

    *contains functions that the /shiny/app.R file calls upon*

- metadata.R 

    *opens all metadata files for the /shiny/app.R file to read*

### /del-git-dev/bruce-hua/bash
### /del-git-dev/connor-coley
### /del-git-dev/sophia-lai

