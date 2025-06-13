# PlasMap
Beta software for plasmid variation detection using high quality references.

This software has two functions.
1) Provide a similar functionality to BRIG, but with an improved visualization of genes on the reference plasmid.
2) Provide an automated analysis to find new genes that may be plasmid associated, but on other contigs. Especially helpful for draft genome assemblies created with Illumina, and some genomic rearrangement.

This project is implemented as relatively simple Python code that relies on BLAST for gene associations. It isn't super sophisticated, and that is by design. Users of the program should be able to easily modify it to fit their systems.

Basic use:

#use an appropriate virutal environment with prerequistes installed.

python3 PlasMap_Input.py ../Path_to_folder/with/all_genbank_files_for_comparison Reference.gb

creates output directory-> ../Path_to_folder/with/all_genbank_files_for_comparison_PlasMap
Fill in PlasMap_input.xlsx - use color names found in PlasMap_color_Ref.xlsx
Place any discriptions you wish to carry over automatically to new plasmid analyses into PlasMap_Defaults.xlsx

#default directory created when running example code above shown for PlasMap directory located.
python3 PlasMap_Image.py ../Path_to_folder/with/all_genbank_files_for_comparison Reference.gb ../Path_to_folder/with/all_genbank_files_for_comparison_PlasMap {font size, integer}

You can re-run the above code after changing and saving PlasMap_Input.xlsx
It will overwrite this analysis.
