# PlasMap
Beta software for plasmid variation detection using high quality references.

This software has two functions:
1) Provide a similar functionality to BRIG, but with an improved visualization of genes on the reference plasmid.

2) Provide an automated analysis to find new genes that may be plasmid associated. This program maps annotated regions of contigs with significant plasmid homology. This can be especially helpful for draft genome assemblies created with Illumina technologies where multiple such contigs may exist.

This project is implemented as relatively simple Python code that relies on BLAST for gene associations.

Basic use:

This program uses annotated genbank files including sequence information as input. Displayed features use the "CDS" tag.

A closed reference plasmid is recommended for the output of this program to make sense.

Use of Miniforge 3 or Anaconda 3 as a package manager is strongly recommended. We tested the program with the following versions, but it will likely work on earlier versions of most packages:
python 3 (3.13.7) environment must have the following packages installed:
pandas (2.3.2), BioPython (1.85), reportlab (4.4.4), matplotlib (3.10.6), Ipython (9.5.0), xlrd (2.0.2), openpyxl (3.1.5), xlsxwriter (3.2.6)

Basic tools included in most Python installations are also required: os, sys, and itertools.

BLAST (2.17.0+) must also be installed to the command line. Any current version should work.

Running the program:

1) Place GenBank files for all plasmids you wish to compare plus the reference genbank

2) Run from the directory containing the python files: python3 PlasMap_Input.py Path_to_folder/with/all_genbank_files_for_comparison Reference.gb

This creates and populates output directory-> Path_to_folder/with/all_genbank_files_for_comparison_PlasMap

3) Fill in PlasMap_input.xlsx - use color names found in PlasMap_color_Ref.xlsx
For features that you do not wish to show annotations, set the "Skip" column to "Yes". You can also skip plasmids whose rings you do not wish to show.

4) Place any discriptions you wish to carry over automatically to new plasmid analyses into PlasMap_Defaults.xlsx. The default file has a small library for a Klebsiella analysis entered as an example.

5) Run from the directory containing the python files: python3 PlasMap_Image.py Path_to_folder/with/all_genbank_files_for_comparison Reference.gb Path_to_folder/with/all_genbank_files_for_comparison_PlasMap {font size, integer}

If you renamed the output folder, you should change the second location here. Replace "{font size, integer}" with a number

You can re-run the above code after changing and saving PlasMap_Input.xlsx. It will overwrite this analysis and replace the output files. This saves on clutter as you develop your analysis.


RUNNING THE DEMO

Either download this repository from the website or use git using this line from a location:

git clone https://www.github.com/worleyjn/PlasMap.git

Ensure you have all dependencies installed. Note: The conda channel "bioconda" is not available on windows, and BLAST will need to be installed separately from your package manager. Additionally, python3 may be called with just "python" in some systems.

Navigate to the PlasMap home directory containing PlasMap_Input.py and PlasMap_Image.py

Run: python3 PlasMap_Input.py Demo/NDM1_FIB_CP014757 NZ_CP014757.gb

There should be a new folder Demo/NDM1_FIB_CP014757_PlasMap

Replace the Excel file in the new folder "PlasMap_input.xlsx" in the new folder with the example in Expected_Results

Run: python3 PlasMap_Image.py Demo/NDM1_FIB_CP014757 NZ_CP014757.gb Demo/NDM1_FIB_CP014757_PlasMap 12

The output directory now has the diagram "Output.pdf", the legend "Legend.svg", the listing of sequence not in the reference "NewGenes.txt", and folder "Contig_diagrams" with a visual aid for results in "NewGenes.txt" showing genetic context. The output, legend, and new genes files should match.
