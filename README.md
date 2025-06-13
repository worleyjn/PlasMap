# PlasMap
Beta software for plasmid variation detection using high quality references.

This software has two functions.
1) Provide a similar functionality to BRIG, but with an improved visualization of genes on the reference plasmid.
2) Provide an automated analysis to find new genes that may be plasmid associated, but on other contigs. Especially helpful for draft genome assemblies created with Illumina, and some genomic rearrangement.

This project is implemented as relatively simple Python code that relies on BLAST for gene associations. It isn't super sophisticated, and that is by design. Users of the program should be able to easily modify it to fit their systems.
