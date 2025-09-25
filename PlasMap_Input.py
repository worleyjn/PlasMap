import os
import sys
import pandas as pd
from pandas import DataFrame
from Bio import SeqIO
from subprocess import run

script, indir, reffile = sys.argv

class feat:
    def __init__(self, ptag, disname, color, skip):
        self.ptag = ptag
        self.disname = disname
        self.color = color
        self.skip = skip

#we can always go home, probably
homedir = os.getcwd()

if not os.path.isdir(indir):
    print(f'INDIR {indir} DOES NOT EXIST')
    sys.exit()
    
if not os.path.isfile(os.path.join(indir, reffile)):
    print(f'{reffile} does not exist in {indir}')

#input gb files
ingbs = [i for i in os.listdir(indir) if (i.endswith('.gb')) and not i.startswith('.') and not i == reffile]

if len(ingbs) == 0:
    print('No valid comparison files.')
    sys.exit()
    
if not os.path.isfile('PlasMap_Defaults.xlsx'):
    print('No PlasMap_Defaults.xlsx default file.')
    
try:
    dxl = pd.ExcelFile('PlasMap_Defaults.xlsx')
    dfeats = pd.read_excel(dxl, sheet_name = 'Features')
    dleg = pd.read_excel(dxl, sheet_name = 'Legend', index_col = 0)
    feat_dict = {}
    for ind in dfeats.index:
        ltag, ptag, disname, color, skip = dfeats.loc[ind]
        if not ptag in feat_dict:
            feat_dict[ptag] = feat(ptag, disname, color, skip)
    
    leg_dict = {}
    for entry in dleg.index:
        leg_dict[entry] = dleg['Color'][entry]
except:
    print('Did not load in default colors from file correctly.')
    sys.exit()

#load in the reference gb record
ref_record = SeqIO.read(os.path.join(indir, reffile), 'gb')

#make an output directory name to start
outdir = '_'.join([indir.rstrip('/'), 'PlasMap'])

#make new dir or clean out contents for new analysis
if not os.path.isdir(outdir):
    os.mkdir(outdir)
    print(f'Output directory: {outdir}')
else:
    #make default behavior to make a new directory
    #with a different name if one already exists
    print(f'Output directory \'{outdir}\' exists, adding suffix')
    num = 1
    while True:
        newoutdir = f'{outdir}{num}'
        if os.path.isdir(newoutdir):
            print(num)
            num+=1
            continue
        else:
            outdir = str(newoutdir)
            os.mkdir(outdir)
            print(f'Output directory: {outdir}')
            break

#make the reference coloring sheet
features = [i for i in ref_record.features if i.type == 'CDS']
locus_tags = []
products = []
genes = []
###Need to add in the defaults!

for feature in features:
    try:
        locus_tag = feature.qualifiers['locus_tag']
        product = feature.qualifiers['product']
        if 'gene' in feature.qualifiers:
            gene = feature.qualifiers['gene']
        else:
            gene = ['',]
    except:
        print(feature, '\nFEATURE ABOVE DID NOT HAVE LOCUS TAG AND OR PRODUCT')
        continue
    locus_tags.append(locus_tag[0])
    products.append(product[0])
    genes.append(gene[0])
    
df = DataFrame(columns = ['Locus_tag', 'Product', 'Display_name', 'Color', 'Skip'])
df['Locus_tag'] = locus_tags
df['Product'] = products
df['Display_name'] = genes
df['Skip'] = ['No' for i in products]
#if the product is in the features dictionary, add in the relevant data
for ind in df.index:
    product = df['Product'][ind]
    if product in feat_dict:
        feat = feat_dict[product]
        disname = feat.disname
        color = feat.color
        skip = feat.skip
        df.loc[ind, ['Display_name', 'Color', 'Skip']] = [disname, color, skip]
    else:
        continue

df.set_index = 'Locus_tag'

#make the comparison genome coloring sheet
dp = DataFrame(columns = ['File', 'Color', 'Skip'])
dp['File'] = ingbs
dp['Skip'] = ['Yes' for i in ingbs]

#make the legend sheet
dl = DataFrame(columns = ['Legend_entry', 'Color'])
entries = []
colors = []
for entry in leg_dict:
    entries.append(entry)
    colors.append(leg_dict[entry])
dl['Legend_entry'] = entries
dl['Color'] = colors

#write out the sheets
sheetnames = ['Features', 'Plasmids', 'Legend']
dfs = [df, dp, dl]
g = pd.ExcelWriter(os.path.join(outdir, 'PlasMap_input.xlsx'), engine='xlsxwriter')
for sheet, d in zip(sheetnames, dfs):
    d.to_excel(g, sheet_name = sheet, index = False)
g.close()

            
#First, set up the blast database from the reference plasmid
blastdbdir = os.path.join(outdir, 'blastdb')
reffasta = reffile.replace('.gb', '.fna')
os.mkdir(blastdbdir)
#navigate to that directory
os.chdir(blastdbdir)
#make a fastafile in the indir
print(os.path.join(homedir, indir, reffile))
try:
    with open(reffasta, 'w') as f:
        for contig in SeqIO.parse(os.path.join(homedir, indir, reffile), 'gb'):
            SeqIO.write(contig, f, 'fasta')
    print('Made reffile')
    #will not work in testing environment
    run(['makeblastdb','-dbtype', 'nucl', '-in', reffasta])
    print('BLASTn database made')
except:
    print(f'Could not parse file {reffile}.')
    sys.exit()

#back in homedir, make query file dir and navigate there
os.chdir(homedir)
qdir = os.path.join(outdir, 'query_files')
os.mkdir(qdir)
os.chdir(qdir)

#make a fasta copy
for qfile in ingbs:
    qfasta = qfile.replace('.gb', '.fna')
    with open(qfasta, 'w') as f:
        for contig in SeqIO.parse(os.path.join(homedir, indir, qfile), 'gb'):
            SeqIO.write(contig, f, 'fasta')
    #do the blast searches
    run(['blastn', '-db', f'../blastdb/{reffasta}', '-query', qfasta, '-outfmt', '6', '-out', qfasta.replace('.fna', '.blast'), '-evalue', '0.001'])

