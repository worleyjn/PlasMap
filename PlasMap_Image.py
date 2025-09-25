import os
import sys
import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from reportlab.lib import colors
from reportlab.lib.units import cm
#this requires reportlab: https://docs.reportlab.com/install/open_source_installation/
from Bio.Graphics import GenomeDiagram
from IPython.core.display import Image
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from itertools import chain

#require the specific outdir here
script, indir, reffile, outdir, contig_label_fontsize = sys.argv

try:
    contig_label_fontsize = int(contig_label_fontsize)
except:
    print('CONTIG FONTSIZE IS NOT AN INTEGER!')
    sys.exit()

ref_record = SeqIO.read(os.path.join(indir, reffile), 'gb')
ref_len = len(ref_record.seq)
print(f'{ref_len} bases in reference plasmid.')

label_file = os.path.join(outdir, 'PlasMap_input.xlsx')
df = pd.read_excel(label_file, sheet_name = 'Features', index_col = 'Locus_tag')
dp = pd.read_excel(label_file, sheet_name = 'Plasmids', index_col = 'File')
dl = pd.read_excel(label_file, sheet_name = 'Legend', index_col = 'Legend_entry')

#hit length cutoff - ignore if the length of the contig is too small
#we're unlikely to get new genetic information, or consistent information
#smaller than this
len_co = 500

#proportional length cutoff
#how much of the contig mappping must not be already covered
#by better scoring hits
#some things may have been covered entirely by better scoring things
#proportion of hit that must be new
prop_co = 0.99

#start the diagram, can add new tracks piecewise
gd_diagram = GenomeDiagram.Diagram(ref_record.id, )
#first track - the genome, decorated
#smalltick interval set impossibly high
#start=1 needed to make large bar at sequence start
gd_ref_track = gd_diagram.new_track(1, name = reffile.replace('.gb', ''),
                                    height = 10,
                                    scale_largetick_interval = 10000.0,
                                    scale_smalltick_interval = 1000000000.0,
                                   start = 1, scale_fontsize = 15,
                                   scale_largeticks = 0.1,
                                   scale_smallticks = 0)
ref_features = gd_ref_track.new_set()
for feature in ref_record.features:
    #looking only at CDS here
    if feature.type != "CDS":
        #Exclude this feature
        continue
    locus_tag = feature.qualifiers['locus_tag'][0]
    if str(df['Skip'][locus_tag]).lower() != 'no':
        name = ''
        color = 'lightgrey'
        label = False
    else:
        name = str(df['Display_name'][locus_tag])
        color = str(df['Color'][locus_tag])
        label = True
    ref_features.add_feature(feature, sigil = 'ARROW', name = name,
                               color = color, label = label,
                               label_size = 10, label_strand = -1,
                            label_position = 'middle')
qdir = os.path.join(outdir, 'query_files')
#add a new track for each blast file
#print out genes that do not align and the contig they are on

g = open(os.path.join(outdir, 'NewGenes.txt'), 'w')

inset_map_pairs = []

for num, gbfile in enumerate(dp.index):
    if str(dp['Skip'][gbfile]).lower() != 'no':
        continue
    color = dp['Color'][gbfile]
    blastfile = gbfile.replace('.gb', '.blast')
    qdict = SeqIO.to_dict(SeqIO.parse(os.path.join(qdir, blastfile.replace('.blast', '.fna')), 'fasta'))
    hits = [i.split('\t') for i in open(os.path.join(qdir, blastfile),'r').read().splitlines()]
    #sort by bitscore
    s_hits = sorted(hits, key = lambda x: float(x[-1]), reverse = True)
    #where on the query sequence have we mapped
    #qmapped[contig] = [list of covered bases]
    qmapped = {}
    smapped = set()
    for hit in s_hits:
        qid, sid, pid, length, mism, gaps, qstart, qend, sstart, send, evalue, bitscore = hit
        qlen = len(qdict[qid].seq)
        #apply length cutoff
        if int(length) < len_co:
            continue
        #print(qid[-4:], length, qlen, pid, sstart, send)
        up, down = sorted([int(i) for i in [sstart, send]])
        #print(up, down)
        bases_hit = set([i for i in range(up, down+1)])
        #apply the bases hit cutoff
        if len(smapped.intersection(bases_hit))/int(length) >= prop_co:
            continue
        smapped = set(list(smapped) + list(bases_hit))
        if not qid in qmapped:
            qmapped[qid] = []
        qup, qdown = sorted([int(i) for i in [qstart, qend]])
        qbases_hit = set([i for i in range(qup, qdown+1)])
        qmapped[qid] = set(list(qmapped[qid]) + list(qbases_hit))
    
    #set up the new track here
    aln_track = gd_diagram.new_track(num+2, name=blastfile.replace('.blast', ''),
                                    scale_ticks=0, scale_smallticks=0)
    aln_features = aln_track.new_set()
    
    #find covered ranges
    s_smapped = list(sorted(smapped))
    beginnings = [s_smapped[0],]+[s_smapped[i] for i in range(1, len(s_smapped)) if not s_smapped[i] - s_smapped[i-1] == 1]
    ends = [s_smapped[i] for i in range(len(s_smapped)-1) if not s_smapped[i+1] - s_smapped[i] == 1] + [s_smapped[-1],]
    for tup in zip(beginnings, ends):
        b, e = tup
        new_feature = SeqFeature(FeatureLocation(b, e))
        aln_features.add_feature(new_feature, sigil = 'BOX', color = color)

    
    #print where to look on contigs with less mapping for interesting genes
    print(gbfile)
    g.write(f'{gbfile}\n')
    gbdict = SeqIO.to_dict(SeqIO.parse(os.path.join(indir, gbfile), 'gb'))
    for qid in qmapped:
        #set of hit positions
        sqid = qmapped[qid]
        print(qid, len(qdict[qid].seq), len(qmapped[qid]))
        g.write(f'{qid}; Length = {len(qdict[qid].seq)}; Mapped bases = {len(qmapped[qid])}\n')
        mapbool = False
        for feature in gbdict[qid].features:
            if not feature.type == 'CDS':
                continue
            else:
                pass
            #print(feature)
            loc1 = int(feature.location.start)
            loc2 = int(feature.location.end)
            uploc, downloc = sorted([loc1, loc2])
            slocs = set(range(uploc, downloc+1))
            if len(sqid.intersection(slocs)) < 0.6*len(slocs):
                print(f'\t{loc1}\t{loc2}\t\t{feature.qualifiers["product"]}')
                g.write(f'\tUp={loc1}\tDown={loc2}\t\t{feature.qualifiers["product"]}\n')
                mapbool = True
             
        #only pass if there are things to map
        if mapbool == False:
            continue
        else:
            pass
        #collecting pairs to map for later
        inset_map_pairs.append((gbfile, qid))            
                
                
    
    print('\n\n')
    g.write('\n\n')
g.close()
gd_diagram.draw(format='circular', circular=True, pagesize=(50*cm,50*cm),
                start=0, end=len(ref_record), circle_core = 0.7)
gd_diagram.write(os.path.join(outdir, 'Output.pdf'), 'pdf')




#make a legend that can be pasted in
color_ref = pd.read_excel('PlasMap_color_Ref.xlsx', index_col = 'name')

fig, ax = plt.subplots()

handles = []
labels = []
for entry in dl.index:
    labels.append(entry)
    colorname = dl['Color'][entry]
    patch = mpatches.Patch(color = f'#{color_ref["hex"][colorname][2:]}', label = entry)
    handles.append(patch)
ax.legend(handles=handles)
legend = plt.legend(handles, labels, framealpha=1, frameon=False, bbox_to_anchor=(2, 2))

def export_legend(legend, filename = 'Legend.svg'):
    fig  = legend.figure
    fig.canvas.draw()
    bbox  = legend.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(os.path.join(outdir, filename), dpi="figure", bbox_inches=bbox)
export_legend(legend)


###Plot each of the contigs with significant alignments
###Plot it good

mp_outdir = os.path.join(outdir, 'Contig_diagrams')

if not os.path.isdir(mp_outdir):
    os.mkdir(mp_outdir)
for pair in list(set(inset_map_pairs)):
    gbfile, contig = pair
    genome = os.path.splitext(gbfile)[0]
    
    blastfile = gbfile.replace('.gb', '.blast')
    qdict = SeqIO.to_dict(SeqIO.parse(os.path.join(qdir, blastfile.replace('.blast', '.fna')), 'fasta'))
    hits = [i.split('\t') for i in open(os.path.join(qdir, blastfile),'r').read().splitlines() if i.startswith(contig)]
    #sort by bitscore
    s_hits = sorted(hits, key = lambda x: float(x[-1]), reverse = True)
    #query/subject start and ends
    #do not sort these, parallel structures
    qses = []
    sses = []
    smapped = set()
        
    for hit in s_hits:
        #keep the code from above to get only quality hits
        qid, sid, pid, length, mism, gaps, qstart, qend, sstart, send, evalue, bitscore = hit
        qlen = len(qdict[qid].seq)
        #apply length cutoff
        if int(length) < len_co:
            continue
        #print(qid[-4:], length, qlen, pid, sstart, send)
        up, down = sorted([int(i) for i in [sstart, send]])
        #print(up, down)
        bases_hit = set([i for i in range(up, down+1)])
        #apply the bases hit cutoff
        if len(smapped.intersection(bases_hit))/int(length) >= prop_co:
            continue
        smapped = set(list(smapped) + list(bases_hit))
        qses.append((int(qstart), int(qend)))
        sses.append((int(sstart), int(send)))
        print(hit)
    qsses = list(zip(qses, sses))
    #put all values in a list
    flat_qses = list(chain(*qses))
    qmin = min(flat_qses)
    qmax = max(flat_qses)
    flat_sses = list(chain(*sses))
    smin = min(flat_sses)
    smax = max(flat_sses)
    #print(qmin, qmax)
    #print('q_hit bounds, s_hit bounds', qsses)
    #print(smin, smax)
    #print(sses)
    #determine if the comparitor needs to be reverse complemented
    #weight by just total sequence
    #direction score - use the subject sequence
    #query sequence always positive - not really using this right now
    #since implementation is a bit of a mess and it's not 
    #completely necessary, it's a nice-to-have for later
    ds = sum([i[1]-i[0] for i in sses])
    
    #Placeholder for reverse-complementation function
    #Keeping it forward - to - forward at this point
    #Keeps it easier for finding contig positions
    #if ds > 0:
        #rc = False
    #else:
        #continue
        #rc = True
    
    #Initialize the diagram
    gd_diagram = GenomeDiagram.Diagram(ref_record.id, fragments = 1, track_size = 0.2)
    
    #Find key positions in the diagram
    #These are overall coordinates   
    
    #different track ticks for different sizes
    #always 1/5th length
    main_track_length = smax - smin
    if main_track_length <= 5000:
        rslti = 500
        rssti = 100
    elif main_track_length <= 10000:
        rslti = 1000
        rssti = 200
    elif main_track_length <= 50000:
        rslti = 5000
        rssti = 1000
    else:
        rslti = 10000
        rssti = 2000
        
        
        
    gd_ref_track = gd_diagram.new_track(2, name=reffile.replace('.gb', ''),
                                    height = 1,
                                    scale_largetick_interval = rslti,
                                    scale_smalltick_interval = rssti,
                                   start = smin, end = smax, scale_fontsize = 15,
                                   scale_largeticks = 0.5, scale_smallticks = 0.2,
                                       scale_fontangle = -45, )
    #add the reference features
    #these can be direct, it's fine
    ref_features = gd_ref_track.new_set()
    for feature in ref_record.features:
        if feature.type != "CDS":
            #Exclude this feature
            continue
        
        locus_tag = feature.qualifiers['locus_tag'][0]
        if str(df['Skip'][locus_tag]).lower() != 'no':
            name = ''
            color = 'lightgrey'
            label = False
        else:
            name = str(df['Display_name'][locus_tag])
            color = str(df['Color'][locus_tag])
            label = True
        ref_features.add_feature(feature, sigil = 'ARROW', name = name,
                                   color = color, label = label,
                                   label_size = contig_label_fontsize, label_strand = 1,
                                label_position = 'middle',
                                label_angle = 20)
        
    #need to shift track over to start at the same place
    #smin will be the starting position
    
    #find the distance markers we should 
    #hit track should start at smin value, then extend for the span of the hits
    
    q_records = SeqIO.to_dict(SeqIO.parse(os.path.join(indir, gbfile), 'gb'))
    contig_length = len(q_records[contig].seq)
        
    gd_hit_track = gd_diagram.new_track(1, name=reffile.replace('.gb', ''),
                                    height = 1,
                                    scale_largetick_interval = rslti,
                                    scale_smalltick_interval = rssti,
                                    start = smin, end = smin+contig_length, scale_largetick_labels = 0,
                                    scale_smalltick_labels = 0,
                                    scale_fontsize = 15, scale_fontangle = -45,
                                        scale_largeticks = 0.5, scale_smallticks = 0.2)

    #need to make new features with the correct modified position
    #recall -> rc = reverse complement
    #begin at position smin
    #new_features = []
    hit_features = gd_hit_track.new_set()
    for feature in q_records[contig].features:
        #skip non-CDS features
        if feature.type != "CDS":
            continue
        #print(feature.location)
        
        fstart = feature.location.start+smin-1
        fend = feature.location.end+smin-1
        strand = feature.location.strand
        locus_tag = feature.qualifiers['locus_tag'][0]
        feature.location = FeatureLocation(start = fstart, end = fend, strand = strand)
        
        product = feature.qualifiers['product'][0]
        found_desc = False
        for locus_tag in df.index:
            if df['Product'][locus_tag] == product:
                name = df['Display_name'][locus_tag]
                color = df['Color'][locus_tag]
                label = True
                found_desc = True
                break
        
        if found_desc == True:
            pass
        else:
            name = product
            color = 'lightgrey'
            label = True
        if type(color) == float:
            color = 'lightgrey'
        if type(name) == float:
            name = ''
        #print(color)
        hit_features.add_feature(feature, sigil = 'ARROW', name = name,
                         color = color, label = label,
                         label_size = contig_label_fontsize, label_strand = -1,
                         label_position = 'middle',
                                label_angle = 160)
        
        #reverse complementation flags
        #if I/we decide/demand to implement
        #if rc == False:
        #   fstart = feature.location.start+smin-1
        #    fend = feature.location.end+smin-1
        #    strand = feature.location.strand
        #    locus_tag = feature.qualifiers['locus_tag'][0]
        #    feature.location = FeatureLocation(start = fstart, end = fend, strand = strand)
        #    if not locus_tag in df.index:
        #        name = ''
        #        color = 'lightgrey'
        #        label = False            
        #    elif str(df['Skip'][locus_tag]).lower() != 'no':
        #        name = ''
        #        color = 'lightgrey'
        #        label = False
        #    else:
        #        name = str(df['Display_name'][locus_tag])
        #        color = str(df['Color'][locus_tag])
        #        label = True
        #    hit_features.add_feature(feature, sigil = 'ARROW', name = name,
        #                           color = color, label = label,
        #                           label_size = 20, label_strand = -1,
        #                        label_position = 'middle')
            #print(fstart, fend, strand)
        #else:
            #continue
            
    #mark the homologous zones
    #(track object, start, end)
    #class Bio.Graphics.GenomeDiagram.CrossLink(featureA, featureB, color=Color(0.564706, 0.933333, 0.564706, 1), border=None, flip=False)
    for aln in qsses:
        #query range, subject range
        qr, sr = aln
        print(qr, sr)
        q1, q2 = qr
        qup = smin+q1-1
        qdown = smin+q2-1
        qtup = (gd_hit_track, qup, qdown)
        s1, s2 = sr
        #bc - box color
        #ec - edge color
        if s1 > s2:
            fliptup = True
            clbc = colors.lightpink
            clec = colors.sienna
        else:
            fliptup = False
            clbc = colors.skyblue
            clec = colors.steelblue
        sup, sdown = sorted(sr)
        stup = (gd_ref_track, sup, sdown)
        #print(qtup, stup)
        cl = GenomeDiagram.CrossLink(qtup, stup, flip = fliptup, color = clbc, border = clec)
        gd_diagram.cross_track_links.append(cl)
        
    
    #adjust smin and smax
    gd_diagram.draw(format='linear',
                    circular=False,
                    pagesize=(20*cm,100*cm),
                    start = smin,
                    end = max([smax, smin+contig_length]))
    
    outfile_base = '_'.join([genome, qid])
    gd_diagram.write(os.path.join(mp_outdir, f'{outfile_base}.svg'), 'svg')
    
    print('\n\n')