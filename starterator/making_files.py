#!/usr/bin/env python

import cPickle
import argparse
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter, A4
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
import PyPDF2
from Bio import SeqIO
import math
import StringIO


parser = argparse.ArgumentParser()
parser.add_argument('-n', '--pham_no', default = -1)
parser.add_argument('-p' , '--phage', default=None)
parser.add_argument('-a', '--one_or_all', default='All')
parser.add_argument('-f', '--pickle_file', default=None)
parser.add_argument('-o', '--output_dir')
parser.add_argument('-d', '--final_dir', default=None)
parser.add_argument('-m', '--make', default=None)

args = parser.parse_args()

output_dir = args.output_dir.strip('"') #
final_path = args.final_dir.strip('"')

def add_pham_no_title(pham_no, first_graph_path, i=''):
    # print i, type(i)
    # print first_graph_path
    packet = StringIO.StringIO()
    can = canvas.Canvas(packet, pagesize=letter)
    width, height = letter
    # print width, height
    can.drawString(280, 750, 'Pham ' + str(pham_no))
    can.save()

    packet.seek(0)
    new_pdf = PyPDF2.PdfFileReader(packet)
    existing_pdf = PyPDF2.PdfFileReader(file(first_graph_path), 'rb')
    output = PyPDF2.PdfFileWriter()

    page = existing_pdf.getPage(0)
    page.mergePage(new_pdf.getPage(0))
    output.addPage(page)

    outputStream = file("%sPham%sGraph%s.pdf" % (output_dir+ pham.phage+ args.one_or_all, pham.pham_no, i),'wb')
    # print outputStream
    output.write(outputStream)
    outputStream.close()

def combine_graphs(phage, pham_no, num_pages):
    merger = PyPDF2.PdfFileMerger()
    for i in xrange(0, num_pages+1):
        graph = open("%sPham%sGraph%d.pdf"  % (output_dir+ phage + args.one_or_all, pham_no, i), "rb")
        merger.append(fileobj=graph)
    merger.write(open("%sPham%sGraph.pdf"  % (output_dir+ phage + args.one_or_all, pham_no), "wb"))

def make_gene_track(gd_diagram, gene_group, num_on_diagram, total):
    """"""
    colors = ['purple', 'red', 'green', 'orange', 'yellow', 'brown'] 
    aligned_gene = gene_group[0]
    gd_gene_track = gd_diagram.new_track(total - num_on_diagram, name='Track %s'  % (num_on_diagram+1), 
                            label=True, greytrack=1)
    gd_seq_set = gd_gene_track.new_set()
    gd_feature_set = gd_gene_track.new_set()

    start_site = aligned_gene.annotations['start']
    start_site_feature = SeqFeature(FeatureLocation(start_site, start_site +1), 
                            strand=None)
    for feature in aligned_gene.features:
        if feature.type == 'seq':
            gd_seq_set.add_feature(feature, color='pink')
    for site in aligned_gene.annotations["possible_starts"]:
        site_color = pham.total_possible_starts.index(site) % len(colors) 
        possible_site = SeqFeature(FeatureLocation(site, site ), strand=None)
        gd_feature_set.add_feature(possible_site, color=colors[site_color], 
            name=str(pham.total_possible_starts.index(site)+1), label=True)
    end_gene_feature = SeqFeature(FeatureLocation(len(aligned_gene), 
                        len(aligned_gene)+1), strand=None)
    gd_feature_set.add_feature(start_site_feature, color="blue", label=True)
    gd_feature_set.add_feature(end_gene_feature, color='purple', label=True)

def graph_start_sites(pham, file_path):
    """
        graphs the alignment, creates a PDF file called {Phage Name}{One or All}Pham{Pham Number}.pdf
    """

    # genes = sorted(pham.genes_in_pham.values())
    genes = pham.groups
    # for group in genes:
    #     print group, group.id
    if pham.phage == None:
        pham.phage = ""
    seq_length = len(genes[0][0].seq)
    if len(genes) > 100:
        for i in xrange(0, int(math.ceil(len(genes)/50.0))):
            gd_diagram = GenomeDiagram.Diagram(pham.pham_no)
            graph_path = "%sPham%sGraph_%s.pdf" % (file_path+pham.phage+ args.one_or_all, pham.pham_no, i)
            graph_path_svg = "%sPham%sGraph_%s.svg" % (file_path+pham.phage+ args.one_or_all, pham.pham_no, i)

            for j in xrange(0, 50):
                if i*50 + j >= len(genes):
                    print i * 50, + j, len(genes)
                    gd_gene_track = gd_diagram.new_track(50-j)
                    gd_feature_set = gd_gene_track.new_set()
                    empty_feature = SeqFeature(FeatureLocation(0,1), 
                            strand=None)
                    gd_feature_set.add_feature(empty_feature, color="black", label=True)
                    print 'blank track added'
                else:
                    gene = genes[i*50 + j][0]
                    make_gene_track(gd_diagram, genes[i*50 + j], j, 50)
            print seq_length, i

            gd_diagram.draw(format="linear", orientation="portrait", pagesize=letter, 
                fragments=1, start=0, end=seq_length)
            gd_diagram.write(graph_path, "PDF")
            gd_diagram.write(graph_path_svg, "SVG")

            add_pham_no_title(pham.pham_no, graph_path, str(i))

        combine_graphs(pham.phage, pham.pham_no, i)
    else:            
        graph_path_svg = "%sPham%sGraph.svg" % (file_path+pham.phage+ args.one_or_all, pham.pham_no)
        graph_path = "%sPham%sGraph_.pdf" % (file_path+pham.phage+ args.one_or_all, pham.pham_no)
        gd_diagram = GenomeDiagram.Diagram(pham.pham_no)
        i = 0
        for gene_group in genes:
            print 'group', i
            make_gene_track(gd_diagram, gene_group, i, len(genes))
            i += 1
        gd_diagram.draw(format="linear", orientation="portrait", pagesize=letter, 
            fragments=1, start=0, end=len(gene_group[0].seq))
        gd_diagram.write(graph_path, "PDF")
        gd_diagram.write(graph_path_svg, "SVG")
        add_pham_no_title(pham.pham_no, graph_path)

        # gd_diagram.write("%s.svg" % (file_path+pham.pham_no), "SVG")
        # gd_diagram.write("%s.eps" % (file_path+pham.pham_no), "EPS")
        # gd_diagram.write("%s.png" % (file_path+pham.pham_no), "PNG")

def make_pham_text(pham, pham_no, output_dir, only_pham=False):
    """
        Creates a PDF Report for the specific pham.
        From Start sites statisitics
        phage specific
    """
    if not args.phage:
        name = "%s%sPham%sText.pdf" % (output_dir, args.one_or_all, pham_no)
    else:
        name = "%s%sPham%sText.pdf" % (output_dir, args.phage + args.one_or_all, pham_no)
    doc = SimpleDocTemplate(name, pagesize=letter)
    story = []
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="paragraph"))
    styles.add(ParagraphStyle(name='Center', alignment=TA_CENTER))
    text = '<font size=14> Pham %s Report </font>' % pham_no
    story.append(Paragraph(text, styles['Center']))
    story.append(Spacer(1, 12))
    groups = pham.groups
    tracks_info = []
    for index in range(len(groups)):
        text = "Track %s : " % (index+1)
        for gene in groups[index]:
            text += gene.id + ', '
        tracks_info.append("<font size=12> "+ u'\u2022'+" %s</font>" % text)
    for line in tracks_info:
        story.append(Paragraph(line, styles["Normal"]))
    story.append(Spacer(1, 12))
    if only_pham:
        output = pham.output_start_sites()
        for line in output:
            if line == '':
                story.append(Spacer(1, 12))
            text = '<font size=12> %s </font>' % line
            # if 'Genes' not in line or '':
            story.append(Paragraph(text, styles['Normal']))
            # else:
            #     story.append(Paragraph(text, styles['Normal']))
        # story.append()
        # story.append(Paragraph("Suggested Starts:", styles["Normal"]))
        suggested_start = pham.output_suggested_starts(all_genes=True)
        story.append(Spacer(1, 12))

        for line in suggested_start:
            text = '<font size=12>%s</font>' % line
            story.append(Paragraph(text, styles["Normal"]))
    else:
        story.append(Paragraph("<font size=14>Suggested Starts: </font>", styles["Normal"]))
        suggested_start = pham.output_suggested_starts()
        story.append(Spacer(1, 12))

        for line in suggested_start:
            text = '<font size=12>%s</font>' % line
            story.append(Paragraph(text, styles["Normal"]))
            # story.append()
        list_candidate_starts = pham.candidate_starts_phage_genes()
        for gene in list_candidate_starts:
            gene_record = pham.genes_in_pham[gene]
            start, stop = gene_record[0].annotations['start'], gene_record[0].annotations['stop']
            candidate_starts = list_candidate_starts[gene]
            # for start in candidate_starts:
            #     text_starts += start
            story.append(Paragraph("<font size = 12> Gene: %s \n Start: %s, Stop: %s </font>" % (gene_record[1].id, start+1, stop), styles["Normal"]))
            story.append(Paragraph("<font size = 12> Candidate Starts for %s: </font>" % (gene_record[1].id), styles["Normal"]))
            story.append(Paragraph("<font size = 12>"+ str(candidate_starts) + "</font>" , styles["Normal"]))
    # suggested_start = pham.output_suggested_starts()
    doc.build(story)

def one_pham_report(pham_no, phage):
    """
        Merges the graphical and text output of one pham into one PDF file
        Pham{Number}Report.pdf
    """
    merger = PyPDF2.PdfFileMerger()
    graph = open("%sPham%sGraph.pdf" % (output_dir+phage+ one_or_all, pham_no), "rb")
    text = open('%s%sPham%sText.pdf' % (output_dir, phage +one_or_all, pham_no), 'rb')
    merger.append(fileobj=graph)
    merger.append(fileobj=text)
    if one_or_all == 'One':
        file_path = "%s%sPham%sReport.pdf" % (final_path, phage, pham_no)
    else:
        file_path = "%s%sPham%sReport.pdf" % (output_dir, phage+one_or_all, pham_no)
    merger.write(open(file_path, 'wb'))
    return file_path

def make_phage_report(phage, pham_genes, seq_length):
    """
        Creates a single PDF file for a phage, combines the pham genome graph,
        the suggested starts page, and graph and text output for each pham in the phage
        {Phage}Report.pdf
    """
    merger = PyPDF2.PdfFileMerger()
    phage_starts = open("%sSuggestedStarts.pdf" % (output_dir+ phage), 'rb')
    phage_genome = open('%sPhamsGraph.pdf' % (output_dir+phage), 'rb')
    merger.append(fileobj=phage_genome)
    merger.append(fileobj=phage_starts)

    for gene_num in sorted(pham_genes.keys()):
        pham = pham_genes[gene_num][1]
        graph = open("%sPham%sGraph.pdf"  % (output_dir+ phage + one_or_all, pham), "rb")
        text = open('%sPham%sText.pdf' % (output_dir + phage + one_or_all, pham), 'rb')
        merger.append(fileobj=graph)
        merger.append(fileobj=text)
    merger.write(open("%s%sReport.pdf" % (final_path)))

def make_fasta_file(genes, fasta_file):
    count = SeqIO.write(genes, fasta_file, 'fasta')



if 'graph' in args.make:
    print args.pickle_file
    pham = cPickle.load(open(args.pickle_file.strip('"'), 'rb'))
    graph_start_sites(pham, output_dir)

if 'start' in args.make:
    make_suggested_starts(args.phage)

if 'genome' in args.make:
    make_pham_genome(args.phage)

if 'text' in args.make:
    print args.pickle_file
    pham = cPickle.load(open(args.pickle_file.strip('"'), 'rb'))
    if args.phage == None:
        make_pham_text(pham, args.pham_no, output_dir, only_pham=True)
    else:
        make_pham_text(pham, args.pham_no, output_dir)

if 'fasta' in args.make:
    # pickle_file = args.file
    genes = cPickle.load(open(args.pickle_file.strip('"'), 'rb'))
    make_fasta_file(genes, (output_dir + '.fasta'))