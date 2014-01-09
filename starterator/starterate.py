#!/usr/bin/env python

from phams import Pham
import phams
import argparse
import utils
from phages import Phage
from genes import Gene
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter, A4
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
import PyPDF2
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter, A4
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline as Blastp
from Bio.Blast.Applications import BlastallCommandline
import MySQLdb
import subprocess
import os, sys
import math
import cPickle
from gi.repository import Gtk, Gdk, GObject
import getpass

def gui():
    GObject.threads_init()
    Gdk.threads_init()
    win = StarteratorWindow()
    win.connect('delete-event', Gtk.main_quit)
    win.show_all()
    Gdk.threads_enter()
    Gtk.main()
    Gdk.threads_leave()

def get_output_one_pham(pham, pham_no, config):
    """
        Creates a PDF Report for the specific pham.
        From Start sites statisitics
    """
    output_dir = config['intermediate_file_dir']
    doc = SimpleDocTemplate("%s%sPham%sText.pdf" % (output_dir, phage+one_or_all, pham_no), pagesize=letter)
    story = []
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="paragraph"))
    styles.add(ParagraphStyle(name='Center', alignment=TA_CENTER))
    text = '<font size=14> Pham %s Report </font>' % pham_no
    story.append(Paragraph(text, styles['Center']))
    story.append(Spacer(1, 12))
    output = pham.output_start_sites()
    for line in output:
        if line == '':
            story.append(Spacer(1, 12))
        text = '<font size=12> %s </font>' % line
        story.append(Paragraph(text, styles['Normal']))
    suggested_start = pham.output_suggested_starts()
    story.append(Spacer(1, 12))
    for line in suggested_start:
        text = '<font size=12>%s</font>' % line
        story.append(Paragraph(text, styles["Normal"]))
    doc.build(story)

def get_arguments():
    parser = argparse.ArgumentParser(prog='starterate.py', usage='Phameratored Phage Report: %(prog)s -p {Phage Name}\n'
            +'One Gene of Phameratored Phage Report:  %(prog)s -p {Phage Name} -n {Pham Number}\n'
            + 'Unphameratored Phage Report:  %(prog)s -p {Phage Name} -u True -f {Path to DNAMaster profile file}\n'
            + 'One Gene of Unphameratored Phage Report:  %(prog)s -p {Phage Name} -u True -s {Start of Gene} -t {Stop of Gene} -o {Orientation of Gene} -g {Number of Gene}')
    parser.add_argument('-n', '--pham_no', default = -1, help='Number of the Pham. For case when want report of a phameratored phage gene.')
    parser.add_argument('-p' , '--phage', default=None, help='The Phamerator database Phage Name. Always needed')
    parser.add_argument('-u', '--unphamed', type=bool, default=False,
                 help='Boolean. If phage has been phameratored: False.'
                +' If phage is unphameratored: True. For use when want report with an unphameratored phage')
    parser.add_argument('-s', '--given_start', type=int, default=-1, 
        help= 'The start of a gene. For case when want report of one gene of an unphameratored phage. ')
    parser.add_argument('-t', '--given_stop', type=int, 
        help= 'The stop of a gene. For case when want report of one gene of an unphameratored phage.')
    parser.add_argument('-o', '--given_orientation',
         help='The orientation of a gene. For case when want report of one gene of an unphameratored phage.')
    parser.add_argument('-g', '--gene_number', default=-1,
         help='The number of a gene. For case when want report of one gene of an unphameratored phage.')
    parser.add_argument('-d', '--profile', help='Path to a DNAMaster profile. For case when want whole report of an unphameratored phage')
    parser.add_argument('-f', '--fasta', help='Path to Fasta File')
    return parser.parse_args()

class UnphameratedGene(object):
    """
        Class to hold the protein sequence record and
        gene sequence record of genes that have not been
        entered into Phamerator
    """
    def __init__(self, gene, protein):
        self.gene = gene
        self.protein = protein

def get_final_file(path_to_file):
    try:
        if os.path.exists(path_to_file) and os.path.isfile(path_to_file):
            return path_to_file, '%sReport.pdf' % phage
        # for root, files, dirs in os.walk(output_dir):
        #     for name in files:
        #         if phage + one_or_all in name:
        #             print 'should wait'
        #             return open_file(path_to_file)
        if args.unphamed == False:
            return phamerated_phage()
        else:
            return unphamerated_phage()
    except Exception:
        raise


def make_report(pham, output_dir, final_dir):
    one_or_all = 'One'
    phage = None
    pham.put_similar_genes_together()
    pickle_file = output_dir + '%sPham%s.p' % (one_or_all, pham.pham_no)
    cPickle.dump(pham, open(pickle_file, 'wb'))
    print "made pickle file", pickle_file

    subprocess.check_call(['python', utils.MAKING_FILES, 
        '-n', pham.pham_no,
        '-a', one_or_all,
        '-f', "\"%s\"" % pickle_file,
        '-o', "\"%s\"" % output_dir,
        '-d', "\"%s\"" % final_dir,
        '-m', 'text'])
    subprocess.check_call(['python', utils.MAKING_FILES, 
        '-n', pham.pham_no,
        '-a', one_or_all,
        '-f', "\"%s\"" % pickle_file,
        '-o', "\"%s\"" % output_dir,
        '-d', "\"%s\"" % final_dir,
        '-m', 'graph'])
    return merge_report(pham.pham_no, output_dir, final_dir)


def merge_report(pham_no, output_dir, final_dir):
    """
        Merges the graphical and text output of one pham into one PDF file
        {Phage}Pham{Number}Report.pdf
    """
    merger = PyPDF2.PdfFileMerger()
    graph = open("%sPham%sGraph.pdf" % (output_dir + 'One', pham_no), "rb")
    text = open('%s%sPham%sText.pdf' % (output_dir, 'One', pham_no), 'rb')
    merger.append(fileobj=graph)
    merger.append(fileobj=text)
    file_path = "%sPham%sReport.pdf" % (final_dir, pham_no)
    merger.write(open(file_path, 'wb'))
    return file_path, "Pham%sReport.pdf" % (pham_no)

def starterate(db, config, info, gui=None, event=None):
    global one_or_all, phage, protein_db, output_dir, final_dir
    one_or_all = 'All' if info['all'] else 'One'
    phage = info['phage']
    # print phage
    protein_db = config['protein_db'] + 'ProteinsDB'
    output_dir = config['intermediate_file_dir']
    final_dir = config['final_file_dir']
    print config["legacy_blast"]
    if info['all'] and info['phamerated']:
        phage = Phage(phage, info['phamerated'], config, db, gui=gui, event=event)
        final_file, short_final = phage.make_report()
    elif info['all'] and not info['phamerated']:
        # phams.update_protein_db(db, config)
        phage = Phage(phage, info['phamerated'], config, db, fasta_file=info['fasta'], profile_file=info['profile'], gui=gui, event=event)
        final_file, short_final = phage.make_report()
    elif not info['all'] and info['phamerated'] and not info['pham']:
        gene = Gene(info['gene_no'], info['phage'], info['phamerated'], config, db)
        print gene
        gene.make_report()
        final_file, s = gene.one_pham_report()
    elif not info['all'] and not info['phamerated'] and not info['pham']:
        # phams.update_protein_db(db, config)
        gene = Gene(info['gene_no'], info['phage'], info['phamerated'], config, db, fasta_file=info['fasta'])
        gene.make_unphamerated_gene(int(info['start']), int(info['stop']), info['orientation'])
        gene.blast_gene()
        print gene
        gene.make_report()
        final_file, s = gene.one_pham_report()
    else:
        # Pham without phage'
        pham = Pham(info['pham'], db, None)
        # self, pham_no, db, phage_name)
        pham.align(config['intermediate_file_dir'], "One")
        final_file, s = make_report(pham, config['intermediate_file_dir'], config['final_file_dir'])

    return final_file



def main():
    args = get_arguments()
    config = utils.get_config()
    db = utils.db_connect(config)
    # config_info['intermediate_file_dir'] = "/home/marissa/"
    # config_info['final_file_dir'] = raw_input('Enter final_file_dir')
    # --Phamerated and only one gene
    if args.gene_number != -1 and args.phage != None and args.unphamed == False:
        gene = Gene(args.gene_number, args.phage, True, config, db)
        print gene
        gene.make_report()
        final_file, s = gene.one_pham_report()

    # --Unphameratored Phage with only one gene
    elif args.given_start > -1 and args.phage != None and args.unphamed == True:
        # given start and stop coordinates and orientation
        one_or_all = 'One'
        given_start = args.given_start
        given_stop = args.given_stop
        given_orientation = args.given_orientation
        gene_name = args.phage + '_' + args.gene_number
        phams.update_protein_db(db, config)
        gene = Gene(args.gene_number, args.phage, args.unphamed, config, db, fasta_file=args.fasta)
        gene.make_unphamerated_gene(given_start, given_stop, given_orientation)
        gene.blast_gene()
        print gene
        gene.make_report()
        final_file, s = gene.one_pham_report()

    # --Phameratored or Unphameratored Phages with all genes
    elif args.pham_no == -1 and args.phage != None and args.unphamed == False:
        phage = Phage(args.phage, True, config, db, gui=None)
        final_file, short_final = phage.make_report()

    elif args.pham_no == -1 and args.phage != None and args.unphamed == True:
        phage = Phage(args.phage, False, config, db, fasta_file=args.fasta, profile_file=args.profile, gui=None)
        final_file, short_final = phage.make_report()

    # clean_up_files()
    # email_final_report(args.email, short_final)


if __name__ == '__main__':
    try:
        print utils.DIR_PATH
        main() 
    except Exception, e:
        # send_error_email(e)
        raise
    finally:
        pass
         # clean_up_files()
