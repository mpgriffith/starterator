#!/usr/bin/env python

import os
from collections import Counter
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation
import copy
import subprocess
import cPickle
from time import sleep
import utils

def find_upstream_stop_site(start, stop, orientation, phage_sequence):
    ahead_of_start = 0
    stop_site_found = False
    stop_codons = ['AGT', 'AAT', 'GAT']
    while not stop_site_found:
        ahead_of_start += 99
        if orientation == 'R':
            sequence = Seq(phage_sequence[stop-1:(start+ahead_of_start)],
                    IUPAC.unambiguous_dna)
            sequence = sequence.reverse_complement()
            if stop < 400:
                return sequence, ahead_of_start
        else:
            if start < ahead_of_start:
                ahead_of_start = start - start % 3
                sequence = Seq(phage_sequence[(start-ahead_of_start):stop],
                     IUPAC.unambiguous_dna)
                return sequence, ahead_of_start
            if stop < start:
                end_sequence = phage_sequence[(start-ahead_of_start):]
                start_sequence = phage_sequence[:stop]
                sequence = Seq(end_sequence+start_sequence, IUPAC.unambiguous_dna)
            else:
                sequence = Seq(phage_sequence[(start-ahead_of_start):stop],
                 IUPAC.unambiguous_dna)
        sequence_ahead_of_start = sequence[:ahead_of_start]
        sequence_ahead_of_start = sequence_ahead_of_start[::-1]
        
        for index in xrange(0, len(sequence_ahead_of_start), 3):
            codon = str(sequence_ahead_of_start[index:index+3])
            if codon in stop_codons:
                new_ahead_of_start = index
                new_sequence = sequence[(ahead_of_start - index):]
                return new_sequence, new_ahead_of_start

def make_gene(geneID, start, stop, orientation, phage_sequence, phage_name):
    """
        Make a SeqRecord of the given gene
    """
    ahead_of_start = 100
    if orientation == 'R':
        temp_start = stop
        stop = start
        start = temp_start
    sequence, ahead_of_start = find_upstream_stop_site(start, stop, orientation, phage_sequence)
    gene = SeqRecord(sequence, id=geneID, name=phage_name)
    gene.annotations["start"] = start
    gene.annotations['stop'] = stop
    gene.annotations["orientation"] = orientation
    gene.annotations['number_before_start'] = ahead_of_start
    gene.annotations['translation'] = sequence.translate()
    return gene

def get_phage_names_seqs(phage_sequences):
    """
        make a dictionary of phage names mapped to phage phageID
        and a dictionary of phage sequences mapped to phage names
    """
    phage_seq = open(phage_sequences, 'r')
    phage_seqs = {}
    phage_names = {}
    next(phage_seq)
    for line in phage_seq:
        temp = line.split('\t')
        phageID = temp[2].strip()
        phage_sequence = temp[1].strip()
        phage_name = temp[0].strip()
        phage_seqs[phage_name] = phage_sequence
        phage_names[phageID] = phage_name
    phage_seq.close()
    return phage_names

class Pham(object):
    """
        Class for Pham
        Arguments for initialization:
            pham_no     :   the number of the pham
            phage_name  :   name of the phage
            db          :   connection to phagesdb.org database
                                (used for obtaining the nucleotide sequence of each gene in a pham)  
        Attributes:
            pham_no                 :   the number of the pham
            phage                   :   the name of phage
            phage_gene              :   list of genes of phage
            genes_in_pham           :   dict of genes with geneID mapping to tuple of gene and aligned gene 
            total_possible_starts   :   list of total possible starts, number corresponds to alignment position
            most_called_start       :   alignment index of most called start
        """
    def __init__(self, pham_no, db, phage_name, root_file_name=None):
        self.pham_no = pham_no
        print pham_no, type(pham_no)
        self.phage = phage_name
        self.phage_gene = []
        self.genes_in_pham = {}
        if root_file_name == None:
            self.file_name = "Pham%s" % (self.pham_no)
        else: 
            self.file_name = root_file_name
        self.load_pham_file(db)
        self.total_possible_starts = []
        self.most_called_start = None
        self.called_start_sites = None
        self.start_stats = None

    def _find_candidate_starts(self, gene_sequence):
        """
            Finds all the possible start site of a gene and returns a list of indexes of start sites
        """
        starts = []
        start_codons = ['ATG', 'GTG', 'TTG']
        for index in xrange(0, len(gene_sequence), 3):
            codon = str(gene_sequence[index:index+3])
            if codon in start_codons:
                starts.append(index)
        return sorted(starts)
    
    def find_aligned_candidate_starts(self, gene_record, alignment_record):
        starts = gene_record.annotations['possible_starts']
        count = -1
        aligned_starts = []
        for index, char in enumerate(alignment_record.seq):
            if char != '-':
                count += 1
            if count in starts and char != '-':
                aligned_starts.append(index)
        return aligned_starts

    def _normalize_gene_id(self, geneID, phage_name):
        """
            Makes the geneID in (mostly) the same format for each gene (PhageName_GeneNo, ex Doom_47) 
        """
        if ' ' in geneID:
           geneID = geneID.replace(' ', '_')
        gene = geneID.split("_")
        if len(gene) > 2:
            gene_no = gene[-1]
        elif len(gene) < 2:
            gene_no = gene[0]
        else:
            gene_no = gene[1]
        new_geneID = phage_name +"_"+ gene_no
        return new_geneID

    def change_amount_ahead_of_gene(self, gene, phage_sequence, ahead_of_start):
        """
            When the length of genes in a pham vary significantly, this is called..
            changes the length of the sequence obtained before the called start of the gene
        """
        start = gene.annotations['start']
        stop = gene.annotations['stop']
        if ahead_of_start > start:
            ahead_of_start = start - 1
        if gene.annotations['orientation'] == 'R':
            sequence = Seq(phage_sequence[stop-1:(start+ahead_of_start)], IUPAC.unambiguous_dna)
            new_sequence = sequence.reverse_complement()
        else:
            new_sequence = Seq(phage_sequence[(start-ahead_of_start):stop], IUPAC.unambiguous_dna)
        gene.annotations['number_before_start'] = ahead_of_start
        gene.seq = new_sequence

    def get_pham_genes_from_db(self, db, pham_no):
        try:
            cursor = db.cursor()
            cursor.execute(
            """SELECT `gene`.`GeneID`,`phage`.`Name`, 
            `Length`, `Start`, `Stop`, `Orientation`
            FROM `gene`
            JOIN `pham` ON `gene`.`GeneID` = `pham`.`GeneID`
            JOIN `phage` ON `gene`.`PhageID` = `phage`.`PhageID`
            WHERE `pham`.`name` =%s; """ ,(pham_no))
            return cursor.fetchall()
        except:
            pass

    def load_pham_file(self, db):
        """
            Takes in the phage_seqs dictionary and the path to the file and returns a list of genes.
            Each gene is a SeqRecord that contains the sequence of the gene beginning 100 before the start site,
            the phage name, the geneID, the location of start and stop coordinates, the orientation, the translation,
            and the candidate starts of the gene
        """
        print db
        genes_from_sql = self.get_pham_genes_from_db(db, self.pham_no)
        print type(genes_from_sql)
        for row in genes_from_sql:
            phage_name = row[1]
            geneID = self._normalize_gene_id(row[0], phage_name)
            length = int(row[2])
            start = int(row[3])
            stop = int(row[4])
            orientation = row[5]
            ahead_of_start = 100
            phage_sequence = utils.get_phage_seq_from_db(db, phage_name)
            gene = make_gene(geneID, start, stop, orientation, 
                            phage_sequence, phage_name)
            gene.annotations["possible_starts"] = self._find_candidate_starts(gene.seq)
            if self.phage == phage_name:
                self.phage_gene.append(gene)
            self.genes_in_pham[geneID] = gene
        return self.genes_in_pham

    def add_alignment(self, alignment):
        for record in alignment:
            gene = self.genes_in_pham[record.id]
            self.genes_in_pham[record.id] = gene, record

    def add_gene(self, gene):
        """
            If genome is not in phamerator, this function will allow it to be added to a pham
        """
        gene.annotations["possible_starts"] = self._find_candidate_starts(gene.seq)
        self.phage_gene.append(gene)

        self.genes_in_pham[gene.id] = gene

    def _make_gaps_features(self):
        """
            Makes gaps and sequences features for each record in the alignment_record
            Useful for the graphical output only
        """
        for gene in self.genes_in_pham.itervalues():
            alignment_record = gene[1]
            gap = False
            gap_count = 0
            seq_count = 0
            for index, char in enumerate(alignment_record.seq):
                if char == '-' and gap == False:
                    gap = True
                    gap_count = 0
                    if seq_count > 0:
                        seq_feature = SeqFeature(FeatureLocation(index-seq_count, index-1),
                            type='seq', strand=None)
                        alignment_record.features.append(seq_feature)
                elif char == '-' and gap == True:
                    gap_count += 1
                elif char != '-' and gap == True:
                    gap = False
                    seq_count = 0
                    gap_feature = SeqFeature(FeatureLocation(index-gap_count-1, index-1),
                            type='gap', strand=None)
                    alignment_record.features.append(gap_feature)
                else: #gap = False and char != '-'
                    seq_count +=1
            if gap == True:
                gap_feature = SeqFeature(FeatureLocation(index-gap_count-1, index), 
                            type='gap', strand=None) 
                alignment_record.features.append(gap_feature)
            else:
                seq_feature = SeqFeature(FeatureLocation(index-seq_count, index), 
                            type='seq', strand=None)
                alignment_record.features.append(seq_feature)
    
    def _add_number_before_start_alignment(self):
        """
            adds the number of nucleotides before start to the alignment_record
            annotations['number_before_start']
        """
        for gene in self.genes_in_pham.itervalues():
            number = gene[0].annotations['number_before_start']
            gene[1].annotations['number_before_start'] = number

    def _find_start_site_alignment(self):
        """
            Finds the position of the start site on the sequence for each alignment record
        """
        for record in self.genes_in_pham.itervalues():
            aligned_record = record[1]
            count = 0
            i = aligned_record.annotations['number_before_start']
            for index, letter in enumerate(aligned_record.seq):
                if letter in ['A', 'G', 'T', 'C']:
                    count += 1
                    i = index
                if count > aligned_record.annotations['number_before_start']:
                    break
            aligned_record.annotations['start'] = i 

    def _possible_starts_alignment(self):
        """
            For each record in the alignment, makes a list of the candidate start sites
            Currently maps to the index of the alignment sequence
            Stored in the annotations['possible_starts']
        """
        for record in self.genes_in_pham.itervalues():
            gene = record[0]
            aligned_gene = record[1]
            aligned_gene.annotations['possible_starts'] = self.find_aligned_candidate_starts(gene, aligned_gene)            

    def _total_possible_starts(self):
        """
            Finds all the possible starts in alignment
        """
        for record in self.genes_in_pham.itervalues():
            aligned_gene = record[1]
            for site in aligned_gene.annotations['possible_starts']:
                if site not in self.total_possible_starts:
                    self.total_possible_starts.append(site)
        self.total_possible_starts = sorted(self.total_possible_starts)

    def coordinate_start_site(self, index, sequence_record, alignment_record):
        """
            Given an index of the alignment and the sequence record and alignment record,
            finds the coordinates of the index on the phage sequence.
            The alignment record and sequence record should have the same id 
            (in other words, should be the same gene for alignment and original gene record)
        """
        current_start_coordinate = sequence_record.annotations["start"]
        alignment_current_start_coord = alignment_record.annotations["start"]
        aligned_seq = alignment_record.seq
        new_start_index = 0
        for i in xrange(0, index):
            if aligned_seq[i] != '-':
                new_start_index += 1
        if sequence_record.annotations["orientation"] == 'R':
            new_start_coords = (current_start_coordinate + 
                alignment_record.annotations['number_before_start'] - new_start_index)
        else:
            new_start_coords = (current_start_coordinate - 
                alignment_record.annotations['number_before_start'] + new_start_index + 1)
        return new_start_coords
   
    def is_equal(self, alignment_1, alignment_2):
        if alignment_1.annotations['start'] != alignment_2.annotations['start']:
            return False
        if alignment_1.annotations['number_before_start'] != alignment_2.annotations['number_before_start']:
            return False
        possible_starts_1 = alignment_1.annotations['possible_starts']
        possible_starts_2 = alignment_2.annotations['possible_starts']
    
        if set(possible_starts_1) != set(possible_starts_2):
            return False
        if len(alignment_1.features) != len(alignment_2.features):
            return False
        alignment_1.features.sort()
        alignment_2.features.sort()
        for feature1, feature2 in zip(alignment_1.features, alignment_2.features):
            if feature1.location.start != feature2.location.start:
                return False
            if feature1.location.end != feature2.location.end:
                return False
            if feature1.type != feature2.type:
                return False
        return True

    def put_similar_genes_together(self):
        genes_list = [x[1][1] for x in self.genes_in_pham.iteritems()]
        genes_list.sort(key=lambda x: x[0])
        i = 0
        groups = []
        grouped = [False for x in genes_list]
        while i < len(genes_list):
            if not grouped[i]:
                alignment_1 = genes_list[i]
                j = i+1
                group = []
                group.append(alignment_1)
                while j <len(genes_list):
                    if not grouped[j]:
                        alignment_2 = genes_list[j]
                        if self.is_equal(alignment_1, alignment_2):
                            grouped[i] = True
                            grouped[j] = True
                            group.append(alignment_2)
                    j += 1
                groups.append(group)
            i += 1


        self.groups = groups
        return self.groups


    def list_of_candidate_starts(self, gene_id):
        gene_record, alignment_record = self.genes_in_pham[gene_id]
        possible_starts_coords = []
        for start in alignment_record.annotations['possible_starts']:
            graph_number = self.total_possible_starts.index(start)+1
            start_coord = self.coordinate_start_site(start, gene_record, alignment_record)
            possible_starts_coords.append((graph_number, start_coord))
        return possible_starts_coords

    def candidate_starts_phage_genes(self):
        phage_genes_candidate_starts = {}
        for gene in self.phage_gene:
            print gene.id
            starts = self.list_of_candidate_starts(gene.id)
            phage_genes_candidate_starts[gene.id] = starts
        return phage_genes_candidate_starts

    def suggest_start_site_coord(self, gene_id):
        """
            Given a gene, finds the coordinate of the most called start if it
            is contained in the phage. Otherwise, returns a list of tuples of all
            start sites and number label.
        """
        most_called_start_aligned_index = self.total_possible_starts[self.most_called_start-1]
        gene_record, alignment_record = self.genes_in_pham[gene_id]
        old_start_coord = gene_record.annotations["start"]
        if gene_record.id in self.start_stats[0]:
            if gene_record.annotations['orientation'] == 'R':
                return (self.most_called_start, old_start_coord)
            else:
                return (self.most_called_start, old_start_coord + 1)
        elif gene_record.id in self.start_stats[1]:
            new_start_coord = self.coordinate_start_site(
                most_called_start_aligned_index, gene_record, alignment_record)
            return (self.most_called_start, new_start_coord)
        else:
            possible_starts = alignment_record.annotations["possible_starts"]
            possible_starts_coords = []
            for start in possible_starts:
                graph_number = self.total_possible_starts.index(start)+1
                new_start_coord = self.coordinate_start_site(start, 
                    gene_record, alignment_record)
                possible_starts_coords.append((graph_number, new_start_coord))
            return possible_starts_coords
    
    def start_sites_stats(self):
        """
            Returns a list of: 
            genes_with_most_called          : list of genes that have the most called start
            genes_have_most_called_start_not: list of genes that have the most called start but not called   
            genes_without_most_called_start : list of genes without the most called start
            genes_start_sites               : dict mapping a start to list of genes with that start called
            genes_possible_starts           : dict mapping a candidate start to a list of genes with the start

        """
        total_possible_starts_set = set(self.total_possible_starts)
        all_start_sites = [record[1].annotations['start'] 
                            for record in self.genes_in_pham.itervalues()]
        all_starts_set = set(all_start_sites)
        genes_start_sites = {}
        genes_possible_starts = {}
        # print self.total_possible_starts
        # print all_starts_set
        # for record in self.alignment:
        #     record_start = record.annotations['start']
        #     print record.id, record_start
        #     print record.annotations['possible_starts']
        #     print record.seq[record_start: record_start + 9]
        # makes dictionary of index of start sites mapped to list of genes with that candidate start site
        # make dictionary of index of start site mapped to list of genes with that start site called
        for site in total_possible_starts_set:
            genes_possible_starts[self.total_possible_starts.index(site)+1] = [record[1].id 
                for record in self.genes_in_pham.itervalues() if site in record[1].annotations["possible_starts"]]

        for site in all_starts_set:
            if site in self.total_possible_starts:
                genes_start_sites[self.total_possible_starts.index(site)+1] = [record[1].id 
                    for record in self.genes_in_pham.itervalues() if site == record[1].annotations['start']]
        all_starts_count = Counter(all_start_sites)
        self.called_starts_count = all_starts_count.most_common()
        self.most_called_start = self.total_possible_starts.index(self.called_starts_count[0][0])+1
        genes_with_most_called = genes_start_sites[self.most_called_start]
        genes_have_most_called_start_not = []
        genes_without_most_called_start = []
        for gene in self.genes_in_pham.itervalues():
            aligned_gene = gene[1]
            if (aligned_gene.id in genes_possible_starts[self.most_called_start] and 
                aligned_gene.id not in genes_start_sites[self.most_called_start]):
                genes_have_most_called_start_not.append(aligned_gene.id)
            elif aligned_gene.id not in genes_possible_starts[self.most_called_start]:
                genes_without_most_called_start.append(aligned_gene.id)
        self.start_stats = [genes_with_most_called, 
                            genes_have_most_called_start_not, 
                            genes_without_most_called_start,
                            genes_start_sites, 
                            genes_possible_starts]
        return self.start_stats

    def output_suggested_starts(self, all_genes=False):
        """
            gives a String list of the output of the suggested start site
            If all_genes is true, finds the suggested start site of all genes in a pham
            Otherwise, finds only the suggested start site of the specific phage gene(s)
        """
        output = []
        if all_genes == False:
            for gene in self.phage_gene:
                coord_suggestion = self.suggest_start_site_coord(gene.id)
                output.append(""+ gene.id + ", " +str(coord_suggestion))
        else:
            for record in self.genes_in_pham.itervalues():
                coord_suggestion = self.suggest_start_site_coord(record[0].id)
                output.append(record[0].id + ", " +str(coord_suggestion))
        return output

    def output_start_sites(self):
        """
            Writes a report of the start sites statistics
            Returns a list of strings with information about
            start sites for this pham  
        """
        output = []
        output.append("Most Called Start: %s (number based on diagram)"
                         % self.most_called_start)
        percent_with_most_called = (float(len(self.start_stats[0])) 
                                    /len(self.genes_in_pham) *100 )
        output.append("Percent with start called: %10.4f%%" 
                        % percent_with_most_called )
        output.append("Genes that call the most called start:")
        s = u'\u2022' + ''
        for gene in self.start_stats[0]:
            s += gene+ ", "
        output.append(s)
        output.append("")
        output.append("Genes that have the most called start but do not call it:")
        s = u'\u2022' + ''
        for gene in self.start_stats[1]:
            s += gene + ", "
        output.append(s)
        output.append('')
        output.append("Genes that do not have the most called start:")
        s = u'\u2022' + ""
        for gene in self.start_stats[2]:
            s += gene + ", "
        output.append(s + '')
        output.append('')
        output.append("Other Starts Called:")
        for start, genes in self.start_stats[3].items():
            if start != self.most_called_start:
                s = ''
                for gene in genes:
                    s += gene + ", "
                output.append(u'\u2022' + str(start) + "\t" +s +'')
                percent = float(len(genes)) / len(self.genes_in_pham) * 100
                output.append("Percent with start called: %10.4f%% \n\t" % percent)
        return output

    def make_fasta_file(self, genes, fasta_file, one_or_all):
        pickle_file = os.path.abspath(fasta_file +'.p')
        print 'python making_files.py -p %s -n %s -a %s -f %s -m fasta' % (self.phage, self.pham_no, one_or_all, pickle_file)
        cPickle.dump(genes, open(pickle_file, 'wb'))
        print "made picke file for fasta", pickle_file
        sleep(20)
        print 'python making_files.py -p %s -n %s -a %s -f %s -m fasta' % (self.phage, self.pham_no, one_or_all, pickle_file)
        subprocess.check_call(['python', utils.MAKING_FILES, 
        '-p', self.phage, 
        '-n', self.pham_no,
        '-a', one_or_all,
        '-f', "\"%s\"" % pickle_file,
        '-o', "\"%s\"" % fasta_file,
        '-m', 'fasta'])
        sleep(1)
   
    def clustalw_alignment(self, fasta_file):
        subprocess.check_call(['clustalw', 
        '-infile=%s' % (fasta_file),
        '-quicktree'])

    def align(self, file_path, one_or_all, already_aligned=False):
        """
            Aligns the pham using ClustalW
            Returns the pham
        """
        print file_path
        file_name = "%s%sgenes" % (file_path, self.file_name)
        print file_name +'.fasta'
        genes = [record for record in self.genes_in_pham.values()]
        count = SeqIO.write(genes, '%s.fasta' % file_name, 'fasta')
        if len(self.genes_in_pham) < 2:
            for record in self.genes_in_pham:
                gene = self.genes_in_pham[record]
                aligned = copy.deepcopy(self.genes_in_pham[record])
                self.genes_in_pham[record] = gene, aligned
        else:
            try: 
                alignment = AlignIO.read(file_name +".aln", "clustal")
            except:
                print 'aligning'
                self.clustalw_alignment(file_name + '.fasta')
                print len(self.genes_in_pham)
                print file_name+'.aln'
                alignment = AlignIO.read(file_name+".aln", "clustal")

            self.add_alignment(alignment)
        self._add_number_before_start_alignment()
        self._find_start_site_alignment()
        self._possible_starts_alignment()
        self._make_gaps_features()
        self._total_possible_starts()
        self.start_sites_stats()