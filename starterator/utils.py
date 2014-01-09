import MySQLdb
import ConfigParser
import getpass
import os, subprocess
import re
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from pkg_resources import resource_string, resource_filename
import shutil

MAKING_FILES = "making_files.py"
ICON_FILE =  "extras/starterator.svg"
CONFIGURATION_FILE = "extras/starterator.config"
DESKTOP_FILE = "extras/starterator.desktop"
HELP_FILES = "Help"
config_file = os.path.abspath(os.path.join(os.environ["HOME"], ".starterator/starterator.config"))
icon_file = os.path.abspath(os.path.join(os.environ["HOME"], ".starterator/starterator.svg"))
desktop_file = os.path.abspath(os.path.join(os.environ["HOME"], ".local/share/applications/", "startertor.desktop"))
help_files = os.path.join(os.environ["HOME"], ".starterator", "Help")

def get_pham_colors(db):
    """
        Function to save the colors of the phams in a dictionary mapping pham name to pham color
    """
    cursor = db.cursor()
    cursor.execute("SELECT `name`, `color` FROM `pham_color`;")
    results = cursor.fetchall()
    pham_colors_ = {}
    for row in results:
        pham_colors_[str(row[0])] = row[1]

    return pham_colors_

def get_pham_no(db, phage_name, gene_number):
    """
        Function that gets the phamily number/name of a gene,
        given a geneID from the Mycobacteriophage protein database
        returns a pham number as a string
    """

    # try:
    cursor = db.cursor()
    #Query the database for genes that match the given gene number and phage name
    like_phage_name = phage_name +'%'
    print like_phage_name, gene_number
    cursor.execute("SELECT `gene`.`GeneID` , `pham`.`name` , `phage`.`PhageID`\n\
    FROM `gene`\n\
    JOIN `pham` ON `gene`.`GeneID` = `pham`.`GeneID`\n\
    JOIN `phage` ON `gene`.`PhageID` = `phage`.`PhageID`\n\
    WHERE `phage`.`Name` LIKE %s \n\
    AND `gene`.`GeneID` LIKE %s \n\
    ESCAPE '!'", (like_phage_name, '%!_'+ str(gene_number)))
    results = cursor.fetchall()
    print results, phage_name, gene_number
    # There should only be one result.
    row = results[0]
    pham_no = row[1]
    # print 'pham_no', pham_no
    # Returns the pham number as a string
    return str(pham_no)
    # except:
    #     pass


def find_phams_of_a_phage(db, phage):
    """
        Given a phage in Phamerator, function finds the phams of the genes
        of that phage.
        Returns a list of phams and the length of the nucleotide sequence of the phage 
    """
    # write sql statement to get phams of a phage
    cursor = db.cursor()
    cursor.execute("""SELECT `pham`.`GeneID`, `pham`.`name`, `phage`.`Name`, `phage`.`SequenceLength`\n\
    from `pham` \n\
    join `gene` on `pham`.`GeneID` = `gene`.`GeneID` \n\
    join `phage` on `phage`.`PhageID` = `gene`.`PhageID`\n\
    where `phage`.`Name` = %s""", (phage))
    results = cursor.fetchall()
    phage_phams = []
    seq_length = results[0][3]
    for row in results:
        print row[0], row[1]
        phage_phams.append([row[0],str(row[1])])
    return phage_phams, seq_length


def get_protein_sequences(db):
    proteins = []
    cursor = db.cursor()
    cursor.execute('SELECT GeneID, translation from gene')
    results = cursor.fetchall()
    for row in results:
        protein = SeqRecord(Seq(row[1], IUPAC.protein), id=row[0]+"+", name=row[0], description=row[0])
        proteins.append(protein)
    return proteins

def update_protein_db(db, config_info):
    proteins = get_protein_sequences(db)
    fasta_file = config_info['protein_db'] + 'Proteins'
    count = SeqIO.write(proteins, fasta_file, 'fasta')
    if not config_info['legacy_blast']:
        blast_db_command = [config_info['blast_dir'] + 'makeblastdb',
                    '-in',"\""+ fasta_file+"\"",
                    "-dbtype","prot", "-title", "Proteins",
                     "-out", "%s"% fasta_file]
        print blast_db_command
    else:
        blast_db_command = [config_info['blast_dir'] + 'formatdb',
                    '-i', "\""+ fasta_file+ "\"",
                    '-o', 'T',
                    "-t", "Proteins"]
        print blast_db_command
    subprocess.check_call(blast_db_command)

def check_protein_db(db):
    cursor = db.cursor()
    cursor.execute('SELECT count(*) from gene')
    results = cursor.fetchall()
    new_count = results[0][0]
    return new_count

def get_phage_seq_from_db(db, phage):
    try:
        cursor = db.cursor()
        cursor.execute("SELECT `Sequence`\n\
        from `phage`\n\
        where `Name` like %s", (phage))
        results = cursor.fetchall()
        return results[0][0]
    except:
        pass

def get_gene_number(geneID):
    """
        Get the gene number as an int given a geneID
    """
    # print geneID
    match = re.search(r'^(\w+)([_-]*\w*)_([a-zA-Z]*)(\d+)+$', geneID)
    gene_num = match.groups()[-1]
    # print 'regular exp:', gene_num
    return int(gene_num)

def desktop_file():
    # if not os.path.exists(os.path.join(os.environ["HOME"], ".starterator", "starterator.svg")):
    shutil.copyfile(ICON_FILE,
            os.path.join(os.environ["HOME"], ".local/share/applications/", "starterator.svg"))
    # if not os.path.exists(os.path.join(os.environ["HOME"], ".local/share/applications/", "starterator.desktop")):
    shutil.copyfile(DESKTOP_FILE,
            os.path.join(os.environ["HOME"], ".local/share/applications/", "starterator.desktop"))

def create_folders():
    if not os.path.exists(os.path.join(os.environ["HOME"], ".starterator")):
        os.mkdir(os.path.join(os.environ["HOME"], ".starterator"))
    if not os.path.exists(os.path.join(os.environ["HOME"], ".starterator", "starterator.config")):
        shutil.copyfile("extras/starterator.config", )
    if not os.path.exists(os.path.join(os.environ["HOME"], ".starterator", "Intermediate Files")):
        os.mkdir(os.path.join(os.environ["HOME"], ".starterator", "Intermediate Files"))
    if not os.path.exists(os.path.join(os.environ["HOME"], ".starterator", "Report Files")): 
        os.mkdir(os.path.join(os.environ["HOME"], ".starterator", "Report Files"))
    if not os.path.exists(os.path.join(os.environ["HOME"], ".starterator", "Proteins")):
        os.mkdir(os.path.join(os.environ["HOME"], ".starterator", "Proteins"))

def config_file():
    shutil.copyfile(CONFIGURATION_FILE,
        os.path.join(os.environ['HOME'], '.starterator', 'starterator.config'))
    config = get_config()
    config["intermediate_file_dir"] = os.path.abspath(os.path.join(os.environ["HOME"], ".starterator", "Intermediate Files"))
    config["final_file_dir"] = os.path.abspath(os.path.join(os.environ["HOME"], ".starterator", "Report Files"))
    config["protein_db"] = os.path.abspath(os.path.join(os.environ["HOME"], ".starterator", "Proteins"))
    write_to_config_file(config)


def set_up():
    create_folders()
    config_file()
    desktop_file()

def write_to_config_file(config_info):
    config = ConfigParser.RawConfigParser()
    config.add_section('Starterator')
    for name in config_info:
        config.set('Starterator', name, config_info[name])
    print 'write'
    with open(config_file, 'wb') as configfile:
        config.write(configfile)

def get_config():
    if not os.path.exists(os.path.join(os.environ["HOME"], ".starterator")):
        set_up()
    config_file = os.path.abspath(os.path.join(os.environ["HOME"], ".starterator/starterator.config"))
    config = ConfigParser.RawConfigParser()
    config.read(config_file)
    print "?", CONFIGURATION_FILE, config
    config_info = dict(config.items('Starterator'))
    return config_info
    # config_info['intermediate_file_dir'] = os.path.abspath(self.config_info['intermediate_file_dir'])+ '/'
    # config_info['final_file_dir'] = os.path.abspath(self.config_info['final_file_dir']) + '/'
    # self.config_info['protein_db'] = os.path.abspath(self.config_info['protein_db']) + '/'


def db_connect(config_info):
    db = MySQLdb.connect(config_info['database_server'], 
            config_info['database_user'],
            config_info['database_password'],
            config_info['database_name'])
    return db
    
def attempt_db_connect(config_info):
    try:
        print 'attempting to connect', config_info
        db = MySQLdb.connect(config_info['database_server'], 
                config_info['database_user'],
                config_info['database_password'],
                config_info['database_name'])
        db.close()
    except:
        config_info['database_server'] = raw_input("Enter Database server: ")
        config_info['database_user'] = raw_input("Enter Database username: ")
        config_info['database_password'] = getpass.getpass('Enter Database password: ')
        config_info['database_name'] = raw_input("Enter database name: ")
        print config_info
        attempt_db_connect(config_info)
    return config_info

def clean_up_files(file_dir):
    for f in os.listdir(file_dir):
        os.remove(os.path.join(file_dir, f))

