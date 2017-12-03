#!/home/model/MD-SOFTWARE/anaconda/bin/python

import os
import re
import urllib2
import argparse
import MDAnalysis as md


from Bio.Seq import Seq
from Bio import pairwise2
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
from Bio import AlignIO,SeqIO,ExPASy,SwissProt
from Bio.Emboss.Applications import NeedleCommandline

from future import print_function

__doc__ = """
------------------------------------------------------------------

A script to renumber protein resisdues based on uniprot identifier

arguments:

    1) required:
        -pdb: protein.pdb

    2) optional:
        -chain: chain ID
        -id:    uniprot identifier ([A-Z0-9]{6})
        -mutagenesis

author : Remi Cuchillo
version: 1.0
------------------------------------------------------------------
"""

#==========================================================================
# Read inputs
#==========================================================================

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__, epilog=" ")

parser.add_argument("-pdb", type=str, default=None, help='pdb file')
parser.add_argument("-id", type=str, default=None, help='uniprot identifier')
parser.add_argument("-chain", type=str, default=None, help='chain selection')
parser.add_argument("-mutagenesis", help='create a pymol session showing known mutations',action='store_true')

# parsing arguments
args = parser.parse_args()

if args.pdb == None:print(__doc__),quit()
else:pdbname=args.pdb
if args.id == None:print(__doc__),quit()
else:accession=args.id

#==========================================================================
# renumber
#==========================================================================

oneletter = {
'ASP':'D','GLU':'E','ASN':'N','GLN':'Q',
'ARG':'R','LYS':'K','PRO':'P','GLY':'G',
'CYS':'C','THR':'T','SER':'S','MET':'M',
'TRP':'W','PHE':'F','TYR':'Y','HIS':'H',
'ALA':'A','VAL':'V','LEU':'L','ILE':'I',
}

# Retrieve pdb to extract sequence
# Can probably be done with Bio.PDB but being able to use the vmd-like selection algebra is nice
u = md.Universe(pdbname)
if args.chain == None:
    chains = []
    for i in dir(u):
        if len(i) < 2:chains.append(i)
else:chains = [args.chain]
if len(chains) > 1:
    chains_str=' '
    for i in chains:
        chains_str = chains_str+i+' '
    print('ERROR: %i chains detected (%s) -- choose one'%(len(chains),chains_str))
    quit()
pdbseq_str=''.join([oneletter[i] for i in u.select_atoms('protein and segid %s'%chains[0]).residues.resnames])
alnPDBseq=SeqRecord(Seq(pdbseq_str,IUPAC.protein),id=pdbname)
SeqIO.write(alnPDBseq,"%s.fasta"%pdbname,"fasta")

# Retrieve reference sequence
handle = ExPASy.get_sprot_raw(accession)
swissseq = SwissProt.read(handle)
refseq=SeqRecord(Seq(swissseq.sequence,IUPAC.protein),id=accession)
SeqIO.write(refseq, "%s.fasta"%accession,"fasta")

# Do global alignment with needle from EMBOSS, stores entire sequences which makes numbering easier
needle_cli = NeedleCommandline(asequence="%s.fasta"%pdbname,bsequence="%s.fasta"%accession,gapopen=10,gapextend=0.5,outfile="%s_needle.out"%pdbname[:-4])
needle_cli()
aln = AlignIO.read("%s_needle.out"%pdbname[:-4], "emboss")

os.remove("%s.fasta"%pdbname)
os.remove("%s.fasta"%accession)

alnPDBseq = aln[0]
alnREFseq = aln[1]
# Initialize per-letter annotation for pdb sequence record
alnPDBseq.letter_annotations["resnum"]=[None]*len(alnPDBseq)
# Initialize annotation for reference sequence, assume first residue is #1
alnREFseq.letter_annotations["resnum"]=range(1,len(alnREFseq)+1)

# Set new residue numbers in alnPDBseq based on alignment
reslist = [[i,alnREFseq.letter_annotations["resnum"][i]] for i in range(len(alnREFseq)) if alnPDBseq[i] != '-']
for [i,r] in reslist:
    alnPDBseq.letter_annotations["resnum"][i]=r

# Set new residue numbers in the structure
newresnums=[i for i in alnPDBseq.letter_annotations["resnum"][:] if i != None]
proteinAtoms = u.select_atoms('protein and segid %s'%chains[0])

proteinAtoms.residues.set_resnums(newresnums)
newresnums_atoms = []
n=0
resID = proteinAtoms[0].resid
for i in proteinAtoms:
    if i.resid == resID:
        newresnums_atoms.append(newresnums[n])
    else:
        resID = i.resid
        n += 1
        newresnums_atoms.append(newresnums[n])
proteinAtoms.set_resids(newresnums_atoms)


other_chains = u.select_atoms('protein and not segid %s'%chains[0])
ligands      = u.select_atoms('not protein')
if len(other_chains) > 0 and len(ligands) > 0:
    u_merge = md.Merge(proteinAtoms,other_chains,ligands)
elif len(other_chains) > 0:
    u_merge = md.Merge(proteinAtoms,other_chains)
elif len(ligands) > 0:
    u_merge = md.Merge(proteinAtoms,ligands)
else:
    u_merge = md.Merge(proteinAtoms)
with md.Writer("./%s_renumbered.pdb"%pdbname[:-4], multiframe=False, bonds='conect', n_atoms=u_merge.atoms.n_atoms) as PDB:
    PDB.write(u_merge.atoms)

if args.mutagenesis:
#mutagenesis
    resid_mutations = []
    known_mutations = []
    effect          = []
    page  = urllib2.urlopen('http://www.uniprot.org/uniprot/%s'%accession)
    lines = page.readlines()
    mutagenesis    = re.compile('.*Pathology and Biotech.*key=Mutagenesis\">([0-9]+)<.*\">([A-Z] .+[A-Z]): <span property=\"text\">(.+)</span> <span.+')
    for i in lines:
        search = re.search(mutagenesis,i)
        if search != None:
            resid_mutations.append(int(search.group(1)))
            known_mutations.append(str(search.group(2)))
            effect.append(str(search.group(3)))
    os.system('mkdir %s_mutagenesis > /dev/null 2>&1'%pdbname[:4])
    os.chdir('%s_mutagenesis'%pdbname[:4])
    f_o=open('mutagenesis.txt','w')
    for i in xrange(len(resid_mutations)):
        f_o.write('%5i %5s %s\n'%(resid_mutations[i],known_mutations[i],effect[i]))
    f_o.close()
    proteinAtoms = u.select_atoms('protein and segid %s'%chains[0])
    proteinAtoms.residues.set_resnums(newresnums)
    other_chains = u.select_atoms('protein and not segid %s'%chains[0])
    ligands      = u.select_atoms('not protein')
    protein = proteinAtoms.residues
    for i in protein:
        if i.resid not in resid_mutations:
            i.bfactors(0.0)
        else:
            i.bfactors(0.5)
    if len(other_chains) > 0 and len(ligands) > 0:
        u_merge = md.Merge(proteinAtoms,other_chains,ligands)
    elif len(other_chains) > 0:
        u_merge = md.Merge(proteinAtoms,other_chains)
    elif len(ligands) > 0:
        u_merge = md.Merge(proteinAtoms,ligands)
    with md.Writer("./%s_mutagenesis.pdb"%pdbname[:4], multiframe=False, bonds='conect', n_atoms=u_merge.atoms.n_atoms) as PDB:
        PDB.write(u_merge.atoms)
