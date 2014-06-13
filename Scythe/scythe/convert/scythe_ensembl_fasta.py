#!/usr/bin/env python

import httplib2
import sys
import os
from ftplib import FTP
import gzip
import getopt

class FastaParser(object):
    def read_fasta(self, fasta, delim = None, asID = 0):
        """read from fasta fasta file 'fasta'
        and split sequence id at 'delim' (if set)\n
        example:\n
        >idpart1|idpart2\n
        ATGTGA\n
        and 'delim="|"' returns ("idpart1", "ATGTGA")
        """
        name = ""
        fasta = open(fasta, "r")
        while True:
            line = name or fasta.readline()
            if not line:
                break
            seq = []
            while True:
                name = fasta.readline()
                name = name.rstrip()
                if not name or name.startswith(">"):
                    break
                else:
                    seq.append(name)
            joinedSeq = "".join(seq)
            line = line[1:]
            if delim:
                line = line.split(delim)[asID]
            yield (line.rstrip(), joinedSeq.rstrip())
        fasta.close()


def getSequencesFromFTP(outdir, release, specieslist=[]):
    path=outdir+os.sep+"fa"
    print(outdir, release, specieslist)
    if len(release)==1:
        release = release[0]
        dirlist = []
        if not os.path.isdir(path):
                os.makedirs(path)
        if (specieslist==[]):
            sys.exit("Error: No set of species selected")

        ftp = FTP('ftp.ensembl.org')
        ftp.login()
        ftp.cwd('pub') #ftp://ftp.ensembl.org/pub/
        ftp.retrlines('LIST', callback=dirlist.append)           # list directory contents
        dirlist = [r for r in dirlist if "release-"+str(release)  in r and "fasta" in r]
        if(len(dirlist )==0):
            sys.exit("Nothing available for release-"+str(release))

        ftprelhome = dirlist[0].split(" ")[-1]
        ftp.cwd(ftprelhome)
        ftp.retrlines('LIST', callback=dirlist.append)
        #pep
        for s in specieslist:
            tmp = [d for d in dirlist if s in d][0]
            spec = tmp.split(" ")[-1]
            ftp.cwd(spec)
            ftp.cwd('pep')
            falist = []
            ftp.retrlines('LIST', callback=falist.append)
            falist = [f for f in falist if "all.fa" in f][0]
            fafile = falist.split(" ")[-1]
            outfaname = spec+".fa.gz"
            outfa = open(outfaname,'wb')
            ftp.retrbinary("RETR "+fafile,outfa.write)
            outfa.close()
            if not os.path.exists(path):
                    os.makedirs(path)
            xtract(outfaname, path)
            ftp.cwd("/pub/"+ftprelhome)
       #nuc
        for s in specieslist:
            tmp = [d for d in dirlist if s in d][0]
            spec = tmp.split(" ")[-1]
            ftp.cwd(spec)
            ftp.cwd('cds')
            falist = []
            ftp.retrlines('LIST', callback=falist.append)
            falist = [f for f in falist if "all.fa" in f][0]
            fafile = falist.split(" ")[-1]
            outfaname = spec+".cds.all.fa.gz"
            outfa = open(outfaname,'wb')
            ftp.retrbinary("RETR "+fafile,outfa.write)
            outfa.close()
            if not os.path.exists(path+"_cds"):
                    os.makedirs(path+"_cds")
            xtract(outfaname, path+"_cds")
            ftp.cwd("/pub/"+ftprelhome)
        ftp.quit()

    else:
        print("Releases:")
        for a in zip(release,specieslist):
            print(a)
        for a in zip(release,specieslist):
            getSequencesFromFTP(outdir, release=[a[0]], specieslist=[a[1]])


def xtract(cfile, outpath = "."):

    if cfile.endswith('.gz'):
        dzf = ".".join(cfile.split(".")[:-1])
        #gzip
        gzf = gzip.open(cfile,'rb')
        content = gzf.read()
        out = open(outpath+os.sep+dzf,'wb')
        out.write(content)
        out.close
        gzf.close()
    else:
        sys.exit("Can't extract "+cfile)

###########################
class Pep(object):
    def __init__(self,pep,gene,length ):
        self.gene = gene
        self.pep  = pep
        self.length = length
        self.isLongest = None

def prepareLocFromFasta(fasta, outpath, specname):
    genes = {}
    longest = {}
    print(fasta, outpath)
    if not os.path.isdir(outpath):
        os.makedirs(outpath)
    if not outpath:
        out = fasta+".loc"
    else:
        out = outpath+os.sep+specname+".loc"
    fp = FastaParser()
    out = open(out,'w')
    for i in fp.read_fasta(fasta):
        tmp = i[0].split(" ")
        geneid = [t for t in tmp if "gene:" in t]
        geneid = geneid[0].split("gene:")[-1]
        proteinid = tmp[0]
        protlen=len(i[1])
        peptide= Pep(proteinid,geneid,protlen)
        try:
            genes[geneid].append(peptide)
            if peptide.length > longest[geneid].length:
                longest[geneid].isLongest =False
                longest[geneid]=peptide
                peptide.isLongest = True
        except KeyError as e:
            genes[geneid] = [peptide]
            longest[geneid] = peptide
            peptide.isLongest=True
    for g in genes:
        print(g)
        firstcol = [w for w in genes[g] if w.isLongest==True ]
        restcol = [w for w in genes[g] if w.isLongest!=True ]
        s = firstcol[0].gene+"\t"+firstcol[0].pep+"\t"
        s+="\t".join([v.pep for v in restcol])
        s+="\n"
        out.write(s)
    out.close()
    print("done")
#######################



def usage():
    """Usage"""
    print("""

    usage: scythe_ensembl_fasta.py -s species1,species2 -r INT

    options:
    -s, --species=STR   comma-separated list of species (eg 'homo_sapiens,gorilla_gorilla')
    -r, --release=NUM   ENSEMBL version (eg '75')
    -d, --dir DIR       output directory [default ./]
    -h, --help          prints this
    """)
    sys.exit(2)

def main():
    specList = None
    release = None
    outdir = "./"
    try:
        opts, args = getopt.gnu_getopt(sys.argv[1:], "s:r:h",["species=","release=","help"])
    except getopt.GetoptError as err:
        sys.stderr.write(str(err))
        usage()
    for o, a in opts:
        if o in ("-s", "--species"):
            specList = a.split(",")
        elif o in ("-r", "--release"):
            release = int(a)
        elif o in ("-h", "--help"):
            usage()
        else:
            assert False, "unhandled option"

    if not (specList):
        usage()
    if not release:
        usage()

    getSequencesFromFTP(outdir=outdir, release=[release], specieslist=specList)

if __name__=="__main__":
    main()

