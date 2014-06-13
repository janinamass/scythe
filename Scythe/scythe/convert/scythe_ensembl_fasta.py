import httplib2
import sys
import os
from ftplib import FTP
import gzip
from helpers.fastahelper import FastaParser

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

if __name__ == "__main__":
    main()
