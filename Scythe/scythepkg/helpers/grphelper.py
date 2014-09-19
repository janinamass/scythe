class GrpParser(object):
    """read .loc files. Return dictionary with locus ID as key and transcript list as value"""
    def readLoc(self, locfiles, speciesNames = None):
        firstLine=True
        res = {}
        for lf in locfiles:
            tmpfile = open(lf,'r')
            for ln in tmpfile:
                if firstLine:
                    firstLine = False
                ln = ln.rstrip()
                tmp = ln.split("\t")
                k,v = str(tmp[0]),tmp[1:]
                res[k]=v
            tmpfile.close()
        return(res)

    def readGroups(self, grpf, locf, speciesNames = None):
            """Iterate over grp return numbered list with orthologous transcripts."""
            grp = open(grpf, 'r')
            locDct = self.readLoc(locfiles =locf)
            traDct ={}
            #locDct has loci as keys, lists of transcripts as values
            for l in grp:
                tmp = l.strip().split("\t")
                grpID = tmp[0]
                tmpLoci=tmp[1:]
                try:
                    tr = [locDct[x] for x in tmpLoci]
                except KeyError as ke:
                    print("key error",ke," @readGroups. Locus does not appear to have transcripts.")
                traDct[grpID] = tr
                yield int(grpID), traDct[grpID]
            grp.close()

    def groupDct(self, grpf, locf,speciesNames = None):
            res = {}
            for sth in self.readGroups(grpf, locf,speciesNames):
                res[sth[0]]=sth[1:][0]
            return res

