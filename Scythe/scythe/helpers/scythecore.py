import re, sys, os, subprocess, string, random
from scythe.helpers.fastahelper import FastaParser
import datetime
import time
class ScytheBase(object):
        """The ScytheBase object."""
        def __init__(self, name):
                self._name = name

        @property
        def name(self):
                return self._name
        @name.setter
        def name(self, value):
                self._name = value
        def __repr__(self):
                return '<Scythe class:%s name:%s dict:%s> ' %(self.__class__.__name__, self.name, self.__dict__.keys())

class ScytheError(Exception):
        def __init__(self, msg):
                self._msg = msg
        @property
        def msg(self):
                return self._msg
        def __str__(self):
                return self._msg

class ScytheSeq(ScytheBase):
        """Store a sequence of type 'dna' or 'aa'."""
        def __init__(self, name = None, species = None, sequence = None, type = None, shortName = None, isReference = None, isSingle = None):
                super().__init__(name)
                if shortName:
                        if len(shortName)>10:
                                raise ScytheError("""'shortName's should be no longer than 10 alphanumeric characters to assure compatibility with software from the phylip package. The recommended length is 6 characters to allow for a 3 character species abbreviation plus separator.
                                 """)
                self._shortName = shortName
                if type:
                        if (type != "dna" and type != "aa"):
                                raise ScytheError("'type' should be 'aa' or 'dna'")
                        else:
                                self._type = type
                #dna or aa
                self._species = species
                #identifier for sequence
                self._sequence = sequence
                #actual dna or amino acid sequence
                self._isReference = isReference
                #this sequence is the reference gene model for its locus
                self._isSingle = isSingle
                #this sequence is the only gene model for its locus
                if self._isSingle:
                    self._isReference = True
                #if it's the only sequence it must be the reference
        @property
        def shortName(self):
                return self._shortName
        @shortName.setter
        def shortName(self, value):
                if len(value)>10:
                        raise ScytheError("""'shortName's should be no longer than 10 alphanumeric characters to assure compatibility with software from the phylip package. The recommended length is 6 characters to allow for a 3 character species abbreviation plus separator.""")
                self._shortName = value

        @property
        def type(self):
                return self._type
        @type.setter
        def type(self, value):
                self._type = value

        @property
        def species(self):
                return self._species
        @species.setter
        def species(self, value):
                self._species = value

        @property
        def sequence(self):
                return self._sequence
        @sequence.setter
        def sequence(self, value):
                self._sequence = value


        @property
        def isSingle(self):
            return self._isSingle
        @isSingle.setter
        def isSingle(self,value):
                self._isSingle = value
                if value == True:
                    self.isReference = True
        @property
        def isReference(self):
            return self._isReference

        @isReference.setter
        def isReference(self,value):
                self._isReference = value

        def toFasta(self, specTF = True, shortTF = False, newlines = None, sep = "|"):
                """Turn ScytheSequence into fasta format.
                Either with species name (specTF = True) or without.
                Warning:
                Use w/ specTF=False and shortTF=True only with globally unique IDs.
                """
                #use w/o spec and w/ short form only with globally unique IDs!
                if shortTF and not self._shortName:
                        raise ScytheError("This Sequence object doesn't have short name (.shortName).")
                if specTF and not self._species:
                        raise ScytheError("This Sequence object doesn't have a species name (.species).")
                if not self._sequence:
                        raise ScytheError("This Sequence object doesn't have a sequence (.sequence). ")

                name = self._name
                if specTF == True:
                        res = ">"+self._species+sep
                else:
                        res = ">"
                if shortTF==True:
                                name = self._shortName
                if newlines:
                        tmp = FastaHelper().insert_newlines(self._sequence, every=newlines)
                        return (res+name+"\n"+tmp+"\n")
                else:
                        return (res+name+"\n"+self._sequence+"\n")

#########################################################################################
class ScytheGroupMap(ScytheBase):
        """Store a dictionary of Group IDs as keys and a list of gene models as values for one species."""
        def __init__(self, name = None, dct = None, gff = None, locfile = None, separator = None, asID = None):
                #The separator is useful in case ids in fasta and gff are different. asID can select the matching part of the ID
                super().__init__(name)
                self._mapping = {}
                self._mappingNames = {}
                self._parents = {}
                self._parentsNames = {}
                self._result = {}
                self._dct = None
                self._gff = None
                self._locfile = None
                if not dct: ##
                        raise ScytheError("Needs a dictionary pointing from a Group ID to a list of CDS.")
                else:
                        self._dct = dct
                if not gff and not locfile:
                        raise ScytheError("Need a gff file to add the other cds of the parent's locus")

                if locfile:
                        self._locfile = locfile
                        locfile = open(locfile, "r")
                        for l in locfile:
                                l = l.rstrip().split("\t")
                                self._mappingNames[l[0]] = l[1:]
                                for i in l[1:]:
                                        self._parentsNames[i] = l[0]
                elif gff:
                        self._gff = gff
                        for e in GFFParser().parse(gff):
                                if e.type == "mRNA":
                                        if e.attrib_dct["Parent"] not in self._mapping:
                                                self._mapping[e.attrib_dct["Parent"]] = [e.attrib_dct["ID"]]
                                                self._mappingNames[e.attrib_dct["Parent"]] = [e.attrib_dct["Name"]]
                                        else:
                                                self._mapping[e.attrib_dct["Parent"]].append(e.attrib_dct["ID"])
                                                self._mappingNames[e.attrib_dct["Parent"]].append(e.attrib_dct["Name"])
                                        #multiple parents? -> ToDo!
                                        self._parents[e.attrib_dct["ID"]] = e.attrib_dct["Parent"]
                                        self._parentsNames[e.attrib_dct["Name"]] = e.attrib_dct["Parent"]

                #often names in gff and fasta are different...try both 'Names' and 'ID' attribute
                tmpdct = self._dct
                res = {}
                if separator:
                        pass
                        #for k in self._dct:
                        #               tmplist = self._dct[k]
                        #               sublist = sum(tmplist,[]) #does that work at all?
                        #               print("sublist,tmplist", sublist,tmplist)
                        #               tmpdct[k] = [a.split(separator)[asID] for a in sublist]
                for k in tmpdct:
                                for sublist in tmpdct[k]:
                                        for e in sublist:
                                                #print(e, self._parents, self._parentsNames)
                                                if e in self._parentsNames:
                                                        res[k] = self._mappingNames[self._parentsNames[e]]
                                                        #print("result matches",e,self._parentsNames[e] )
                if not res:
                        for k in tmpdct:
                                for e in tmpdct[k]:
                                        for etmp in e:
                                                if etmp in self._parents:
                                                        #print(etmp, self._parents, self._parentsNames)
                                                        res[k] = self._mappingNames[self._parentsNames[etmp]]
                if not res:
                        print("WARNING: No orthogroups with species or non matching identifiers between group and loc ",locfile)
                        #raise ScytheError("You might want to check whether your identifiers match between gff and groups.")
                self._result = res
        @property
        def mapping(self):
                return self._mapping
        @property
        def mappingNames(self):
                return self._mappingNames
        @property
        def parents(self):
                return self._parents
        @property
        def parentsNames(self):
                return self._parentsNames
        @property
        def gff(self):
                return self._gff
        @property
        def dct(self):
                return self._dct
        @property
        def result(self):
                return self._result
        def free(self):
                self._dct = {}
                self._parents = {}
                self._parentsNames = {}
                self._mapping = {}
                self._mappingNames = {}

class ScytheGroup(ScytheBase):
        """Merge ScytheGroupMaps together"""
        def __init__(self, name = None, groupMaps = None):
                self._groups = {}
                self._names = []
                self._name = name
                self._groups = AutoVD()
                if groupMaps:
                        for g in groupMaps:
                                self._names.append(g.name)
                        for gm in groupMaps:
                                for k in gm.result:
                                        self._groups[k][gm.name] = gm.result[k]

        @property
        def groups(self):
                return self._groups

class AutoVD(dict):
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value


class ScytheSpec(ScytheBase):
    _knownFormats = ["loc"]
    def __init__(self, source = None, format = None, name= None, fasta = None, type = None):
        if format and format not in self._knownFormats:
            raise ScytheError(format+" not (yet) implemented; known formats: "+",".join(self._knownFormats))
        super().__init__(name)
        self._source =source
        self._format = format
        self._loci = {}
        self._cds = {}
        self._sequences = {}
        self._referenceForm = set()
        self._singleForm = set()
        self._defForm = set()
        self._fasta = fasta
        self._type = type

    def __iter__(self):
        for i in sorted(self._loci.keys()):
            yield (i, self._loci[i])

    @property
    def loci(self):
        return self._loci
    @loci.setter
    def loci(self, value):
        if not isinstance(value, dict):
            raise ScytheError("Should be dictionary")
        else:
            self._loci = dict(value)
    @property
    def cds(self):
        return self._cds
    @cds.setter
    def cds(self, value):
        if not isinstance(value, dict):
            raise ScytheError("Should be dictionary")
        else:
            self._cds = dict(value)

    @property
    def sequences(self):
        return self._sequences
    @sequences.setter
    def sequences(self, value):
        if not isinstance(value,ScytheSequence):
            raise ScytheError("Needs to be a ScytheSequence object")
        self._sequences = value
    @property
    def singleForm(self):
        return self._singleForm

    @property
    def referenceForm(self):
        return self._referenceForm

    @property
    def numCDS(self):
        return self._numCDS

    @property
    def fasta(self):
        return self._fasta
    @fasta.setter
    def fasta(self, value):
        """Source of .fasta file goes here """
        self._fasta = value
    @property
    def type(self):
        return self._type
    @type.setter
    def type(self, value):
        """ Type for .fasta goes here (pep/dna) """
        self._fasta = value
    def fillLociCDS(self):
        if self._format == "loc":
            f = open(self._source, "r")
            cnt = 0
            for i in f:
                i = i.rstrip()
                i = i.split("\t")
                loc, cds = i[0], i[1:]
                self._loci[loc] = cds
                if len(cds) == 1:
                    self._singleForm.add(cds[0])
                    self._referenceForm.add(cds[0])
                for c in cds:
                    cnt +=1
                    self._cds[c] = loc
        f.close()

    def fillSequences(self, fasta=None, type=None, sep = None, asID=None):
        """fill locus's sequence from fasta file 'fasta'.
        'type' should be 'pep' or 'dna'.
        If your identifiers have unnecessary additional
        suffixes, cut these off after 'sep'.
        """
        #print("# fill sequences -> ",self._name,self._cds)
        #todo asid
        if not fasta:
            fasta = self.fasta
        if not type:
            type = self.type
        for seq_record in FastaParser().read_fasta(fasta, sep):
            #print("sth is there...", seq_record)
            if seq_record[0] in self._cds:
                #print("works", seq_record)
                temp = ScytheSeq(name = seq_record[0], type = type, species = self._name, sequence = seq_record[1])
                self._sequences[seq_record[0]] = temp
    def fillDefForm(self, sep = "\t"):
        found = False
        """Add default gene models from 2nd column in .loc file.
        These forms will be saved in self._defForm.
        """
        lfile = open(self._source,"r")
        for l in lfile:
                l = l.rstrip()
                try:
                        l = l.split(sep)
                        if len(l) == 2:
                            self._sequences[l[1]].isSingle = True
                        l = l[1]
                except IndexError as e:
                        print (e, l)
                self._defForm.add(l)
                self._sequences[l].isReference=True
###########
class ScytheFrame(object):
    """Provide the general framework with directories etc and
    functionality to call external programs.
    """
    def __init__(self, path="."):
        self._path = path
        self._sr = self._path+"/ScytheResults/"
        self._srfa = self._sr+"/fasta/"
        self._srofa = self._sr+"/orthogroups_fasta/"
        self._infolog = self._sr+os.sep+"info.log"
        self._debuglog = self._sr+os.sep+"debug.log"
        self._errorlog = self._sr+os.sep+"error.log"
        ##########test for /dev/shm ######
        rdev = self.testRamDev()
        if rdev:
                self._fat = rdev
        #no /dev/shm, mac os
        else:
                self._fat = self._path+"/ScytheFastaTemp/"

    @property
    def path(self):
        return self._path
    @property
    def fat(self):
        return self._fat
    @property
    def srfa(self):
        return self._srfa

    def mkAllDirs(self):
        if not os.path.isdir(self._path):
            os.makedirs(self._path)
        if not os.path.isdir(self._sr):
            os.makedirs(self._sr)
        if not os.path.isdir(self._srfa):
            os.makedirs(self._srfa)
        if not os.path.isdir(self._srofa):
            os.makedirs(self._srofa)
        if not os.path.isdir(self._fat):
            os.makedirs(self._fat)
    def mkLogFiles(self):
        try:
                infohandle = open(self._infolog,'w')
                infohandle.write(datetime.datetime.now().strftime("%D_%T")+"\t started")
                infohandle.close()
                debughandle = open(self._debuglog,'w')
                debughandle.write(datetime.datetime.now().strftime("%D_%T")+"\t started")
                debughandle.close()
                errorhandle = open(self._errorlog,'w')
                errorhandle.close()
        except IOError as e:
                print(e)
                sys.exit(1)
    def writeLog(self,type="info", s=""):
        timestamp = time.time()
        file = ""
        s = str(timestamp)+"\t"+datetime.datetime.now().strftime("%D_%T")+"\t"+s
        if (type == "info"):
                file = self._infolog
        elif (type == "error"):
                file = self._errorlog
        elif (type == "debug"):
                file = self._debuglog
        else:
                print('type can only be "info","error" or "debug"')
        try:
                handle = open(file, 'a+')
                handle.write(s)
                handle.write("\n")
                handle.close()
        except IOError as e:
                print(e)
                print(self._sr,type)

    def testRamDev(self):
        rand = ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(2))
        ramdevpath = "/dev/shm/ScytheFastaTemp"+"_"+rand+"/"
        if not os.path.isdir("/dev/shm/"):
                return(False)
        if not os.path.isdir(ramdevpath):
                os.makedirs(ramdevpath)
        else:
                return(self.testRamDev())
        return(ramdevpath)



    def callNeedleAll(self, asequence, bsequence, outfile, gapOpen, gapExtend, aformat = "score", stype1="-sprotein1", stype2="-sprotein2" ,  stdout = False):
        """
        call 'needleall' from the emboss package
        stype1: either '-sprotein1' or '-snucleotide1'
        stype2: either '-sprotein2' or '-snucleotide2'
        """
        cmd = 'needleall -auto -asequence {0} -bsequence {1} -aformat {2} {3} {4} -outfile {5} -gapopen {6} -gapextend {7}'.format(
                                                                            asequence,
                                                                            bsequence,
                                                                            aformat,
                                                                            stype1,
                                                                            stype2,
                                                                            outfile,
                                                                            gapOpen,
                                                                            gapExtend
                                                                            )


        if stdout:
                cmd = 'needleall -auto -asequence {0} -bsequence {1} -aformat {2} {3} {4} -gapopen {5} -gapextend {6} -stdout'.format(
                                                                            asequence,
                                                                            bsequence,
                                                                            aformat,
                                                                            stype1,
                                                                            stype2,
                                                                            gapOpen,
                                                                            gapExtend
                                                                            )
        #self.writeLog("debug", cmd)
        ret = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, close_fds=True)
        #print("RET arrived")
        #retval = ret.wait()
        #ret.stdout.close()
        return(ret)

    #----------------------------------------------#


    def cleanUp(self):
        for f in os.listdir(self._fat):
                print(f)
                os.remove(os.sep.join([self._fat,f]))
