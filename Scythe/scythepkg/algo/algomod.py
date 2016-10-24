from itertools import chain
import logging
logger = logging.getLogger("algomod")
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = logging.FileHandler('scy_algomod.log')
fh.setLevel(logging.DEBUG)
# create console handler with a higher log level
ch = logging.StreamHandler()
ch.setLevel(logging.ERROR)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
ch.setFormatter(formatter)
# add the handlers to the logger
logger.addHandler(fh)
logger.addHandler(ch)

class AlgoHandler(object):
    def __init__(self):
        pass

    def getSingleGMS(self, sequenceDct):
        """Return list of species that only have one gene model"""
        res = [s.species for s in sequenceDct.values() if s.isSingle]
        return(res)

    def initCollProc(self, scoringDct, sequenceDct):
        """Initialize sets of species. """
        "speciesID -> [gmIDs]"
        species2id = {}

        #species:
        processed = set()
        unprocessed = set()

        #seq:
        coll = set()
        uncoll = set()

        for fkey in scoringDct:
            uncoll.add(fkey)
            unprocessed.add(sequenceDct[fkey].species)
            append_add(species2id, sequenceDct[fkey].species, fkey, unique = True)
            for skey in scoringDct[fkey]:
                uncoll.add(skey)
                unprocessed.add(sequenceDct[skey].species)
                append_add(species2id, sequenceDct[skey].species, skey, unique=True)
        return(processed, unprocessed, coll, uncoll, species2id)

    def getMaxSeed(self, scoringDct, sequenceDct, uncoll):
        """Find the best scoring pair. In case of a tie, prefer default models. """
        globMax = -1
        globMaxIds = None
        globMaxSps = None
        if scoringDct =={}:
            raise Warning("ScoringDct empty")
            return(None,None)
        defaultForms = [s.name for s in sequenceDct.values() if s.isReference]
        for u in uncoll: #uncollected species
            tmpspec = sequenceDct[u].species
            if scoringDct[u]:
                skl=list(scoringDct[u].keys())
                scorel=list(scoringDct[u].values())
                tmpMax = max(scorel)
                if tmpMax > globMax:
                    globMax = tmpMax
                    globMaxIds = (skl[scorel.index(tmpMax)], u)
                    tieList = [n for (n, e) in enumerate(scorel) if e == tmpMax]
                    if len(tieList) >1:
                        pass
                    #favoring defaults:
                    fav = [(skl[n],u)  for n in tieList if (skl[n] in defaultForms and u in defaultForms)]
                    fav1 = [(skl[n],u)  for n in tieList if (skl[n] in defaultForms or u in defaultForms)]
                    if fav:
                        globMaxSps = (sequenceDct[fav[0][0]].species, sequenceDct[fav[0][1]].species)
                        globMaxIds = (fav[0][0],fav[0][1])
                    elif fav1:
                        fav1tmp = fav1[0] #tuple
                        globMaxSps = (sequenceDct[fav1tmp[0]].species, sequenceDct[fav1tmp[1]].species)
                        globMaxIds =  (fav1tmp[0], fav1tmp[1])
                    else:
                        globMaxSps = (sequenceDct[globMaxIds[0]].species, sequenceDct[u].species)
                        globMaxIds = (globMaxIds[0], u)

                else:
                    pass
            else:
                sys.stderr.write(scoringDct, u)
                raise Exception("Sth wrong with ScoringDct")
        return(globMaxIds, globMaxSps)

    def sl_ref(self, scoringDct = {}, sequenceDct = {}, referenceAlgo = True):
        """Reference algorithm"""
        if not scoringDct:
            raise EmptyScoringDctException("sth wrong during sl_ref")
        if not sequenceDct:
            raise EmptyScoringDctException("sth wrong during sl_ref")
        singleGMSpec = self.getSingleGMS(sequenceDct)
        processed, unprocessed, coll, uncoll, species2id  = self.initCollProc(scoringDct, sequenceDct)
        if referenceAlgo :
        ### add single GM ###
            for sgms in set(singleGMSpec):
                assert len(species2id[sgms]) == 1
                processed.add(sgms)
                unprocessed.remove(sgms)
                coll.add(species2id[sgms][0])
                uncoll.remove(species2id[sgms][0])
        #####################

        while(unprocessed):
            logger.debug("unprocessed: {}".format(unprocessed))
            if not coll:
                logger.debug("coll: {} ".format(coll))
                max_seqid,max_specid = self.getMaxSeed(scoringDct = scoringDct, sequenceDct= sequenceDct, uncoll = uncoll)
                logger.debug("max_seqid: {}, max_specid: {} ".format(max_seqid, max_specid))
                for j,k in  zip(max_seqid,max_specid):
                    logger.debug("add/rm j: {}, coll/uncoll:{}, k: {}, processed: {}, unprocessed: {} ".format(
                        j,k,coll,uncoll,processed,unprocessed))
                    coll.add(j)
                    uncoll.remove(j)
                    processed.add(k)
                    unprocessed.remove(k)
                    logger.debug("add/rm j: {}, coll/uncoll:{}, {}, k: {}, processed: {}, unprocessed: {} ".format(
                        j,coll,uncoll,k,processed,unprocessed))
            for c in coll:
                logger.debug("c: {} ".format(c))
                cmax = -1
                cmaxid = None
                cmaxsp = None
                for up in unprocessed:#spec
                    logger.debug("up: {} ".format(up))
                    for uc in species2id[up]:#seq
                        if uc in scoringDct: # is in distance dict
                            try:
                                tmp = int(scoringDct[uc][c])
                                logger.debug("tmp score: {} ".format(tmp))
                            except TypeError as e:
                                tmp = -1
                            if (tmp > cmax or (tmp == cmax and sequenceDct[uc].isReference)) :
                                # tie resolved
                                cmaxid = uc
                                cmaxsp = up
                                cmax=tmp
                                logger.debug("cmaxid {}, cmaxsp {}, cmax {} ".format(cmaxid, cmaxsp, cmax))
                        elif c in scoringDct:
                            try:
                                tmp = int(scoringDct[c][uc])
                            except TypeError as e:
                                    tmp = -1
                            if (tmp >cmax or (tmp == cmax and sequenceDct[uc].isReference)):
                                # todo ??
                                cmax = int(avd[c][uc])
                                cmaxid = uc
                                cmaxsp = up
                                cmax=tmp
                                logger.debug("cmaxid {}, cmaxsp {}, cmax {} ".format(cmaxid, cmaxsp, cmax))
            logger.debug("uncoll {}, unprocessed {}, coll {} , processed {}".format(uncoll, unprocessed, uncoll, coll, cmaxsp))
            logger.debug("species2id {}".format(species2id))
            if unprocessed:
                uncoll.remove(cmaxid)
                unprocessed.remove(cmaxsp)
                coll.add(cmaxid)
                processed.add(cmaxsp)
        return(sequenceDct, coll, species2id)

    def sl_glob(self, scoringDct = {},sequenceDct = {}):
        print("debug", "sl_glob")
        if not scoringDct:
            raise EmptyScoringDctException("sth wrong during sl_glob")
        if not sequenceDct:
            raise EmptyScoringDctException("sth wrong during sl_glob")
        return(self.sl_ref(scoringDct = scoringDct, sequenceDct = sequenceDct, referenceAlgo=True))

    def mx_sum(self,scoringDct, sequenceDct):
        if not scoringDct:
            raise EmptyScoringDctException("sth wrong during mx_sum")
        if not sequenceDct:
            raise EmptyScoringDctException("sth wrong during mx_sum")

        processed, unprocessed, coll, uncoll, species2id  = self.initCollProc(scoringDct = scoringDct, sequenceDct = sequenceDct)
        pairwise,allkeys = getPairwiseAsTuples(scoringDct)
        actual = pairwise.copy()
        lot = actual.keys()
        keylengths = [len(key) for key in lot]
        while(max(keylengths) < len(unprocessed)): #proxy for num of spec
            newdct = adddyn(lot, pairwise, actual, allkeys, sequenceDct)
            actual = newdct
            lot = newdct.keys()
            keylengths = [len(key) for key in lot]
        tupkeys = []
        scores = []
        for k,v in actual.items():
            tupkeys.append(k)
            scores.append(v)
        tmax=max(scores)
        tieList = [(tupkeys[n],e,n) for (n, e) in enumerate(scores) if e == tmax]
#    #prefer defaults
        defaults = [u for u in uncoll if sequenceDct[u].isReference]
        tiedef = {}
        for tie in tieList:
            tiedef[tie[0]]=sum([1 for x in tie[0] if x in defaults])
        maxDefaults = max(tiedef.values())
        proc = [k for k, v in tiedef.items() if v == maxDefaults]
        first = proc[0]
        return(sequenceDct,set(first),species2id)

class EmptyScoringDctException(Exception):
    pass

class EmptySequenceDctException(Exception):
    pass


def append_add(dct, key, val, unique = False):
    """Append new (unique) value to the list of 'key'.
    If 'key' does not exist, create new list.
    """
    if key in dct:
        if unique:
            if val in dct[key]:
                pass
            else:
                dct[key].append(val)
        else:
            dct[key].append(val)
    else:
        dct[key] = [val]

def getPairwiseAsTuples(avd):
    firstdim = avd.keys()
    secdim = []
    pwdict={}
    listofkeys =[]
    for f in firstdim:
         if f not in listofkeys:
             listofkeys.append(f)
         for s in avd[f].keys():
            if s not in listofkeys:
                 listofkeys.append(s)
            pwdict[(f,s)]=  avd[f][s]
    return(pwdict.copy(), listofkeys)


def adddyn(listoftuples, pairwisedist, actualdict, allkeys, sequencesdct):
    tmdct = actualdict.copy()
    lenot = [len(t) for t  in listoftuples]
    maxlen = max(lenot)
    listoftuples = [t for t in listoftuples if len(t) == maxlen]
    for t in listoftuples:#start with pairwise dist (i,j)
         #print("TUPLE ",t)
         specdone=[sequencesdct[e].species for e in t]
         for k in allkeys: #should be single keys
            #print("ALLKEYS k",k)
            if sequencesdct[k].species not in specdone:
                #if k not in t: #and species not in etc
                newtup=tuple(chain.from_iterable([t,[k]]))
                if newtup not in tmdct:
                    #print("NEWTUP",newtup)
                    specdone.append(sequencesdct[k].species)
                    addscore=0
                    for l in t:
                        try:
                            addscore = addscore+pairwisedist[(l,k)]
                        except KeyError as ke:
                            try:
                                addscore = addscore+pairwisedist[(k,l)]#both di i,j an j,i
                            except KeyError:
                                sys.stderr.write("KEYERROR",ke)
                    tmdct[newtup] = tmdct[t]+addscore
    return(tmdct.copy())

