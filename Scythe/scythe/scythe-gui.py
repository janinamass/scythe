#!/usr/bin/env python3
import os
import sys
import re

import tkinter as tk
import configparser
import multiprocessing

from tkinter import filedialog
from tkinter import messagebox
from tkinter import ttk
from tkinter import OptionMenu
from tkinter import Scale
from tkinter import Listbox

import scythe as scythe

#from gui.dialogs import ScytheConvertDialogLoc
#from gui.dialogs import ScytheConvertDialogGrp

import helpers.mergeSubsets as mergeSubsets
import helpers.ensembl_ortho_mysql as ensembl_ortho_mysql
import helpers.ensembl2grp as ensembl2grp
import helpers.ensembl as ensembl

import multiprocessing
import queue

wd = os.path.join(os.path.dirname(__file__))

root=tk.Tk()
root.title("Scythe GUI 0.1.0")
try:
    root.iconbitmap("@"+wd+os.sep+"gui"+os.sep+"scy.xbm")
except tk.TclError as e:
    print(e)
    sys.stderr.write(str(e))
    root.iconbitmap(None)

global SCYTHE_PROCESS
SCYTHE_PROCESS = None

global CURRENTCONFIG
CURRENTCONFIG = configparser.ConfigParser()
BACKUPCONFIG = configparser.ConfigParser()
######################################
###       config files            ####
############ labels ##################
CF_MODE = "Mode"
CF_MODE_use_ensembl = "use_ensembl"
CF_MODE_use_local_files = "use_local_files"
CF_PATHS = "Paths"
CF_PATHS_fasta_directory = "fasta_directory"
CF_PATHS_loc_directory = "loc_directory"
CF_PATHS_grp_file = "grp_file"
CF_PATHS_output_directory = "output_directory"
CF_CLEANUP = "Cleanup"
CF_CLEANUP_clean_up_directories = "clean_up_directories"
CF_RUN="Run_options"
CF_RUN_num_CPU ="num_CPU"
CF_PENALTIES = "Penalties"
CF_PENALTIES_gap_open_cost = "gap_open_cost"
CF_PENALTIES_gap_extend_cost="gap_extend_cost"
CF_PENALTIES_substitution_matrix="substitution_matrix"
CF_ALGORITHM = "Algorithm"
CF_ALGORITHM_use_global_max ="use_sl_glob"
CF_ALGORITHM_use_default="use_sl_ref"
CF_ALGORITHM_use_global_sum="use_mx_sum"
CF_FASTAHEADER="Fasta_header"
CF_FASTAHEADER_delimiter = "fasta_header_delimiter"
CF_FASTAHEADER_part = "fasta_header_part"
############################

##################options###########################
OPTIONS = {}
#dropdown menus
yn =["yes","no"]
for o in [CF_ALGORITHM_use_global_max,CF_ALGORITHM_use_default,
          CF_ALGORITHM_use_global_sum, CF_CLEANUP_clean_up_directories,
          CF_MODE_use_local_files,CF_MODE_use_ensembl]:
    OPTIONS[o]=yn

MAXCPU=multiprocessing.cpu_count

SECTIONS = [CF_MODE,CF_PATHS, CF_CLEANUP, CF_RUN, CF_PENALTIES,CF_ALGORITHM, CF_FASTAHEADER]

#setup dummy config
MAXCONFIG = configparser.ConfigParser()
for i in SECTIONS:
    MAXCONFIG.add_section(i)
for i in [CF_MODE_use_ensembl,CF_MODE_use_local_files]:
    MAXCONFIG.set(CF_MODE,i,"unset")
for i in [CF_PATHS_fasta_directory,CF_PATHS_loc_directory,CF_PATHS_grp_file,CF_PATHS_output_directory ]:
    MAXCONFIG.set(CF_PATHS,i,"unset")
for i in [CF_CLEANUP_clean_up_directories]:
    MAXCONFIG.set(CF_CLEANUP,i,"yes")
for i in [CF_RUN_num_CPU]:
    MAXCONFIG.set(CF_RUN,i,"1")

for i in [CF_PENALTIES_gap_open_cost,CF_PENALTIES_gap_extend_cost,CF_PENALTIES_substitution_matrix]:
    if i == CF_PENALTIES_gap_open_cost:
        MAXCONFIG.set(CF_PENALTIES,i,"10")
    if  i == CF_PENALTIES_gap_extend_cost:
         MAXCONFIG.set(CF_PENALTIES,i,"0.5")
    if i == CF_PENALTIES_substitution_matrix:
        MAXCONFIG.set(CF_PENALTIES,i,"EBLOSUM62")

for i in [CF_ALGORITHM_use_global_max,CF_ALGORITHM_use_default,CF_ALGORITHM_use_global_sum   ]:
    MAXCONFIG.set(CF_ALGORITHM,i,"unset")

for i in [CF_FASTAHEADER_delimiter, CF_FASTAHEADER_part]:
    if i == CF_FASTAHEADER_delimiter:
        MAXCONFIG.set(CF_FASTAHEADER,i,'" "')
    if i == CF_FASTAHEADER_part:
        MAXCONFIG.set(CF_FASTAHEADER,i,"0")
####

##########################
# /config files
##########################
class ConfigHandler():
    def __init__(self):
        self._currentconfig = configparser.ConfigParser()
    @property
    def currentconfig(self):
        return self._currentconfig
    @currentconfig.setter
    def currentconfig(self,config):
        self._currentconfig=config

    def reset(self):
        self.initConfCurrent()
        print("full reset")
#q'n'd
    def initConfCurrent(self):
        global CURRENTCONFIG
        global MAXCONFIG
        for i in MAXCONFIG.sections():
            try:
                CURRENTCONFIG.add_section(i)
            except configparser.DuplicateSectionError as e:
                 pass
            for j in MAXCONFIG.options(i):
                #print(i,j)
                CURRENTCONFIG.set(i,j,MAXCONFIG.get(i,j))

    def backupConf(self):
        global CURRENTCONFIG
        global BACKUPCONFIG
        for i in CURRENTCONFIG.sections():
            try:
                BACKUPCONFIG.add_section(i)
            except configparser.DuplicateSectionError as e:
                 pass
            for j in CURRENTCONFIG.options(i):
                BACKUPCONFIG.set(i,j,CURRENTCONFIG.get(i,j) )

    def backupConfTo(self, newconf):
        global CURRENTCONFIG
        for i in CURRENTCONFIG.sections():
            try:
                newconf.add_section(i)
            except configparser.DuplicateSectionError as e:
                 pass
            for j in CURRENTCONFIG.options(i):
                newconf.set(i,j,CURRENTCONFIG.get(i,j) )

    def setCurrentConf(self, newconf):
        global CURRENTCONFIG
        for i in newconf.sections():
            try:
                CURRENTCONFIG.add_section(i)
            except configparser.DuplicateSectionError as e:
                 pass
            for j in newconf.options(i):
                CURRENTCONFIG.set(i,j,newconf.get(i,j) )

    def restoreConf(self):
         global BACKUPCONFIG
         global CURRENTCONFIG
         for i in BACKUPCONFIG.sections():
             try:
                 CURRENTCONFIG.add_section(i)
             except configparser.DuplicateSectionError as e:
                 pass
             for j in  BACKUPCONFIG.options(i):
                 CURRENTCONFIG.set(i,j,BACKUPCONFIG.get(i,j))

class ScytheConfigEditor():
    def __init__(self):

        global CURRENTCONFIG
        global MAXCONFIG
        global CF_MODE
        self.confighandler = ConfigHandler()
        self.confighandler.backupConf()
        tmpconfig= configparser.ConfigParser()
        self.confighandler.backupConfTo(tmpconfig)
        top = tk.Toplevel()
        top.title("Set configuration")
        nb = ttk.Notebook(top)
        b_config_ok = tk.Button(top, text="OK", command=top.destroy)
        b_config_ok.bind('<ButtonRelease-1>',self.onSetConfigOK)
        b_config_apply = tk.Button(top, text="Apply", command=self.onSetConfigApply)
        b_config_cancel = tk.Button(top, text="Cancel", command=top.destroy)
        b_config_cancel.bind('<ButtonRelease-1>',self.onSetConfigCancel())


        fr_paths = tk.Frame(nb,width=200, height=100)
        fr_penalties = tk.Frame(nb,width=200, height=100)
        fr_mode = ttk.Frame(nb,width=200, height=100)
        fr_cleanup = ttk.Frame(nb,width=200, height=100)
        fr_run = ttk.Frame(nb,width=200, height=100)
        fr_algorithm = ttk.Frame(nb,width=200, height=100)
        fr_fastaheader = ttk.Frame(nb,width=200, height=100)

        #######labels########################
        self.txt_sec=[]
        self.txt_subsec={}
        for section in MAXCONFIG.sections():
            #print( "["+section +"]\n")
            self.txt_sec.append(section)
            for opt in MAXCONFIG.options(section):
                try:
                    self.txt_subsec[section].append(opt)
                except KeyError as e:
                    self.txt_subsec[section]=[opt]
        lab_sec=[]
        lab_subsec={}
        dd_subsec={}
        self.var_subsec={}
        for t in self.txt_sec:
            lab_sec.append(tk.Label(fr_paths,text = t))
        for t in self.txt_subsec:
            #print(t,self.txt_subsec[t])
            for u in self.txt_subsec[t]:
                if t == CF_MODE:
                    fr = fr_mode
                elif t == CF_PATHS:
                    fr = fr_paths
                elif t == CF_CLEANUP:
                    fr = fr_cleanup
                elif t == CF_RUN:
                    fr = fr_run
                elif t == CF_PENALTIES:
                    fr = fr_penalties
                elif t == CF_ALGORITHM:
                    fr = fr_algorithm

                elif t == CF_FASTAHEADER:
                    fr = fr_fastaheader
                    #print("fastaheader_fr")
                ################################
                else:
                    sys.stderr.write("No such section:",t)
                try:
                    lab_subsec[t].append(tk.Label(fr,text = u))
                    self.var_subsec[t].append(tk.StringVar(fr))
                    if u in OPTIONS:
                        dd_subsec[t].append(OptionMenu(fr,self.var_subsec[t][-1],*OPTIONS[u]))
                    else:
                        dd_subsec[t].append("")
                except KeyError as e:
                    try:
                        lab_subsec[t]=[tk.Label(fr,text = u)]
                        self.var_subsec[t]=[tk.StringVar(fr)]
                        if u in OPTIONS:
                            dd_subsec[t] = [OptionMenu(fr,self.var_subsec[t][-1],*OPTIONS[u])]
                        else:
                            dd_subsec[t] = [""]
                    except KeyError as e:
                        sys.stderr.write(str(e))
                        dd_subsec[t].append("")

        for t in lab_subsec:
            r=0
            c=0
            for i in  lab_subsec[t]:
                #print(i.cget("text"))
                i.grid(row=r,column=c, sticky=tk.E)
                r+=1
                #print(r,i.cget("text"))
        for t in dd_subsec:
            c=1
            r=0
            for i in dd_subsec[t]:
                #print(i)
                if i is not "":
                    i.grid(row=r,column=c,sticky=tk.N)
                r+=1
                #print(r)
        ######################################
        self.st_submat = tk.StringVar()
        self.st_fasta_header_delimiter = tk.StringVar()
        self.st_fasta_header_part = tk.StringVar()
        #cpu_count starts at 0 for one cpu
        self.sc_config_numthreads = Scale(fr_run, from_=1, to=multiprocessing.cpu_count(), orient=tk.HORIZONTAL)
        self.sc_config_numthreads.grid(row=0, column=1, sticky=tk.E)
        en_config_gapopen=tk.Entry(fr_penalties, textvariable=self.var_subsec[CF_PENALTIES][0])
        en_config_gapextend=tk.Entry(fr_penalties,textvariable=self.var_subsec[CF_PENALTIES][1] )
        self.en_config_fasta_header_delimiter= tk.Entry(fr_fastaheader,textvariable=self.st_fasta_header_delimiter,width=6 )
        self.en_config_fasta_header_part= tk.Entry(fr_fastaheader,textvariable=self.st_fasta_header_part ,width=6 )


        self.om_config_submat=tk.OptionMenu(fr_penalties, self.st_submat, *["EBLOSUM62","EDNAFULL"])
        self.om_config_submat.grid(row=2,column=1 )
        en_config_gapopen.grid(row=0, column=1)
        en_config_gapextend.grid(row=1, column=1)

        self.en_config_fasta_header_delimiter.grid(row=0, column=1)
        self.en_config_fasta_header_part.grid(row=1,column=1)
        nb.add(fr_penalties, text=CF_PENALTIES)
        nb.add(fr_cleanup, text=CF_CLEANUP)
        nb.add(fr_run, text=CF_RUN)
        nb.add(fr_algorithm, text=CF_ALGORITHM)
        nb.add(fr_fastaheader, text=CF_FASTAHEADER)

        nb.grid()
        b_config_cancel.grid(row=1, column=0, sticky=tk.E,padx=115)
        b_config_apply.grid(row=1, column=0, sticky=tk.E,padx=50)
        b_config_ok.grid(row=1, column=0, sticky=tk.E)
        self.setFieldsFromConfig()

    def onSetConfigApply(self):
        self.setConfigFromFields()

    def onSetConfigOK(self,event):
        self.setConfigFromFields()

    def onSetConfigCancel(self):
        self.confighandler.restoreConf()
        #print("RESTORED-->CURRENTCONF set")
        #print("Config CANCEL")

    def setConfigFromFields(self):
        tempconf = configparser.ConfigParser()
        self.confighandler.backupConfTo(tempconf)
        #get all values from fields
        #penalties
        tempconf.set(CF_PENALTIES,CF_PENALTIES_gap_open_cost,self.var_subsec[CF_PENALTIES][0].get() )
        tempconf.set(CF_PENALTIES, CF_PENALTIES_gap_extend_cost,self.var_subsec[CF_PENALTIES][1].get())
        tempconf.set(CF_PENALTIES, CF_PENALTIES_substitution_matrix,self.st_submat.get())
        tempconf.set(CF_ALGORITHM, CF_ALGORITHM_use_global_max,self.var_subsec[CF_ALGORITHM][0].get())
        tempconf.set(CF_ALGORITHM, CF_ALGORITHM_use_default,self.var_subsec[CF_ALGORITHM ][1].get())
        tempconf.set(CF_ALGORITHM, CF_ALGORITHM_use_global_sum,self.var_subsec[CF_ALGORITHM][2].get())
        tempconf.set(CF_RUN, CF_RUN_num_CPU,str(self.sc_config_numthreads.get()))
        #CLEANUP
        tempconf.set(CF_CLEANUP, CF_CLEANUP_clean_up_directories, self.var_subsec[CF_CLEANUP][0].get())
        #Fasta header
        tempconf.set(CF_FASTAHEADER, CF_FASTAHEADER_delimiter, self.var_subsec[CF_FASTAHEADER][0].get())
        tempconf.set(CF_FASTAHEADER, CF_FASTAHEADER_part, self.var_subsec[CF_FASTAHEADER][1].get())
        tempconf.set(CF_FASTAHEADER, CF_FASTAHEADER_part,self.st_fasta_header_part.get())
        tempconf.set(CF_FASTAHEADER, CF_FASTAHEADER_delimiter,self.st_fasta_header_delimiter.get())

        self.confighandler.setCurrentConf(tempconf)

    def setFieldsFromConfig(self):
        #penalties
        self.var_subsec[CF_PENALTIES][0].set(CURRENTCONFIG.get(CF_PENALTIES,self.txt_subsec[CF_PENALTIES][0]))
        self.var_subsec[CF_PENALTIES][1].set(CURRENTCONFIG.get(CF_PENALTIES,self.txt_subsec[CF_PENALTIES][1]))
        self.st_submat.set(CURRENTCONFIG.get(CF_PENALTIES, CF_PENALTIES_substitution_matrix))
        #output
        #cleanup
        self.var_subsec[CF_CLEANUP][0].set(CURRENTCONFIG.get(CF_CLEANUP,self.txt_subsec[CF_CLEANUP][0]))
        #run
        #slider
        #algo
        self.var_subsec[CF_ALGORITHM][0].set(CURRENTCONFIG.get(CF_ALGORITHM,self.txt_subsec[CF_ALGORITHM][0]))
        self.var_subsec[CF_ALGORITHM][1].set(CURRENTCONFIG.get(CF_ALGORITHM,self.txt_subsec[CF_ALGORITHM][1]))
        self.var_subsec[CF_ALGORITHM][2].set(CURRENTCONFIG.get(CF_ALGORITHM,self.txt_subsec[CF_ALGORITHM][2]))

        self.var_subsec[CF_FASTAHEADER][0].set(CURRENTCONFIG.get(CF_FASTAHEADER,self.txt_subsec[CF_FASTAHEADER][0]))
        self.var_subsec[CF_FASTAHEADER][1].set(CURRENTCONFIG.get(CF_FASTAHEADER,self.txt_subsec[CF_FASTAHEADER][1]))
        self.st_fasta_header_part.set(CURRENTCONFIG.get(CF_FASTAHEADER, CF_FASTAHEADER_part))
        self.st_fasta_header_delimiter.set(CURRENTCONFIG.get(CF_FASTAHEADER, CF_FASTAHEADER_delimiter))

def logged(f):
    def wrapped(*args, **kargs):
        print ("%s called..." % f.__name__)
        try:
            return f(*args, **kargs)
        finally:
            print ("..Done.")
    return wrapped

class Infobox():
    @logged
    def todo(self):
        message="Todo..."
        messagebox.showinfo(title="Todo...", message = message )
    @logged
    def about(self):
        message="Scythe GUI v0.1.0 (June 2014)\nJ. Mass\nThis is under construction."
        messagebox.showinfo(title="About Scythe", message = message )

    def bepatient(self):
        message="Scythe is running\n This may take some time.\n"
        messagebox.showinfo(title="Running", message = message )

    def saveConfig(self):
        formats = [('Scythe configuration','*.scy')]
        tmp= tk.filedialog.asksaveasfilename(parent=self.parent,filetypes=formats ,title="Save configuration as...")

    def showConfig(self):
        global CURRENTCONFIG
        tmp = ConfigHandler()
        tmp.currentconfig=CURRENTCONFIG
        message = ""
        for section in tmp.currentconfig.sections():
            message += "["+section +"]\n"
            for option in tmp.currentconfig.options(section):
                message += " "+ option+ "="+ tmp.currentconfig.get(section, option)+"\n"
        top = tk.Toplevel(root)
        top.title("Configuration")
        txt = tk.Text(top)
        scrollv = tk.Scrollbar(top, command=txt.yview)
        txt.insert(tk.INSERT,message)
        txt.configure(yscrollcommand=scrollv.set, state=tk.DISABLED, background="black", foreground="green" )
        txt.grid(row=0, column=0)
        scrollv.grid(row=0, column=1)


class ScytheMenu(tk.Frame):

    def __init__(self, parent, arg = None):
        tk.Frame.__init__(self, parent)
        self.parent = parent
        self.initGUI()
        self.confighandler = ConfigHandler()
        self.confighandler.initConfCurrent()
        self.scythewizard= ScytheWizard(self.parent)
        self.configEditor = None
        if arg:
            print(arg)
            self.loadConfigArg(arg)

    def initGUI(self):
        menubar = tk.Menu(self.parent)
        self.parent.config(menu=menubar)

        fileMenu = tk.Menu(menubar)
        fileMenu.add_command(label="New run...", command=self.onNewRun)
        #fileMenu.add_command(label="Convert files...", command=self.onConvertFiles)
        fileMenu.add_command(label="Load configuration...", command=self.onLoadConfig)
        fileMenu.add_command(label="Save configuration...", command=self.onSaveConfig)
        fileMenu.add_command(label="Exit", command=self.onExit)

        optionsMenu = tk.Menu(menubar)
        optionsMenu.add_command(label="Show configuration...", command=self.onShowOptions)
        optionsMenu.add_command(label="Set configuration...", command=self.onSetOptions)

        helpMenu = tk.Menu(menubar)
        helpMenu.add_command(label="About...", command=self.onAbout)

        menubar.add_cascade(label="File", menu=fileMenu)
        menubar.add_cascade(label="Options", menu=optionsMenu)
        menubar.add_cascade(label="Help", menu=helpMenu)

    def onExit(self):
        self.quit()
    def onNewRun(self):
        self.scythewizard=ScytheWizard(self.parent)
        self.confighandler.reset()

    def loadConfigArg(self,arg):
        cfg = self.confighandler.currentconfig.read(arg)
        global CURRENTCONFIG
        CURRENTCONFIG=self.confighandler.currentconfig
        self.scythewizard.st_fastaDir.set(self.confighandler.currentconfig.get(CF_PATHS,CF_PATHS_fasta_directory) )
        self.scythewizard.st_locDir.set(self.confighandler.currentconfig.get(CF_PATHS,CF_PATHS_loc_directory) )
        self.scythewizard.st_grpFile.set(self.confighandler.currentconfig.get(CF_PATHS,CF_PATHS_grp_file) )
        self.scythewizard.st_outDir.set(self.confighandler.currentconfig.get(CF_PATHS,CF_PATHS_output_directory) )
        if self.confighandler.currentconfig.get(CF_MODE,CF_MODE_use_ensembl)=='yes':
            self.scythewizard.cb_use_ensembl.select()
        elif self.confighandler.currentconfig.get(CF_MODE,CF_MODE_use_local_files)=='yes':
            self.scythewizard.cb_use_local.select()
        self.scythewizard.cb_use_local.configure(state=tk.NORMAL)
        self.scythewizard.cb_use_ensembl.configure(state=tk.NORMAL)
        self.onSetOptions()

    def onLoadConfig(self):
        formats = [('Scythe configuration','*.scy')]
        tmp = tk.filedialog.askopenfilename(parent=self.parent,filetypes=[('Scythe configuration','*.scy')],title="Load configuration...")
        cfg = self.confighandler.currentconfig.read(tmp)
        global CURRENTCONFIG
        CURRENTCONFIG=self.confighandler.currentconfig
        if self.confighandler.currentconfig.get(CF_MODE,CF_MODE_use_local_files)=='yes':
            self.scythewizard.st_fastaDir.set(self.confighandler.currentconfig.get(CF_PATHS,CF_PATHS_fasta_directory) )
            self.scythewizard.st_locDir.set(self.confighandler.currentconfig.get(CF_PATHS,CF_PATHS_loc_directory) )
            self.scythewizard.st_grpFile.set(self.confighandler.currentconfig.get(CF_PATHS,CF_PATHS_grp_file) )
            self.scythewizard.st_outDir.set(self.confighandler.currentconfig.get(CF_PATHS,CF_PATHS_output_directory) )
            self.scythewizard.cb_use_local.select()
            self.scythewizard.cb_use_local.configure(state=tk.NORMAL)
            self.scythewizard.cb_use_ensembl.configure(state=tk.DISABLED)
            self.scythewizard.ent_fastaDir.configure(state=tk.NORMAL)
            self.scythewizard.ent_locDir.configure(state=tk.NORMAL)
            self.scythewizard.ent_grpFile.configure(state=tk.NORMAL)
            self.scythewizard.ent_outDir.configure(state=tk.NORMAL)
        elif self.confighandler.currentconfig.get(CF_MODE,CF_MODE_use_ensembl)=='yes':
            self.scythewizard.st_outDir.set(self.confighandler.currentconfig.get(CF_PATHS,CF_PATHS_output_directory) )
            self.scythewizard.cb_use_ensembl.select()
            self.scythewizard.cb_use_ensembl.configure(state=tk.NORMAL)
            self.scythewizard.cb_use_local.configure(state=tk.DISABLED)
            self.scythewizard.ent_fastaDir.configure(state=tk.DISABLED)
            self.scythewizard.ent_locDir.configure(state=tk.DISABLED)
            self.scythewizard.ent_grpFile.configure(state=tk.DISABLED)
            self.scythewizard.ent_outDir.configure(state=tk.NORMAL)

        return tmp

    def onSaveConfig(self):
        self.scythewizard.setConfigFromFields()
        formats = [('Scythe configuration','*.scy')]
        tmp= tk.filedialog.asksaveasfilename(parent=self.parent,filetypes=formats ,title="Save configuration as...")
        global CURRENTCONFIG
        try:
            out = open(tmp,"w")
            CURRENTCONFIG.write(out)
            out.close()
        except Error as e:
            print (e)
        return tmp

    def onAbout(self):
        Infobox().about()

    def onShowOptions(self):
        self.scythewizard.setConfigFromFields()
        Infobox().showConfig()

    def onSetOptions(self):
        self.configEditor = ScytheConfigEditor()

class EnsemblSelector(tk.Listbox):
        """Show available species on ENSEMBL"""
        lb = None
        top= None

        speclist = None
        rellist = None
        itemlist = None
        outdir = None
        def __init__(self, outdir):
            print("")
            speclist = []
            rellist = []
            itemlist = []
#todo ensembl specInfo
            data = ensembl.specInfo()
            print(data)
            self.data = data
            self.outdir = outdir
            top = tk.Toplevel(root)
            #adjust width and height
            top.title("Select Species from Ensembl")
            self.top=top

            lb = Listbox(self.top,selectmode='multiple',exportselection=0 ,width=40, height=30,)
            self.lb= lb
            #print(data)
            for d in data["species"]:
                print(d["name"])
                if not d["name"].startswith("Ancestral"):
                    speclist.append(d["name"])
                    itemlist.append(d["name"]+'_core_'+str(d["release"]))
                    rellist.append(d["release"])
            self.itemlist=itemlist
            self.speclist=speclist
            self.rellist = rellist
            self.b_ensOK = tk.Button(self.top, text="OK", command=self.onEnsOK)
            self.b_ensQuit = tk.Button(self.top, text="Cancel", command=self.onEnsQuit)
            self.b_ensOK.grid(row=1, column=0,sticky="E", padx=60)
            self.b_ensQuit.grid(row=1, column=0,sticky="E", padx=0)
            self.prepRun(itemlist)

        def fileExists(self,filename):
            try:
                tmp = open(filename,'r')
                tmp.close()
            except IOError as e:
                return(False)
            return(True)

        def onEnsOK(self):
            res = tk.messagebox.showinfo("Be patient.", "Download data. This could take some time.")
            if not res:
                self.top.destroy()
            specs,rel = self.readListBox()
            print("wait...", self.outdir, specs, rel)
            self.b_ensOK.configure(state=tk.DISABLED)
            self.b_ensQuit.configure(state=tk.DISABLED)
            fapath = self.outdir+os.sep+"fa"
            tmp = [s for s in specs if not self.fileExists(fapath+os.sep+s+".fa")]
            present = [s for s in specs if self.fileExists(fapath+os.sep+s+".fa")]
            for u in present:
                print("already there: "+fapath+os.sep+u+".fa")
            if (len(tmp)>0):
                ensembl.getSequencesFromFTP(self.outdir, rel, tmp)
            locpath = self.outdir+os.sep+"loc"
            print("fasta download done: {}".format(self.outdir))
            for i in specs:
                print(i)
                if not self.fileExists(locpath+os.sep+i+".loc"):
                    try:
                        ensembl.prepareLocFromFasta(fapath+os.sep+i+".fa",locpath+os.sep,i  )
                    except IOError as e:
                        sys.stderr.write(str(e))
                        print("Warning: No such fasta: ",fapath+os.sep+i+".fa")
                else:
                    print("already there: "+locpath+os.sep+u+".loc")
                    ###TODO deal with different releases: Throw warning, has to be done manually

            grpstring =""
            for i in specs:
                grpstring+=i[0:2]
            grpfile = self.outdir+os.sep+grpstring+"_tmp.grp"
            if self.fileExists(grpfile):
                 print("alredy there:", grpfile)
            else:
                listoftsv=ensembl_ortho_mysql.fetchOrthoFromMySQL(specieslist = specs, release=rel[0])
                ensembl2grp.readTsvFiles(listoftsv=listoftsv, outfile=grpfile)
                #####update
            outNoSubsets=self.outdir+os.sep+grpstring+".full.grp"
            out = self.outdir+os.sep+grpstring+".shared_by_all.grp"
            numspec= mergeSubsets.filterGroups(grpfile,None, None, False)
            mergeSubsets.mergeSubsets(grpfile,outNoSubsets, True)
            mergeSubsets.filterGroups(outNoSubsets, out, numspec, False)
            #destroy window
            self.top.destroy()

            CURRENTCONFIG.set(CF_PATHS,CF_PATHS_fasta_directory, fapath+os.sep)
            CURRENTCONFIG.set(CF_PATHS,CF_PATHS_loc_directory, locpath+os.sep)
            CURRENTCONFIG.set(CF_PATHS,CF_PATHS_grp_file, out)
            ScytheWizard(root).prepRun(reloadFields=False)

        def onEnsQuit(self):
            self.top.destroy()

        def prepRun(self, itemlist):
            #print("EnsemblSelectorPrepRun")
            for item in itemlist:
                self.lb.insert(tk.END, item)
            self.lb.grid(row=0, column=0)

        def readListBox(self):
            items = self.lb.curselection()
            #print(items)
            print(self.speclist)
            selecteditemsspec = [self.speclist[int(item)] for item in items]
            selecteditemsrel = [self.rellist[int(item)] for item in items]

            print(selecteditemsspec, selecteditemsrel)
            return(selecteditemsspec, selecteditemsrel)

class ScytheWizard(tk.Tk):
    def __init__(self, parent):
        self.parent = parent
        self.initWizard()
        self.confighandler = ConfigHandler()
    def quit(self):
        print("terminate")
        if SCYTHE_PROCESS:
            SCYTHE_PROCESS.terminate()
        root.destroy()
    def setConfigFromFields(self):
        tempconf = configparser.ConfigParser()
        self.confighandler.backupConfTo(tempconf)
        tempconf.set(CF_PATHS, CF_PATHS_output_directory,self.ent_outDir.get())
        tempconf.set(CF_PATHS, CF_PATHS_fasta_directory,self.ent_fastaDir.get())
        tempconf.set(CF_PATHS, CF_PATHS_loc_directory,self.ent_locDir.get())
        tempconf.set(CF_PATHS, CF_PATHS_grp_file,self.ent_grpFile.get())

        self.confighandler.setCurrentConf(tempconf)
        print(CURRENTCONFIG)

    def prepRun(self, reloadFields=True):
        global SCYTHE_PROCESS
        scythe.VERBOSE=False

        #update config one more time
        if reloadFields:
            self.setConfigFromFields()
            print("Read entry fields one more time", CURRENTCONFIG.get(CF_PATHS, CF_PATHS_output_directory))

        ###############################################
        #check whether ensembl or local is checked
        useEnsembl= CURRENTCONFIG.get(CF_MODE, CF_MODE_use_ensembl)
        useLocal = CURRENTCONFIG.get(CF_MODE, CF_MODE_use_local_files)
        #outdir to
        outdir = CURRENTCONFIG.get(CF_PATHS,CF_PATHS_output_directory)
        #catch unset outdir
        if not outdir:
            outdir = os.getcwd()
        cleanUp="yes"
        scythe.GLOBMAX = False
        scythe.GLOBSUM = False

        if useEnsembl == "yes" and reloadFields: #has just come back from EnsemblSelector
            ens = EnsemblSelector(outdir)
            print("Will download from ENSEMBL")
        else:

            if  useLocal == "yes":
                fastaDir =  CURRENTCONFIG.get(CF_PATHS,CF_PATHS_fasta_directory)
                locDir =  CURRENTCONFIG.get(CF_PATHS,CF_PATHS_loc_directory)
                grpFile = CURRENTCONFIG.get(CF_PATHS,CF_PATHS_grp_file)

                print("Will use local files")
                print(fastaDir,locDir,grpFile,outdir)
                ################
            if CURRENTCONFIG.get(CF_ALGORITHM,CF_ALGORITHM_use_global_max)!="yes":
                scythe.GLOBMAX = False
            else:
                scythe.GLOBMAX = True
            if CURRENTCONFIG.get(CF_ALGORITHM,CF_ALGORITHM_use_global_sum)!="yes":
                scythe.GLOBSUM = False
            else:
                scythe.GLOBSUM = True

            if CURRENTCONFIG.get(CF_CLEANUP,CF_CLEANUP_clean_up_directories) !="yes":
                cleanUp = False
            else:
                cleanUp = True


        groups= CURRENTCONFIG.get(CF_PATHS,CF_PATHS_grp_file)
        namesList = None
        faDir = CURRENTCONFIG.get(CF_PATHS,CF_PATHS_fasta_directory)
        if not faDir.endswith(os.sep):
            faDir+=os.sep
        inDir = faDir
        outDir = CURRENTCONFIG.get(CF_PATHS,CF_PATHS_output_directory)
        if not outDir:
            outDir = os.getcwd()
        if not outDir.endswith(os.sep):
            outDir+=os.sep
        locDir = CURRENTCONFIG.get(CF_PATHS,CF_PATHS_loc_directory)
        if not locDir.endswith(os.sep):
            locDir+=os.sep
        fastaList = os.listdir(faDir)
        delim = CURRENTCONFIG.get(CF_FASTAHEADER,CF_FASTAHEADER_delimiter).strip('"')
        try:
            asID = int(CURRENTCONFIG.get(CF_FASTAHEADER,CF_FASTAHEADER_part))
        except ValueError as e:
            print(e)
            asID=None
        gapOpen = CURRENTCONFIG.get(CF_PENALTIES,CF_PENALTIES_gap_open_cost)
        gapExtend = CURRENTCONFIG.get(CF_PENALTIES,CF_PENALTIES_gap_extend_cost)
        numCPU = int(CURRENTCONFIG.get(CF_RUN,CF_RUN_num_CPU))
        faFileList = os.listdir(faDir)
        namesList = os.listdir(faDir)
        namesList = [n[0:3] for n in namesList]
        try:
            assert len(set(namesList))==len(namesList)
        except AssertionError as e:
            if reloadFields:
                pass
            elif len(namesList)==0:
                tk.messagebox.showwarning("Error","Check your fasta directory ({}), it seems to be empty.".format(faDir))
                sys.exit(1)
            elif len(namesList)<0:
                tk.messagebox.showwarning("Error","Please make the first three letters of your fasta file names in {} unique.\nSorry for the inconvenience.".format(faDir))
                sys.exit(1)

        #run scythe
        if (not reloadFields or useLocal=="yes") and len(set(namesList)) >0:
            self.progbar = ttk.Progressbar(root, mode='indeterminate')
            self.progbar.grid(column=0, row=5, sticky = "W")
            self.q = multiprocessing.Queue()
            argsdct = {
                    "groups":groups,
                    "delim":delim,
                    "asID":asID,
                    "faFileList":faFileList,
                    "namesList":namesList,
                    "cleanUp":cleanUp,
                    "inDir":inDir,
                    "outDir":outDir,
                    "gapOpen":gapOpen,
                    "gapExtend":gapExtend,
                    "locDir":locDir,
                    "faDir":faDir,
                    "numCPU": numCPU,
                    }
###
            self.p = scythe.ThreadedScythe(self.q, argsdct)
            self.p.start()
            self.progbar.start()
            tk.messagebox.showinfo("Be patient.", "Scythe is running. This could take some time.")
            SCYTHE_PROCESS = self.p
            self.process_queue()

    def process_queue(self):
        if not self.p.is_alive():
            self.progbar.stop()
            tk.messagebox.showinfo("Done.","Done.")
        try:
            msg = self.q.get(0)
            self.progbar.stop()
            #tk.messagebox.showinfo("Process done.",msg)
            sys.exit(0)
        except queue.Empty:
            #print(self.p, self.p.is_alive(), self.q)
            self.parent.after(100, self.process_queue)

    def debug(self,n):
        pass


    def initWizard(self):
        global SCYTHE_PROCESS
        #Labels
        self.lab_fastaDir = tk.Label(text="Fasta Directory")
        self.lab_locDir = tk.Label(text=".loc Directory")
        self.lab_grpFile = tk.Label(text=".grp File")
        self.lab_outDir = tk.Label(text="Output Directory")

        #Ints
        self.int_ensembl = tk.IntVar()
        self.int_local = tk.IntVar()

        #Strings
        self.st_fastaDir = tk.StringVar()
        self.st_fastaDir.set("")
        self.st_locDir= tk.StringVar()
        self.st_locDir.set("")
        self.st_grpFile= tk.StringVar()
        self.st_grpFile.set("")
        self.st_outDir= tk.StringVar()
        self.st_outDir.set("")

        #Entries
        self.ent_fastaDir = tk.Entry(root, width = 30,
                                     textvariable = self.st_fastaDir, state = tk.DISABLED)
        self.ent_locDir = tk.Entry(root, width = 30,
                                     textvariable = self.st_locDir, state = tk.DISABLED)
        self.ent_grpFile = tk.Entry(root, width = 30,
                                     textvariable = self.st_grpFile, state = tk.DISABLED)
        self.ent_outDir = tk.Entry(root, width = 30,
                                     textvariable = self.st_outDir, state = tk.DISABLED)
        #Buttons
        self.b_loadConfig = tk.Button()
        self.b_saveConfig = tk.Button()
        self.b_fastaDir = tk.Button(text="open...",command=self.askFastaDir, state=tk.DISABLED)#self.askopenfilename)
        self.b_locDir = tk.Button(text="open...",command=self.askLocDir, state=tk.DISABLED)#self.askopenfilename)
        self.b_grpFile = tk.Button(text="open...",command=self.askGrpFile, state=tk.DISABLED)#self.askopenfilename)
        self.b_outDir = tk.Button(text="open...",command=self.askOutDir, state=tk.DISABLED)#self.askopenfilename)

        self.b_next = tk.Button(root, text="Next...", command = self.prepRun)
        self.b_quit = tk.Button(root, text="Quit", command = self.quit)
        ######

        #Checkbuttons
        self.cb_use_ensembl = tk.Checkbutton(root, text='use ENSEMBL',variable=self.int_ensembl,
                                         command=self.useEnsembl)
        self.cb_use_local = tk.Checkbutton(root, text='use local files',variable=self.int_local,
                                         command=self.useLocal, state=tk.NORMAL)

        #add to grid
        self.cb_use_ensembl.grid(row=0, column=0 )
        self.cb_use_local.grid(row=0, column=1)

        self.lab_fastaDir.grid(row=1, column=0)
        self.ent_fastaDir.grid(row=1, column=1)
        self.b_fastaDir.grid(row=1, column=2)

        self.lab_locDir.grid(row=2, column=0)
        self.ent_locDir.grid(row=2, column=1)
        self.b_locDir.grid(row=2, column=2)


        self.lab_grpFile.grid(row=3, column=0)
        self.ent_grpFile.grid(row=3, column=1)
        self.b_grpFile.grid(row=3, column=2)

        self.lab_outDir.grid(row=4, column=0, sticky="E")
        self.ent_outDir.grid(row=4, column=1, sticky="E")
        self.b_outDir.grid(row=4, column=2, sticky="E")

        self.b_next.grid(row=5, column=1, sticky="E")
        self.b_quit.grid(row=5, column=2, sticky="W")

    def useLocal(self):
        global CURRENTCONFIG
        if self.int_local.get() ==1:
            CURRENTCONFIG.set(CF_MODE, CF_MODE_use_local_files, "yes")
            CURRENTCONFIG.set(CF_MODE, CF_MODE_use_ensembl, "no")
            self.ent_fastaDir.config(state=tk.NORMAL)
            self.ent_locDir.config(state=tk.NORMAL)
            self.ent_grpFile.config(state=tk.NORMAL)
            self.ent_outDir.config(state=tk.NORMAL)

            self.b_fastaDir.config(state=tk.NORMAL)
            self.b_locDir.config(state=tk.NORMAL)
            self.b_grpFile.config(state=tk.NORMAL)
            self.b_outDir.config(state=tk.NORMAL)

            self.cb_use_ensembl.config(state=tk.DISABLED)

        else:
            CURRENTCONFIG.set(CF_MODE, CF_MODE_use_local_files, "no")
            CURRENTCONFIG.set(CF_MODE, CF_MODE_use_ensembl, "no")

            self.ent_fastaDir.config(state=tk.DISABLED)
            self.ent_locDir.config(state=tk.DISABLED)
            self.ent_grpFile.config(state=tk.DISABLED)
            self.ent_outDir.config(state=tk.DISABLED)

            self.b_fastaDir.config(state=tk.DISABLED)
            self.b_locDir.config(state=tk.DISABLED)
            self.b_grpFile.config(state=tk.DISABLED)
            self.b_outDir.config(state=tk.DISABLED)

            self.cb_use_ensembl.config(state=tk.NORMAL)

    def askFastaDir(self):
        tmp= filedialog.askdirectory()
        self.ent_fastaDir.config(state=tk.NORMAL)
        self.st_fastaDir.set(tmp)
        if self.ent_locDir.get()=="":
            if os.path.isdir(os.path.split(tmp)[0]+os.sep+"loc"):
                self.st_locDir.set(os.path.split(tmp)[0]+os.sep+"loc")
                CURRENTCONFIG.set(CF_PATHS, CF_PATHS_loc_directory, os.path.split(tmp)[0]+os.sep+"loc")
        CURRENTCONFIG.set(CF_PATHS, CF_PATHS_fasta_directory,tmp)
        return tmp

    def askLocDir(self):
        tmp= filedialog.askdirectory()
        print(tmp)
        self.ent_locDir.config(state=tk.NORMAL)
        self.st_locDir.set(tmp)
        if self.ent_locDir.get()=="":
            if os.path.isdir(os.path.split(tmp)[0]+os.sep+"fa"):
                self.st_locDir.set(os.path.split(tmp)[0]+os.sep+"fa")
                CURRENTCONFIG.set(CF_PATHS, CF_PATHS_fasta_directory,os.path.split(tmp)[0]+os.sep+"fa")
        CURRENTCONFIG.set(CF_PATHS, CF_PATHS_loc_directory,tmp)
        return filedialog.tmp()

    def askGrpFile(self):
        tmp= filedialog.askopenfilename()
        self.ent_grpFile.config(state=tk.NORMAL)
        self.st_grpFile.set(tmp)
        CURRENTCONFIG.set(CF_PATHS, CF_PATHS_grp_file,tmp)
        return tmp

    def askOutDir(self):
        tmp= filedialog.askdirectory(mustexist=False)
        self.ent_outDir.config(state=tk.NORMAL)
        self.st_outDir.set(tmp)
        CURRENTCONFIG.set(CF_PATHS, CF_PATHS_output_directory,tmp)

    def useEnsembl(self):
        CURRENTCONFIG.set(CF_MODE, CF_MODE_use_local_files, "no")
        CURRENTCONFIG.set(CF_MODE, CF_MODE_use_ensembl, "yes")
        if self.int_ensembl.get() ==1:
            self.ent_fastaDir.config(state=tk.DISABLED)
            self.ent_locDir.config(state=tk.DISABLED)
            self.ent_outDir.config(state=tk.NORMAL)

            self.b_fastaDir.config(state=tk.DISABLED)
            self.b_locDir.config(state=tk.DISABLED)
            self.b_outDir.config(state=tk.NORMAL)
            self.cb_use_local.config(state=tk.DISABLED)
        else:
            self.cb_use_local.config(state=tk.NORMAL)
            self.b_outDir.config(state=tk.DISABLED)
            self.ent_outDir.config(state=tk.DISABLED)
            CURRENTCONFIG.set(CF_MODE, CF_MODE_use_local_files, "no")
            CURRENTCONFIG.set(CF_MODE, CF_MODE_use_ensembl, "no")

try:
    #use config file
    arg = sys.argv[1]
except IndexError as e:
    #from scratch
    app=ScytheMenu(root)
else:
    app=ScytheMenu(root,arg)

root.mainloop()
