#!/usr/bin/env python

debug = 0
from Load_SOB import *

#outputDir = '/data2/users/jabeen/DATA_2/SB-All'
#resultsdir = '/data2/users/jabeen/DATA_2/SB-All'

# The DATA-2/SB-All/Selected/ will be generated automatically
# when running the source file.

# Change output directory (outputDir) to '/data/users/{username}/Data-2/SB-All'.
outputDir = 'data'

# Change results directory (resultsDir) to '/data/users/{username}/Data-2/SB-All'.
resultsdir = 'data'

parser = ArgumentParser()
parser.add_argument('--baseDirMuG',      default=None,           dest='baseDirMuG',         required=False, help='Path to muon base directory')
#stitparser.add_argument('--baseDirElG',      default=None,           dest='baseDirElG',         required=False, help='Path to electron base directory')
parser.add_argument('--baseDirElG',      default=None,           dest='baseDirElG',         required=False, help='Path to electron base directory')
parser.add_argument('--outputDir',       default=None,           dest='outputDir',          required=False, help='Output directory to write histograms')
parser.add_argument('--data',            default=False,          dest='data',               required=False, help='Use data or MC')
parser.add_argument('--batch',           default=None,          dest='batch',              required=False, help='Supress X11 output')

options = parser.parse_args()

_TREENAME = 'UMDNTuple/EventTree'
_FILENAME = 'tree.root'
#_XSFILE   = 'cross_sections/photon17.py'
_XSFILE   = 'WG_Analysis/Plotting/cross_sections/photon17.py'
_LUMI     = 36000
#_BASEPATH = '/home/jkunkle/usercode/Plotting/LimitSetting/'
#_SAMPCONF = 'Modules/Resonance2017.py'
_SAMPCONF = 'WG_Analysis/Plotting/Modules/Resonance2017.py'
#_SAMPCONF = 'Modules/Resonance.py'



if options.batch:
    ROOT.gROOT.SetBatch(True)
if options.outputDir is not None :
    if not os.path.isdir( options.outputDir ) :
        os.makedirs( options.outputDir )

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(1)


#ROOT.gROOT.SetBatch(True)

# if no option is given, here are the default directories to read
#if options.baseDirMuG is None: options.baseDirMuG = "/data2/users/kakw/Resonances2017/LepGamma_mug_2019_09_15/"
if options.baseDirMuG is None: options.baseDirMuG = "/data2/users/kakw/Resonances2017/LepGamma_mug_2019_10_28/"
#if options.baseDirElG is None: options.baseDirElG = "/data2/users/kakw/Resonances2017/LepGamma_elg_2019_09_15/"
if options.baseDirElG is None: options.baseDirElG = "/data2/users/kakw/Resonances2017/LepGamma_elg_2019_10_28/"
#options.baseDirElG = "/data/users/friccita/WGammaNtuple/LepGamma_elg_2019_04_11/"
#=========Provide the string for base cuts to compare the sigbificanse against in the final SOB plots
#sigstr_BCuts = ["selbase_el_gtmet25_phpt80_elpt40_elidTight_phidTight_invZ10", "selbase_el_gtmet25_phpt80_elpt40_elidTight_phidTight_invZ10", "selbase_el_gtmet25_phpt80_elpt40_elidTight_phidTight_invZ10"]
sigstr_BCuts = "selbase_el_gtmet25_phpt80_elpt40_elidTight_phidTight_invZ10"
    
signal_name = ["MadGraphResonanceMass200_width0p01"
               ,"MadGraphResonanceMass250_width5"
               ,"MadGraphResonanceMass300_width5"
               ,"MadGraphResonanceMass400_width5"
               ,"MadGraphResonanceMass450_width5"
               ,"MadGraphResonanceMass500_width5"
               ,"MadGraphResonanceMass700_width0p01"
               ,"MadGraphResonanceMass800_width0p01"
               ,"MadGraphResonanceMass800_width5"
               ,"MadGraphResonanceMass900_width0p01"
               ,"MadGraphResonanceMass900_width5"
               ,"MadGraphResonanceMass1000_width0p01"
               ,"MadGraphResonanceMass1200_width0p01"
               ,"MadGraphResonanceMass1400_width0p01"
               ,"MadGraphResonanceMass1400_width5"
               ,"MadGraphResonanceMass1800_width5"
               ,"MadGraphResonanceMass2000_width0p01"
               ,"MadGraphResonanceMass2000_width5"
               ,"MadGraphResonanceMass2200_width0p01"
               ,"MadGraphResonanceMass2200_width5"
               ,"MadGraphResonanceMass2600_width0p01"
               ,"MadGraphResonanceMass2800_width0p01"
               ,"MadGraphResonanceMass3500_width5"
               ,"MadGraphResonanceMass4000_width5"]

#sigstr = ["MadGraphResonanceMass250_width5", "MadGraphResonanceMass1000_width0p01", "MadGraphResonanceMass2200_width0p01"]

def main() :
    if options.outputDir: f1 = ROOT.TFile("%s/output.root"%(options.outputDir),"RECREATE")

    #sampManMuG= SampleManager( options.baseDirMuG, _TREENAME, filename=_FILENAME, xsFile=_XSFILE, lumi=_LUMI )
    sampManElG= SampleManager( options.baseDirElG, _TREENAME, filename=_FILENAME, xsFile=_XSFILE, lumi=_LUMI )
    
    #sampManMuG.ReadSamples( _SAMPCONF )
    sampManElG.ReadSamples( _SAMPCONF )
    
    # Set values for cuts.
    cut_phpt = [(phpt70,"phpt70"), (phpt80,"phpt80"), (phpt90,"phpt90"), (phpt100,"phpt100"), (phpt105,"phpt105"), (phpt110,"phpt110"),  (phpt120,"phpt120"), (phpt130,"phpt130"),  (phpt140,"phpt140"),  (phpt150,"phpt150"),  (phpt160,"phpt160"),  (phpt170,"phpt170"),  (phpt180,"phpt180"),  (phpt190,"phpt190"),  (phpt200,"phpt200"),  (phpt210,"phpt210"),  (phpt220,"phpt220"),  (phpt230,"phpt230"),  (phpt240,"phpt240"),  (phpt250,"phpt250")]
    
    cut_elid = [(elidTight,"elidTight"), (elidMedium,"elidMedium"), (elidLoose,"elidLoose")]
    cut_phid = [(phidTight,"phidTight"), (phidMedium,"phidMedium"), (phidLoose,"phidLoose")]
    cut_z = [ (invZ10, "invZ10"), (invZ15, "invZ15"), (invZ20,"invZ20")]
    
    cut_met = [(gtmet25,"gtmet25"), (gtmet30,"gtmet30"), (gtmet40,"gtmet40"), (gtmet50,"gtmet50"),  (gtmet60,"gtmet60"), (gtmet70,"gtmet70"),  (gtmet80,"gtmet80"), (gtmet90,"gtmet90"),  (gtmet100,"gtmet100"),  (gtmet110,"gtmet110"),  (gtmet120,"gtmet120"),  (gtmet130,"gtmet130"),  (gtmet140,"gtmet140"),  (gtmet150,"gtmet150"),  (gtmet160,"gtmet160"),  (gtmet170,"gtmet170"),  (gtmet180,"gtmet180"), (gtmet190,"gtmet190"),  (gtmet200,"gtmet200"),(gtmet210,"gtmet210")]
    
    cut_elpt = [(elpt30,"elpt30"), (elpt40,"elpt40"), (elpt50,"elpt50"), (elpt60,"elpt60"), (elpt70,"elpt70"),  (elpt80,"elpt80"),  (elpt90,"elpt90"), (elpt100,"elpt100"), (elpt110,"elpt110"), (elpt120,"elpt120"), (elpt130,"elpt130"), (elpt140,"elpt140"), (elpt150,"elpt150"), (elpt160,"elpt160")]
    
    #el_pt0_selbase_el_gtmet30_phpt60_elpt160_elidLoose_phidLoose_invZ20_.pdf.log
    
    #Uncomment these instead for debugging selbase_el_gtmet25_phpt80_elpt40_elidTight_phidTight_invZ10
    #cut_met = [(gtmet30,"gtmet30")]
    #cut_elpt = [(elpt160,"elpt160")]
    #cut_phpt = [(phpt60,"phpt60")]
    #cut_elid = [(elidLoose,"elidLoose")]
    #cut_phid = [(phidLoose,"phidLoose")]
    #cut_z = [ (invZ20, "invZ20")]
    
    selarray = [[(selbase_el,"selbase_el"),],  cut_met, cut_phpt, cut_elpt, cut_elid, cut_phid, cut_z]
    #variables to be plotted
    vararray = [ #("el_n",        (10,0,10),      "num of electrons"), ## variable name, x axis range, x axis label
        ("el_pt[0]",    (50,0,500),     "p_{T}(e, leading)"),
        ("el_eta[0]",    (10,-5.0, 5.0),     "#eta (e, leading)"),
                    #("ph_n",        (10,0,10),      "num of photons"), 
        ("ph_pt[0]",    (50,0,500),     "p_{T}(#gamma, leading)"),
        ("ph_eta[0]",    (10,-5.0,5.0),     "#eta (#gamma, leading)"),
        ("met_pt",      (50,0,500),     "MET"),
        ("met_phi",    (20,-pi,pi),    "MET #phi")
        ] 
   
    
    #    legend_config = {'legendLoc':"Double","legendTranslateX":0.3}
    hist_config = {"logy":1,"blind":True, "weight": "PUWeight*NLOWeight"}
   
    # ========steps for finding optimal cuts ================ 

# Do the following steps when first running the file:
    # a) Set getyields to 1. This will apply the initial cuts and save all relevant background 
    #    and signal numbers into log files.
    
    # b) Set makesob to 0. It is necessary to create the relevant log files to be able to make
    #    the plots.
    
    # c) Uncomment lines 114-119. This will replace the extensive signal name list with a much 
    #    shorter version, for the purpose of debugging. 
    
    # d) Run the ./run_MakeVarDist_and_SB.py --batch BATCH command to make the log files.

    # e) Reset getyields and makesob in main to 0 and 1, respectively. This will generate the plots 
    #    based on the recently created log files. 

    #STEP 1 -  apply cuts and save final yields for all background and signal samples
    getyields = 1
    
#STEP 2  get numbers from above saved log files, make tables and calculate significanse. Finally make SOB plots for every signal sample 
    makesob = 1

#STEP 3 get the cut strings from tables saved in the above step and plot variables corresponding to those cuts
    drawvars = 1


    if (getyields):
        vararray = [ ("el_pt[0]",    (50,0,200),     "p_{T}(e, leading)")] 
    
        makeplots(0, sampManElG, vararray, resultsdir, selarray, hist_config, {}, "")
        # first 0 means dont save the plots as we just need logfiles
                #legend_config = {'legendLoc':"Double","legendTranslateX":0.3}
                #hist_config = {"blind":True, "weight": "PUWeight*NLOWeight"}
                #makeplots(vararray, selarray, hist_config, legend_config)    
    
    if(makesob):
        vararray = [ ("el_pt[0]",    (50,0,200),     "p_{T}(e, leading)")] 
    
        for j in range(len(signal_name)):
            makesob_plots(0,  sigstr_BCuts, sampManElG,resultsdir,vararray, signal_name[j], selarray, hist_config,{}, "" )# first param 0 to use already existing table rather than running reading from all the log files again
#            makesob_plots(1,  sigstr_BCuts, sampManElG,resultsdir,vararray, signal_name[j], selarray, hist_config,{}, "" )# first param 0 to use already existing table rather than running reading from all the log files again

    n = 0
    if(drawvars):
        
        for j in range(len(signal_name)):
            if (n == 0):
                legend_config = {'legendLoc':"Double","legendTranslateX":0.3}
            else:
                legend_config = {'legendLoc':"Double","legendTranslateX":0.95}
            n = n +1
            hist_config = {"logy":1,"blind":True, "weight": "PUWeight*NLOWeight"}
            
            make_selected_plots(sigstr_BCuts,sampManElG,resultsdir, vararray, signal_name[j],  selarray, hist_config, {}, "log")
       
            hist_config = {"blind":True, "weight": "PUWeight*NLOWeight"}
   
        
            make_selected_plots(sigstr_BCuts,sampManElG,resultsdir, vararray, signal_name[j],  selarray, hist_config, {}, "")
       




    if options.outputDir:
       ## write and close root file
        f1.Write()
        f1.Close()

    

main()
