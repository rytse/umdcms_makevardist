#####!/usr/bin/env python
import os
import sys
import subprocess
# import matplotlib
# matplotlib.use('agg')
# import matplotlib.pyplot as plt
# import pylab as pl
# import pandas as pd
import math

# Import modules from ROOT.

import ROOT, os, re, string, csv, math, sys
from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH1D, TH2F, TF1, TGaxis, TPad, TPaveLabel, TPaveText, TLegend, \
    TLatex, THStack, TLine, TMath, TGraph, TGraphErrors, TMultiGraph, TStyle, TLeaf
from ROOT import gROOT, gBenchmark, gRandom, gSystem, Double, gPad, gStyle
from ROOT import TCanvas, TGraph, gROOT, gPad
from array import array
from itertools import product

ROOT.PyConfig.IgnoreCommandLineOptions = True
import numpy as np
from math import pi
import selection_defs as defs
# from FitManager import FitManager
from SampleManager import SampleManager

from argparse import ArgumentParser
from itertools import product

debug = 0
# electron base selection:
# ph_n==1 && el_n==1 && mu_n==0 && ph_IsEB[0]&& &&abs(el_eta[0])<2.1)*PUWeight*NLOWeight
# el_pt30_n==1 &&
# met_pt>25&&
# && abs(m_lep_ph-91)>15
# && ph_hasPixSeed[0]==0&1
# &&ph_pt[0]>80
# &&el_passTight[0]
# ph_passMedium[0]
# &&el_pt[0]>40
# with bad tracker portions excluded

# Set cuts for event processing.

selbase_el_excludephi = 'ph_n>=1&&el_n==1&&!( ph_eta[0]<0&&ph_phi[0]>2.3&&ph_phi[0]<2.7)&&!(ph_phi[0]>1.2&&ph_phi[0]<1.5)'
selbase_mu = 'mu_pt30_n==1&& mu_n==1'
# selbase_el = 'ph_n>=1&&el_n==1'

selbase_el = "ph_n==1&&el_n==1&&mu_n==0&&ph_IsEB[0]&&abs(el_eta[0])<2.1&&ph_hasPixSeed[0]==0&1"

## general event selections
nocut = ("", "")

# ltmet25 = 'met_pt<25'
gtmet25 = 'met_pt>25'  # analysis cut
gtmet30 = 'met_pt>30'
gtmet35 = 'met_pt>35'
gtmet40 = 'met_pt>40'
gtmet45 = 'met_pt>45'
gtmet50 = 'met_pt>50'
gtmet55 = 'met_pt>55'
gtmet60 = 'met_pt>60'
gtmet65 = 'met_pt>65'
gtmet70 = 'met_pt>70'
gtmet75 = 'met_pt>75'
gtmet80 = 'met_pt>80'
gtmet85 = 'met_pt>85'
gtmet90 = 'met_pt>90'
gtmet95 = 'met_pt>95'
gtmet100 = 'met_pt>100'
gtmet105 = 'met_pt>105'
gtmet110 = 'met_pt>110'
gtmet115 = 'met_pt>115'
gtmet120 = 'met_pt>120'
gtmet125 = 'met_pt>125'
gtmet130 = 'met_pt>130'
gtmet135 = 'met_pt>135'
gtmet140 = 'met_pt>140'
gtmet145 = 'met_pt>145'
gtmet150 = 'met_pt>150'
gtmet155 = 'met_pt>155'
gtmet160 = 'met_pt>160'
gtmet165 = 'met_pt>165'
gtmet170 = 'met_pt>170'
gtmet175 = 'met_pt>175'
gtmet180 = 'met_pt>180'
gtmet190 = 'met_pt>190'
gtmet195 = 'met_pt>195'
gtmet200 = 'met_pt>200'
gtmet210 = 'met_pt>210'
# cut_met = [(gtmet25,"gtmet25"), (gtmet20,"gtmet20"), (gtmet30,"gtmet30"), (gtmet35,"gtmet35"), (gtmet40,"gtmet40"),(gtmet45,"gtmet45"), (gtmet50,"gtmet50"), (gtmet55,"gtmet55"), (gtmet60,"gtmet60"), (gtmet65,"gtmet65"), (gtmet70,"gtmet70"), (gtmet75,"gtmet75"), (gtmet80,"gtmet80"),(gtmet85,"gtmet85"), (gtmet90,"gtmet90"), (gtmet95,"gtmet95"), (gtmet100,"gtmet100"), (gtmet105,"gtmet105"), (gtmet110,"gtmet110"), (gtmet115,"gtmet115"), (gtmet120,"gtmet120"), (gtmet125,"gtmet125"), (gtmet130,"gtmet130"), (gtmet135,"gtmet135"), (gtmet140,"gtmet140"), (gtmet145,"gtmet145"), (gtmet150,"gtmet150"), (gtmet155,"gtmet155"), (gtmet160,"gtmet160"), (gtmet165,"gtmet165"), (gtmet170,"gtmet170"), (gtmet175,"gtmet175"), (gtmet180,"gtmet180"), (gtmet190,"gtmet190"), (gtmet195,"gtmet195"), (gtmet200,"gtmet200"),(gtmet205,"gtmet205")]


elpt30 = "el_pt[0]>30"
elpt40 = "el_pt[0]>40"  # analysis cut
elpt50 = "el_pt[0]>50"
elpt60 = "el_pt[0]>60"
elpt65 = "el_pt[0]>65"
elpt70 = "el_pt[0]>70"
elpt75 = "el_pt[0]>75"
elpt80 = "el_pt[0]>80"
elpt85 = "el_pt[0]>85"
elpt90 = "el_pt[0]>90"
elpt100 = "el_pt[0]>100"
elpt110 = "el_pt[0]>110"
elpt120 = "el_pt[0]>120"
elpt130 = "el_pt[0]>130"
elpt140 = "el_pt[0]>140"
elpt150 = "el_pt[0]>150"
elpt160 = "el_pt[0]>160"

phpt60 = "ph_pt[0]>60"
phpt70 = "ph_pt[0]>70"
phpt80 = "ph_pt[0]>80"  # analysis cut
phpt90 = "ph_pt[0]>90"
phpt100 = "ph_pt[0]>100"
phpt105 = "ph_pt[0]>105"
phpt110 = "ph_pt[0]>110"
phpt115 = "ph_pt[0]>115"
phpt120 = "ph_pt[0]>120"
phpt130 = "ph_pt[0]>130"
phpt140 = "ph_pt[0]>140"
phpt150 = "ph_pt[0]>150"
phpt160 = "ph_pt[0]>160"
phpt170 = "ph_pt[0]>170"
phpt180 = "ph_pt[0]>180"
phpt190 = "ph_pt[0]>190"
phpt200 = "ph_pt[0]>200"
phpt210 = "ph_pt[0]>210"
phpt220 = "ph_pt[0]>220"
phpt230 = "ph_pt[0]>230"
phpt240 = "ph_pt[0]>240"
phpt250 = "ph_pt[0]>250"
# cut_phpt = [(phpt70,"phpt70"), (phpt80,"phpt80"), (phpt90,"phpt90"), (phpt100,"phpt100"), (phpt105,"phpt105"), (phpt110,"phpt110"),  (phpt120,"phpt120"), (phpt130,"phpt130"),  (phpt140,"phpt140"),  (phpt150,"phpt150"),  (phpt160,"phpt160"),  (phpt170,"phpt170"),  (phpt180,"phpt180"),  (phpt190,"phpt190"),  (phpt200,"phpt200"),  (phpt210,"phpt210"),  (phpt220,"phpt220"),  (phpt230,"phpt230"),  (phpt240,"phpt240"),  (phpt250,"phpt250")]


elidTight = 'el_passTight[0]'
elidMedium = 'el_passMedium[0]'
elidLoose = 'el_passLoose[0]'

phidTight = 'ph_passTight[0]'
phidMedium = 'ph_passMedium[0]'
phidLoose = 'ph_passLoose[0]'

invZ15 = 'abs(m_lep_ph-91)>15'  ## analysis inverse Z resonance cut
invZ10 = 'abs(m_lep_ph-91)>10'  ##
invZ20 = 'abs(m_lep_ph-91)>20'  ##
inZ = 'abs(m_lep_ph-91)<15'  ## Z resonance cut
gtZ = '(m_lep_ph-91)>15'  ## greater than Z mass
ltZ = '(91-m_lep_ph)>15'  ## less than Z mass


# cut_z = [ (invZ10, "invZ10"), (invZ15, "invZ15")]
# selbase_el_EB = selbase_el + ' && ph_IsEB[0]' #leading photon in barrel

# For testing uncomment these instaed:

# sigstr_BCuts = ["selbase_el_gtmet25_phpt80_elpt40_elidTight_phidTight_invZ10", "selbase_el_gtmet25_phpt80_elpt40_elidTight_phidTight_invZ10", "selbase_el_gtmet25_phpt80_elpt40_elidTight_phidTight_invZ10"]

# sigstr = ["MadGraphResonanceMass250_width5", "MadGraphResonanceMass1000_width0p01", "MadGraphResonanceMass2200_width0p01"]


def makeplots(dostack, samples, vararray, dirname, selprearray=None, hist_config=None, legend_config=None, extratag=""):
    """
    Sample from each combination of cuts specified in `selprearray`, saving the samples to log files and histogram plots
    for each cut. Note that this dumps data for every signal (independent variable).

    :param dostack: whether to save the plot files
    :param samples: `SampleManager` containing the (full-dimensional) simulated data
    :param vararray: variables of interest (dependent variables) to log/plot
    :param dirname: directory to save logs and/or histograms to
    :param selprearray: zipped list of cuts (see `makeselection` for format)
    :param hist_config: dict of `ROOT` plotting options
    :param legend_config: legend plotting options
    :param extratag: log file ending string (optional)
    """
    # print 'in makeplots'
    if (debug):
        print selprearray
    chlist = "[]{}()"  # chars used in SampleManager cut names that are excluded in output that humans need to read
    if selprearray is None or not isinstance(selprearray, (tuple, list)):
        return
    else:
        # Make array a single tuple/product. ach
        selarray = product(*selprearray)

    if debug:
        print str(selarray)
    for sel in selarray:  # iterate over every possible set of cuts
        # print "sel in sellarray" + str(sel)
        selection, name = makeselection(sel)  # format cut string
        if debug:
            print "selectionl = " + str(selection)
        # print  "name = " + str(name)
        # print str(selection)
        nl = 0
        for var in vararray:
            # Open output files
            savename = var[0] + "_" + name + "_" + extratag + ".pdf"
            for ch in chlist:  # strip out the human unreadable stuff ("[]{}()")
                savename = savename.replace(ch, "")
            sys.stdout = open(dirname + '/' + savename + '.log', 'w')

            # Configure plot options
            hist_config["xlabel"] = var[2]
            # print var[0], selection, savename
            if nl == 0:
                legend_config = {'legendLoc': "Double", "legendTranslateX": 0.3}
            else:
                legend_config = {'legendLoc': "TopLeft", "legendTranslateX": 0.95, "legendTranslateY": 0.0}
            nl = nl + 1

            # Draw from the SampleManger "DB" given the chosen cuts
            samples.Draw(var[0], selection, var[1], hist_config, legend_config)

            # Prints out number of events for every signal and background
            # sample (i.e., all combinations of signals and backgrounds). ach
            samples.print_stack_count()
            if dostack:  # plot hist
                samples.SaveStack(savename, dirname, 'base')


# Python manipulation. Adds && to selection array. ach
def makeselection(sel):
    """
    Converts a cut string from zipped format into a `SampleManager`-friendly format.

    We get each cut that we plot from a generator (`product(*selprearray)`) that yields zipped lists of tuples of the
    following form:

        `sel = [(<cut range>, <cut variable name>), (<cut range>, <cut variable name>), ...]`

    This function unzips this into a list of cut ranges and a list of cut variable names, and reformats these strings in
    the following way to conform to `SampleManager`'s expected format:
        * Prepend "_" to names
        * Prepend "&&" to selections and strip whitespace

    :param sel: zipped list of cut tuples generated by `product(*selprearray)`
    :return: (<list of cut ranges>, <list of cut variable names>) formatted the way `SampleManager` expects selections
    """
    # print 'in make selection sel ='
    if debug:
        print sel
    sel, name = zip(*sel)
    # print "sel =" + str(sel)
    # print "name = " + str(name)
    sel = [s for s in sel if s != ""]
    name = [n for n in name if n != ""]
    sel = "&&".join(sel)
    sel.replace(" ", "")
    name = "_".join(name)
    # print sel
    # print name
    return sel, name


# Calculates significance.
def signif(samples, vararray, signal, dirname, selprearray=None, hist_config=None, legend_config=None, extratag=""):
    """
    Load samples logged by `makeplots` and extract the number of background and signal events, as well as their
    respective error estimates.

    These values summary values are are dumped to table files for each variable and each cut, as well as summarized in
    one big table placed in the `Selected` directory.

    :param samples: `SampleManager` containing the (full-dimensional) simulated data
    :param vararray: variables of interest (dependent variables) to aggregate counts and error estimates for
    :param signal: name of signal (independent variable) to count
    :param dirname: name of directory to save counts to
    :param selprearray: zipped list of cuts (see `makeselection` for format)
    :param hist_config: dict of `ROOT` plotting options
    :param legend_config: legend plotting options
    :param extratag: log file ending string (optional)
    :return:
        - xname - "human readable" cut names
        - total_b - total number of background events
        - total_b_err - error bars on number of background events
        - total_s - total number of signal events
        - total_s_err - error bars on total number of signal events
        - sig_cut_for_plots - cut values
    """
    # print "in signf ==========================="
    # print selprearray
    chlist = "[]{}()"
    data = {}  # TODO remove unused

    # Total background.
    total_b = []

    # Total background error.
    total_b_err = []

    # Total signal.
    total_s = []

    # Total signal error.
    total_s_err = []

    # total_s_list  = [[] for i in range(len(sigstr))]

    # Significance cut for plots (i.e., the primary cut array)
    sig_cut_for_plots = []
    xname = []  # cutname to appear on x axis

    if selprearray is None or not isinstance(selprearray, (tuple, list)):
        return
    else:
        selarray = product(*selprearray)

    # print signal
    numselection = -1
    for sel in selarray:
        numselection = numselection + 1
        # print "Number of selection" + str(numselection)
        # if (numselection > 200):
        #    continue
        selection, name = makeselection(sel)
        # print selection
        # print name
        for var in vararray:
            savename = var[0] + "_" + name + "_" + extratag + ".pdf"
            # print savename
            for ch in chlist:
                savename = savename.replace(ch, "")
            os.system("grep +/- " + str(dirname + '/' + savename + '.log  > tmp.f'))
            tf = os.stat('tmp.f')
            x = tf.st_size
            # print x

            # Load the logfile corresponding to the specified cut and variable of interest and dump the values of
            # interest into table files
            if os.path.isfile(dirname + '/' + savename + '.log') and x:
                # f_table = dirname + '/' + savename + "_table"
                f_table = dirname + "/" + signal + savename + "_table"
                sigrep = "Make " + signal + " hist"  # string representing the signal file
                lookfor = ["TOTAL", sigrep, signal]  # headers to look for in logfiles
                # print lookfor

                # TODO don't open f_table if you aren't going to use it ...
                with open(f_table, 'w+') as f:  # W+ creates the new file for writing
                    nline = 0  # TODO remove unused
                    for i in range(len(lookfor)):
                        # print lookfor[i]

                        # grep to dump column headers into the table
                        command = "grep \"" + str(lookfor[i]) + "\" " + str(  # grep command used to find the file
                            dirname + '/' + savename + '.log |  grep +/- >> ') + str(f_table)
                        os.system(command)

                        # grep to dump signal values into the table
                        if lookfor[i] == sigrep:
                            command = "grep \"" + sigrep + "\" " + str(dirname + '/' + savename + '.log  >> ') + str(
                                f_table)
                            # print command
                            os.system(command)

                f = open(f_table)  # after all the data is dumped into the file, load it up
            else:
                continue

            # except OSError:
            # continue

            # Reload all that data you just dumped into the table file into arrays
            lines = f.readlines()
            # print "read first line of the table we just saved" + str(lines[0])
            for line in lines:
                listline = line.split()
                # print line

                # Grab cut names
                if 'Make' in line.split():
                    # print line
                    sigcut = str(line.split('(', 1)[1])  # this selects the part after '(
                    sigcut.replace(" ", "")
                    sig_cut_for_plots.append(sigcut)
                    # print sigcut

                # Grab background event counts and their uncertainty values
                if 'TOTAL' in line and '+/-' in line.split():
                    # print line
                    # print "line type = " + str(type(line))
                    # y_total = float(listline[1])
                    # y_total_err = float(listline[3])
                    total_b.append(float(listline[1]))
                    total_b_err.append(float(listline[3]))
                    # print total_b, total_b_err

                # Grab signal event counts and their uncertainty values
                if signal in line.split() and '+/-' in line.split():
                    total_s.append(float(listline[1]))
                    total_s_err.append(float(listline[3]))
                    # print "y_sig[i] ============" + str(y_sig[i])
            # print "appending ===" + str(name)
            xname.append(name)
            f.close()
    # print '%20s  %15s  = %5s      Backgd=  %5s' %( name,   sigstr[i],  str(round(y_sig[i],2)),str( round(total_b[0],2))  )
    # for n in range(len(sigstr)):

    # Log the background, background error est, signal, signal error est values for each cut + variable of interest to
    # a table in the "Selected" subdir
    seldir = dirname + "/Selected"
    if not os.path.exists(seldir):  # If the directory does not exist, the operating system will auto generate it
        os.makedirs(seldir)
    ftable = seldir + "/" + signal + "TOTAL_table"
    with open(ftable, 'w+') as f:
        f.write("cutname total_b total_b_err total_s total_s_err  cut_array_used\n")
        # f.write("xname  total_b total_b_err total_s_list total_s_err_list sigstr sig_cut_for_plots\n")
        # print sigstr[n]

        print "\n\n\n\nDEBUGGING CURRENT CUT\n"
        print "signal: " + str(signal)
        print "xname: " + str(xname)
        print "total_b: " + str(total_b)
        print "total_b_err: " + str(total_b_err)
        print "total_s: " + str(total_s)
        print "total_s_err: " + str(total_s_err)
        print "sig_cut_for_plots: " + str(sig_cut_for_plots)
        print "\n\n\n\n"

        for j in range(len(xname)):
            # print xname[j]
            # print('%15s ----  %5s ----   %5s --- %5s --- %5s --- %5s' %(str(xname[j]),  str(total_b[j]), str(total_b_err[j]) , str(total_s[j]) , str(total_s_err[j]) , sig_cut_for_plots[j]))
            f.write('%15s   %5f    %5f  %5f  %5f  %5s' % (
                xname[j], total_b[j], total_b_err[j], total_s[j], total_s_err[j], sig_cut_for_plots[j]))
    # row_format = '{:>15}' * len(xname)
    # for v in zip(*xname):
    #    if (debug): print (row_format.format(*v))
    f.close()
    return xname, total_b, total_b_err, total_s, total_s_err, sig_cut_for_plots


# Make signal over background plots.

def makesob_plots(calcsig, sigstr_BCuts, samples, dirname, vararray, signal, selarray=None, hist_config=None,
                  legend_config=None, extratag=""):
    """
    Filter cuts by a significance threshold, plot signal-over-background (SOB) for each cut that passed the
    threshold, and save signal and background counts and errors for each "significant" cut into a table.

    SOB is defined as $\frac{S}{\sqrt{B}}$.

    :param calcsig: whether to recalculate tables of signal and background event counts and errors using `signif`
    :param sigstr_BCuts: list of background cut names
    :param samples: `SampleManager` containing the (full-dimensional) simulated data
    :param dirname: directory to save logs and histograms to
    :param vararray: variables of interest (dependent variables) to aggregate counts and error estimates for
    :param signal: name of signal (independent variable) to count
    :param selarray: list of cuts (see `makeselection` for format)
    :param hist_config: dict of `ROOT` plotting options
    :param legend_config: legend plotting options
    :param extratag: log file ending string (optional)
    """
    # Set up ROOT plot
    c1 = TCanvas('c1', 'Significance', 200, 10, 1500, 700)
    c1.SetBottomMargin(0.4)
    c1.SetRightMargin(0.18)
    c1.SetLeftMargin(0.1)
    # p1 = TPad("p1","p1",0.0,0.0,1.0,1.0,0)#TPad(name,     title,     xlow,     ylow,    xup,    yup, color )
    # p1.SetBottomMargin(1.7)
    # p1.GetFrame().SetBorderSize(5)
    # p1.SetRightMargin(0.01)#percent of the pad height
    # p1.SetLeftMargin(0.1)
    # p1.SetBottomMargin(0.5)
    # p1.Draw()

    # get the cut set name and corresponding significance
    xname_all = []  # cut names for x labels
    # signame = []
    sig_all = []
    sig_err_all = []
    bkg_all = []
    bkg_err_all = []
    sigstr_BaseCuts = sigstr_BCuts  # base cut string
    basecut_listvalue = 0  # index number of the base cut in the xname_all list
    sig_cut_plots = []

    # xname, ysig, signame = signif(sampManElG,vararray, selarray, hist_config,{}, "log")

    # If calcsig != 0 ach
    # If recalculating from scratch, call signif to count background and signal events
    if calcsig:
        xname_all, bkg_all, bkg_err_all, sig_all, sig_err_all, sig_cut_plots = signif(samples, vararray, signal,
                                                                                      dirname, selarray, hist_config,
                                                                                      {}, "")
    # If loading background and signal events from a file already logged by signif, then load them
    else:
        ftable = dirname + "/Selected/" + signal + "TOTAL_table"
        n = 0
        with open(ftable, 'r') as f:
            lines = f.readlines()
            # print "read first line of the table we just saved" + str(lines[1])
            for line in lines:
                n = n + 1
                listline = line.split()
                # print line
                if n > 1:
                    # print listline[5]

                    # Add background and significance to their respective
                    # arrays for the plot generation.
                    xname_all.append(str(listline[0]))
                    bkg_all.append(float(listline[1]))
                    bkg_err_all.append(float(listline[2]))
                    sig_all.append(float(listline[3]))
                    sig_err_all.append(float(listline[4]))
                    sig_cut_plots.append(str(listline[5]))
                    # print sig_all[0]

        f.close()

    # if (debug):
    # print "The bigger lists : sig_all[0]= " + str(sig_all[0])+ "bkg_all[1]=" + str(bkg_all[1])
    print "size of bigger array is " + str(len(sig_all))
    # check the list of lists is correct - update: it is ot a list of lists anymore but leaving for future.
    # row_format = '{:<4}' * len(sig_all)
    # for v in zip(*sig_all):
    #    #if (debug):
    #   print (row_format.format(*v))
    # print sig_all

    # Cuts and associated values that are deemed significant by Punzi's criterion
    xname = []
    # signame = []
    sig = []  # it is alist of lists for significance for all the signal samples
    sig_err = []  # it is alist of lists for significance for all the signal samples
    # bkg = [[] for i in range(len(signame))]
    bkg = []
    bkg_err = []
    sig_cut_for_plots = []  # exact string of selection cuts to be used later to make plots
    histall = []

    sig_max = 0.
    bkg_max = 0.

    # for j in range(len(signame)):# loop on signal samples
    sig_max = sig_all[0]
    bkg_max = bkg_all[0]
    signif_max = sig_max / bkg_max
    # if (debug):
    print "signif_max  " + str(signif_max)
    # sig[j].append(sig_all[j][0])
    sig.append(sig_max)
    bkg.append(bkg_all[0])
    sig_err.append(sig_err_all[0])
    bkg_err.append(bkg_err_all[0])
    xname.append(xname_all[0])
    sig_cut_for_plots.append(sig_cut_plots[0])
    # if (debug):
    # print "j   " + str(j) + "   " + str(sig_all[j][0]/bkg_all[0])
    for i in range(1, len(xname_all)):  # loop on cuts sets
        signif_cut = sig_all[i] / math.sqrt(bkg_all[i]) # compute variance assuming
        # print sigstr_BaseCuts
        # print "cut i : " + str(xname_all[i])
        if xname_all[i] == sigstr_BaseCuts:   # if the cut is a base cut, print
            print "Base cut for this signal is : " + xname_all[i]
            # print  " ith bin in xname_all[i]  " + str(i)
            basecut_listvalue = i
            # if (debug):
            print "i   " + str(i) + "  " + xname_all[i] + "    " + str(signif_cut)
        if signif_cut > (1.015 * signif_max): # append significant events to list of relevant counts
            print "i   " + str(i) + "  " + xname_all[i] + "    " + str(signif_cut)
            signif_max = signif_cut
            # if (debug):
            print " new signif_max  " + str(signif_max) + "   " + str(i)
            sig.append(sig_all[i])
            bkg.append(bkg_all[i])
            sig_err.append(sig_err_all[i])
            bkg_err.append(bkg_err_all[i])
            xname.append(xname_all[i])
            sig_cut_for_plots.append(sig_cut_plots[i])
            # signame.append(signame_all[i+1])
    if debug:
        print sig_all[1]
        print sig[1]
        print "How many cuts ?   " + str(len(xname[0]))

    # n= int(len(xname))
    # nsig= int(len(signame))
    nhist = ["_sig", "_bkg", "_bkgsqrt", "_sob", "_sosqrtb"]

    # Bin 0 contains the underflow
    # Bin 1 contains the first bin with low-edge ( xlow INCLUDED).
    # The second to last bin (bin# nbins) contains the upper-edge (xup EXCLUDED).
    # The Last bin (bin# nbins+1) contains the overflow.

    # Make the SOB plot for ALL cuts
    histName_All = signal + "_h_SOB_ALL"
    nall = int(len(xname_all))
    histNameAll = ROOT.TH1F(histName_All, "", nall + 1, 0, nall + 1)
    for i in range(len(xname_all)):  # loop on cuts sets
        histNameAll.SetBinContent(i + 1, sig_all[i] / (math.sqrt(bkg_all[i])))
        # histNameAll.SetBinError(i+1, sig_err[i])#
    histNameAll.Draw()
    histNameAll.GetYaxis().SetTitle('S/#\sqrt B')
    histNameAll.GetXaxis().SetTitle('Cut combination')
    savename_all = str(dirname + "/Selected/" + signal + "_SOsqrtB_ALL")
    c1.SaveAs(savename_all + ".png")
    c1.SaveAs(savename_all + ".pdf")

    c1.Modified()
    c1.Update()

    #################################################
    # for i in range(len(signame)):
    if debug:
        print "Histogram " + signal
    n = int(len(xname))
    for j in range(len(nhist)):
        histName = signal + "_h" + nhist[j]
        # if (debug):
        # print histName
        histall.append(ROOT.TH1F(histName, "", n + 1, 0, n + 1))
        # if (debug):
        # print"histogram name is = " + histall[j].GetName()
    # for j in range(len(signame)):# loop on signal samples
    if debug:
        print signal
    n = int(len(xname))
    bscut = basecut_listvalue
    print " basecut_listvalue = " + str(basecut_listvalue)
    for k in range(len(nhist)):
        for i in range(len(xname)):  # loop on cuts sets
            if (debug):
                print "signal  and cut set  " + str(j) + "  " + str(i)
                print xname[i] + "has significance" + str(ysig[i])  # first index is signal sample, second index is cut
                print str(i) + "   sig[i]   " + str(sig[i])
                print nhist[k]
            # sig = round(float(ysig[i]),2)
            # if (debug):
            # print histall[k].GetName()

            histall[k].GetXaxis().SetBinLabel(i + 1, xname[i])
            # print "bscut  =" + str(bscut)
            # print "length of xname_all[j] arra  =" + str(len(xname_all))
            # print "xname_all[j][bscut]    " + xname_all[bscut]
            histall[k].GetXaxis().SetBinLabel(n + 1,
                                              "X" + xname_all[bscut])  # notice it is from the initial, bigger cutlist

            # print sig_cut_plots[i]
            # print "base cut value = " + str(bscut) + str(sig_all[j][bscut])
            if nhist[k] == '_sig':
                histall[k].SetBinContent(i + 1, sig[i])
                histall[k].SetBinError(i + 1, sig_err[i])  #
                histall[k].SetBinContent(n + 1, sig_all[bscut])
                histall[k].SetBinError(n + 1, sig_err_all[bscut])  #

                # histall[j][k].GetXaxis().SetBinLabel(i+1, xname[j][i])
                # print "n+1 bin label = " + histall[j][k].GetXaxis().GetBinLabel(n+1)
                # if (debug):
                # print "n+1 bin content " + str(histall[k].GetBinContent(n+1))
            if nhist[k] == '_bkg':
                histall[k].SetBinContent(i + 1, bkg[i])
                histall[k].SetBinError(i + 1, bkg_err[i])  #
                histall[k].SetBinContent(n + 1, bkg_all[bscut])
                histall[k].SetBinError(n + 1, bkg_err_all[bscut])  #
                if (debug):
                    print histall[k].GetBinContent(i + 1)
            if nhist[k] == '_bkgsqrt':
                percenterr = (bkg_err[i] / bkg[i]) / 2
                berr = percenterr * math.sqrt(bkg[i])
                berr_bscut = ((bkg_err_all[bscut] / bkg_all[bscut]) / 2) * math.sqrt(bkg_all[bscut])
                histall[k].SetBinContent(i + 1, math.sqrt(bkg[i]))
                histall[k].SetBinError(i + 1, berr)  #
                histall[k].SetBinContent(n + 1, math.sqrt(bkg_all[bscut]))
                histall[k].SetBinError(n + 1, berr_bscut)  #
                if debug:
                    print histall[k].GetBinContent(i + 1)
                    print Sigh.GetXaxis().GetBinLabel(i + 1)

    gStyle.SetOptStat(0);
    # for j in range(len(signame)):
    for k in range(len(nhist)):
        if nhist[k] == "_sob":
            histall[k] = doratio(histall[0], histall[1])
            histall[k].GetYaxis().SetTitle('S/B')
            if debug:
                print "sob" + str(histall[k].GetBinContent(1))
        if nhist[k] == "_sosqrtb":
            histall[k] = doratio(histall[0], histall[2])
            histall[k].GetYaxis().SetTitle('S/#\sqrt B')
        if nhist[k] == "_sig":
            histall[k].GetYaxis().SetTitle('Signal')
        if nhist[k] == "_bkg":
            histall[k].GetYaxis().SetTitle('Total Background')
        if nhist[k] == "_bkgsqrt":
            histall[k].GetYaxis().SetTitle('Sqrt Total Background')

        histall[k].SetLineColor(k + 2)
        histall[k].SetLineWidth(4)
        # histall[k].GetXaxis().LabelsOption("v")
        histall[k].GetXaxis().SetLabelSize(0.025)  # percentage of pad size
        # histall[k].GetXaxis().SetTitleSize(0.10)
        # histall[k].GetXaxis().SetTitleOffset(1.0)
        histall[k].GetYaxis().SetTitleSize(0.06)
        histall[k].GetYaxis().SetTitleOffset(0.8)
        histall[k].GetYaxis().SetLabelSize(0.04)
        histall[k].GetYaxis().CenterTitle()
        histall[k].GetYaxis().SetNdivisions(510, True)
        if debug:
            print histall[k].GetXaxis().GetBinLabel(i + 1)

    maxy = 0.0
    miny = 0.0
    # for j in range(len(signame)):
    sn = signal
    # sob bin with max value
    sob_max_bin = 0
    for k in range(len(nhist)):
        maxvalue = histall[k].GetBinContent((histall[k].GetMaximumBin()))
        minvalue = histall[k].GetBinContent((histall[k].GetMinimumBin()))
        maxy = maxvalue
        miny = minvalue
        if maxvalue > maxy:
            maxy = maxvalue
        if minvalue > miny:
            miny = minvalue

        histall[k].SetMaximum(maxvalue + 0.5 * maxvalue)
        histall[k].SetMinimum(minvalue - 0.5 * minvalue)
        if k == 4:
            sob_max_bin = histall[k].GetMaximumBin()
        # sigstr_BaseCut[]
        # = histall[k].SetBinContent((n+a, histall[k].GetXaxis().GetBinLabel(i+1)))

    # for j in range(len(signame)):
    sn = signal

    for k in range(len(nhist)):
        if k == 2:
            continue
        c1.Modified()
        c1.Update()
        if k == 1:
            gPad.SetLogy()
        else:
            gPad.SetLogy(0)
        legend = TLegend(0.1, 0.7, 0.5, 0.89, " ")  # x1,y1,x2,y2
        maxcut = histall[k].GetXaxis().GetBinLabel(sob_max_bin)
        maxvalue = round(histall[k].GetBinContent(sob_max_bin), 2)
        # legend. AddEntry(histall[j][k] , signame[j] + ' (' + maxcut + ')' , "")
        print signal
        # p1.cd()
        histall[k].Draw()
        c1.Modified()
        c1.Update()
        BC = 0
        BC = round(histall[k].GetBinContent(n + 1), 2)
        print str(c1.GetUxmax()) + "  " + str(BC)
        # p1.Update();
        l = ROOT.TLine(c1.GetUxmin(), BC, c1.GetUxmax(), BC)

        l.SetLineColor(4)
        l.SetLineWidth(2)
        l.Draw()

        latex = TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.04)
        latex.SetTextAlign(12)  # 12 center, 31align right
        latex.DrawLatex(0.2, 0.85, signal)
        latex.DrawLatex(0.2, 0.75, 'Max Cut ' + maxcut)
        latex.DrawLatex(0.2, 0.65, 'Base Value = ' + str(BC) + ' Max Value  ' + str(maxvalue))
        latex.SetTextAlign(31)  # align right
        latex.DrawLatex(0.25, 0.93, " 36 fb^{-1} at #sqrt{s} = 13 TeV")
        latex.DrawLatex(0.8, 0.93, "CMS Internal")

        # legend.SetShadowColor(0);
        # legend.SetFillColor(0);
        # legend.SetLineColor(0);
        # legend.Draw('same')
        savename = str(dirname + "/Selected/" + signal + str(nhist[k]))
        c1.SaveAs(savename + ".png")
        c1.SaveAs(savename + ".pdf")
        # c1.SaveAs(savename + ".C")
        # c1.SaveAs(savename + ".root")
        c1.Modified()
        c1.Update()

        # for n in range(len(signame)):
        ftable = dirname + "/Selected/" + signal + "_selected_table"
        with open(ftable, 'w+') as f:
            f.write("cutname total_b total_b_err total_s total_s_err  cut_array_used\n")
            for j in range(len(xname)):
                # print signal
                print('%15s   %5s    %5s  %5s  %5s  %5s' % (
                    str(xname[j]), str(bkg[j]), str(bkg_err[j]), str(sig[j]), str(sig_err[j]), sig_cut_for_plots[j]))
                # f.write('%15s   %5s    %5s  %5s  %5s  %5s' %(str(xname[j]),  str(bkg[j]), str(bkg_err[j]) , str(sig[j]) , str(sig_err[j]) , sig_cut_for_plots[j]))
                f.write('%15s   %5f    %5f  %5f  %5f  %5s\n' % (
                    xname[j], bkg[j], bkg_err[j], sig[j], sig_err[j], sig_cut_for_plots[j]))
        f.close()
    print "Total number of cut combinations used for plotting: " + str(len(xname_all))


# Make plots of the log files in the selected folder.

def make_selected_plots(sigstr_BCuts, samples, dirname, vararray, signal, selarray=None, hist_config=None,
                        legend_config=None, extratag=""):
    """
    Plot histograms and dump logs for just the cuts decided on by `makesob_plots`.

    :param sigstr_BCuts: list of background cut names
    :param samples: `SampleManager` containing the (full-dimensional) simulated data
    :param dirname: directory to save logs and histograms to
    :param vararray: variables of interest (dependent variables) to aggregate counts and error estimates for
    :param signal: name of signal (independent variable) to plot
    :param selarray: list of cuts (see `makeselection` for format)
    :param hist_config: dict of `ROOT` plotting options
    :param legend_config: legend plotting options
    :param extratag: log file ending string (optional)
    """
    seldir = dirname + "/Selected/"
    fcuts = seldir + signal + "_selected_table"
    with open(fcuts, 'r') as f:
        lines = f.readlines()
        n = 0
        for line in lines:
            # print line
            n = n + 1
            listline = line.split()
            # print n
            if n > 1 and line.strip():  # first line is title and if line is not empty
                # print line
                # print listline[0]
                sigcut1 = str(line.split('ph_n==', 1)[1])  # this selects the part after first occurance of 'ph'
                sigcut = str(sigcut1.split('*', 1)[0])  # this selects the part after first occurance of '*'
                sigcut = "(ph_n==" + sigcut
                cutname = signal + listline[0]
                # pair = sigcut, listline[0]
                # print sigcut
                # listline = line.split()
                # print sigcut
                selarray = [([sigcut, cutname],)]
                makeplots(1, samples, vararray, seldir, selarray, hist_config, {}, extratag)

    f.close()


def doratio(h1, h2):
    # print "Inbside do ratio()"
    hratio = h1.Clone("hratio")
    hratio.Divide(h2)
    # print hratio.GetBinContent(1)
    '''
    hratio.SetMarkerStyle(20)
    hratio.SetMarkerSize(1.1)
    hratio.SetStats(0)
    hratio.SetTitle("")
    hratio.GetYaxis().SetTitle("ratio")
    hratio.SetLineColor(ROOT.kBlack)
    hratio.SetLineWidth(2)
    hratio.GetYaxis().SetTitleSize(0.10)
    hratio.GetYaxis().SetTitleOffset(0.6)
    hratio.GetYaxis().SetLabelSize(0.10)
    hratio.GetXaxis().SetLabelSize(0.10)
    hratio.GetXaxis().SetTitleSize(0.10)
    hratio.GetXaxis().SetTitleOffset(1.0)
    hratio.GetYaxis().SetRangeUser(0.5,1.5)
    #hratio.GetYaxis().UnZoom()
    hratio.GetYaxis().CenterTitle()
    hratio.GetYaxis().SetNdivisions(506, True)
    '''
    return hratio


def cmsinternal(pad):
    pad.cd()
    tex = ROOT.TLatex(0.18, 0.93, "CMS Internal")
    tex.SetNDC()
    tex.SetTextSize(0.05)
    tex.SetLineWidth(2)
    tex.Draw()
    return tex


def ratioline(hratio):
    left_edge = hratio.GetXaxis().GetXmin()
    right_edge = hratio.GetXaxis().GetXmax()

    oneline = ROOT.TLine(left_edge, 1, right_edge, 1)
    oneline.SetLineStyle(3)
    oneline.SetLineWidth(2)
    oneline.SetLineColor(ROOT.kBlack)
    oneline.Draw()
    return oneline


# Formatting for histograms.
def hformat(h1, color):
    h1.SetLineColor(color)
    h1.SetMarkerColor(color)
    # h1.SetMarkerStyle(20)
    h1.SetMarkerSize(1.1)
    h1.SetStats(0)
    h1.SetLineWidth(2)
    h1.GetYaxis().SetTitleSize(0.05)
    h1.GetYaxis().SetTitleOffset(1.15)
    h1.GetYaxis().SetLabelSize(0.05)
    h1.GetXaxis().SetLabelSize(0.05)
    h1.GetXaxis().SetTitleSize(0.05)
    h1.GetXaxis().SetTitleOffset(0.8)
