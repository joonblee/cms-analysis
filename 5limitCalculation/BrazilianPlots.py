import ROOT
from ROOT import TFile, TTree, TCanvas, TGraph, TMultiGraph, TGraphErrors, TLegend
import CMS_lumi, tdrstyle
import subprocess # to execute shell command
ROOT.gROOT.SetBatch(ROOT.kTRUE)

# CMS style
CMS_lumi.cmsText = "CMS"
CMS_lumi.extraText = "Preliminary"
CMS_lumi.cmsTextSize = 0.65
CMS_lumi.outOfFrame = True
tdrstyle.setTDRStyle()


# CREATE datacards
def createDataCardsThetaB(labels,values):

    datacard_lines = [ "# automatic generated counting experiment",
                       "imax 1  number of channels",
                       "jmax 8  number of backgrounds",
                       "kmax 10  number of nuisance parameters (sources of systematical uncertainties)",
                       "------------",
                       "bin bin1"
                      ]

		### systematic uncertainties ###
    datacard_lines_unc = [ "------------",
                           "lumi       lnN   1.11     1.11   1.11   1.11   1.11   1.11   1.11   1.11   1.11  ",
                           "xs_sig     lnN   1.16      -      -      -      -      -      -      -      -    ",
                           "xs_higgs   lnN    -       1.50    -      -      -      -      -      -      -    ",
                           "xs_ggZZ    lnN    -        -     1.50    -      -      -      -      -      -    ",
                           "xs_qqZZ    lnN    -        -      -     1.50    -      -      -      -      -    ",
                           "xs_DY      lnN    -        -      -      -     1.50    -      -      -      -    ",
                           "xs_ZG      lnN    -        -      -      -      -     1.50    -      -      -    ",
                           "xs_RB_I    lnN    -        -      -      -      -      -     1.50    -      -    ",
                           "xs_RB_II   lnN    -        -      -      -      -      -      -     1.50    -    ",
                           "bg_others  lnN    -        -      -      -      -      -      -      -     1.30  "
                          ]

    # make datacards for differents values of theta_B
    for label, theta_B in zip(labels,values):

        rates = {"data":0, "sig":0, "higgs":0, "gg->ZZ":0, "qq->ZZ":0, "DY":0, "ZG":0, "RB_I":0, "RB_II":0, "Others":0}

				### You should fill here as data/bkg yields & sig rate ###
        if theta_B ==  200:
          rates = {"data":        81, "sig":    15.519, "higgs":         0, "gg->ZZ":    6.0873, "qq->ZZ":    44.807, "DY":         0, "ZG":    7.4084, "RB_I":    2.0518, "RB_II":    6.7751, "Others":   0.71378}
        elif theta_B ==  210:
          rates = {"data":        77, "sig":    15.696, "higgs":         0, "gg->ZZ":    6.0727, "qq->ZZ":    45.577, "DY":    2.5622, "ZG":    6.1011, "RB_I":    2.1452, "RB_II":    6.3472, "Others":   0.72177}
        elif theta_B ==  220:
          rates = {"data":        76, "sig":    15.863, "higgs": 0.0004352, "gg->ZZ":    5.6943, "qq->ZZ":    43.218, "DY":    2.5622, "ZG":    3.9221, "RB_I":    2.1338, "RB_II":    6.9812, "Others":   0.73797}
        elif theta_B ==  230:
          rates = {"data":        82, "sig":    16.021, "higgs": 0.0004352, "gg->ZZ":    5.2887, "qq->ZZ":     40.69, "DY":    5.1244, "ZG":     2.179, "RB_I":    2.8089, "RB_II":    7.0729, "Others":   0.77778}
        elif theta_B ==  240:
          rates = {"data":        78, "sig":     16.17, "higgs": 0.0004352, "gg->ZZ":    4.7847, "qq->ZZ":    37.539, "DY":    5.1244, "ZG":    3.9221, "RB_I":    2.1783, "RB_II":     6.401, "Others":   0.79954}
        elif theta_B ==  250:
          rates = {"data":        72, "sig":    16.309, "higgs":         0, "gg->ZZ":    4.3171, "qq->ZZ":    34.044, "DY":    5.1244, "ZG":    3.0505, "RB_I":   0.82207, "RB_II":    4.7365, "Others":   0.82217}
        elif theta_B ==  260:
          rates = {"data":        55, "sig":    16.439, "higgs":         0, "gg->ZZ":    3.8671, "qq->ZZ":    31.469, "DY":         0, "ZG":    3.4863, "RB_I":     0.826, "RB_II":    3.7358, "Others":   0.85065}
        elif theta_B ==  270:
          rates = {"data":        55, "sig":    16.561, "higgs":         0, "gg->ZZ":    3.4623, "qq->ZZ":      28.9, "DY":    2.5622, "ZG":     2.179, "RB_I":   0.67479, "RB_II":    3.0927, "Others":   0.79014}
        elif theta_B ==  280:
          rates = {"data":        50, "sig":    16.674, "higgs":         0, "gg->ZZ":    3.1199, "qq->ZZ":    26.358, "DY":    2.5622, "ZG":    1.7432, "RB_I":   0.82254, "RB_II":    3.3057, "Others":   0.82193}
        elif theta_B ==  290:
          rates = {"data":        49, "sig":    16.778, "higgs":         0, "gg->ZZ":    2.7768, "qq->ZZ":    24.715, "DY":    2.5622, "ZG":    3.4863, "RB_I":   0.29602, "RB_II":    2.8283, "Others":   0.93513}
        elif theta_B ==  300:
          rates = {"data":        36, "sig":    16.875, "higgs":         0, "gg->ZZ":    2.4683, "qq->ZZ":    22.498, "DY":    2.5622, "ZG":    4.3579, "RB_I":   0.55753, "RB_II":    2.9411, "Others":   0.90994}
        elif theta_B ==  310:
          rates = {"data":        35, "sig":    16.963, "higgs":         0, "gg->ZZ":    2.2084, "qq->ZZ":    21.261, "DY":    2.5622, "ZG":    4.3579, "RB_I":   0.49662, "RB_II":    2.3151, "Others":   0.87832}
        elif theta_B ==  320:
          rates = {"data":        32, "sig":    17.044, "higgs":         0, "gg->ZZ":    1.9811, "qq->ZZ":    19.767, "DY":    2.5622, "ZG":    3.4863, "RB_I":   0.82503, "RB_II":    1.9308, "Others":   0.84386}
        elif theta_B ==  330:
          rates = {"data":        32, "sig":    17.117, "higgs":         0, "gg->ZZ":    1.7827, "qq->ZZ":    18.009, "DY":         0, "ZG":    2.6147, "RB_I":   0.97498, "RB_II":    1.8252, "Others":   0.75641}
        elif theta_B ==  340:
          rates = {"data":        31, "sig":    17.183, "higgs":         0, "gg->ZZ":    1.5911, "qq->ZZ":    17.259, "DY":    2.5622, "ZG":    1.7432, "RB_I":    0.9255, "RB_II":    1.2687, "Others":   0.78713}
        elif theta_B ==  350:
          rates = {"data":        31, "sig":    17.241, "higgs":         0, "gg->ZZ":    1.4339, "qq->ZZ":    16.015, "DY":    2.5622, "ZG":    1.3074, "RB_I":   0.66336, "RB_II":    1.0438, "Others":   0.75213}
        elif theta_B ==  360:
          rates = {"data":        27, "sig":    17.292, "higgs":         0, "gg->ZZ":    1.2948, "qq->ZZ":    14.805, "DY":    2.5622, "ZG":    1.3074, "RB_I":   0.63604, "RB_II":     1.063, "Others":   0.69257}
        elif theta_B ==  370:
          rates = {"data":        24, "sig":    17.337, "higgs":         0, "gg->ZZ":    1.1814, "qq->ZZ":    13.885, "DY":    2.5622, "ZG":    1.7432, "RB_I":   0.40778, "RB_II":   0.99992, "Others":   0.71836}
        elif theta_B ==  380:
          rates = {"data":        22, "sig":    17.374, "higgs":         0, "gg->ZZ":    1.0927, "qq->ZZ":    13.196, "DY":         0, "ZG":     2.179, "RB_I":   0.40736, "RB_II":   0.87996, "Others":   0.70521}
        elif theta_B ==  390:
          rates = {"data":        22, "sig":    17.405, "higgs":         0, "gg->ZZ":    1.0135, "qq->ZZ":    12.243, "DY":         0, "ZG":    3.0505, "RB_I":   0.27497, "RB_II":   0.96372, "Others":   0.70112}
        elif theta_B ==  400:
          rates = {"data":        19, "sig":     17.43, "higgs":         0, "gg->ZZ":   0.93697, "qq->ZZ":    11.425, "DY":         0, "ZG":    3.0505, "RB_I":   0.34097, "RB_II":   0.96305, "Others":   0.64299}
        elif theta_B ==  420:
          rates = {"data":        20, "sig":     17.46, "higgs":         0, "gg->ZZ":   0.80459, "qq->ZZ":    10.039, "DY":         0, "ZG":    1.7432, "RB_I":   0.21626, "RB_II":   0.84617, "Others":   0.62215}
        elif theta_B ==  440:
          rates = {"data":        19, "sig":    17.467, "higgs":         0, "gg->ZZ":   0.68615, "qq->ZZ":    8.8153, "DY":         0, "ZG":   0.87158, "RB_I":   0.31365, "RB_II":   0.73992, "Others":   0.54924}
        elif theta_B ==  460:
          rates = {"data":        14, "sig":    17.451, "higgs":         0, "gg->ZZ":   0.60012, "qq->ZZ":    8.0243, "DY":         0, "ZG":         0, "RB_I":   0.48223, "RB_II":   0.51121, "Others":   0.49013}
        elif theta_B ==  480:
          rates = {"data":         8, "sig":    17.415, "higgs":         0, "gg->ZZ":   0.53049, "qq->ZZ":     7.301, "DY":         0, "ZG":         0, "RB_I":   0.42163, "RB_II":   0.37929, "Others":   0.45877}
        elif theta_B ==  500:
          rates = {"data":         7, "sig":    17.358, "higgs":         0, "gg->ZZ":   0.46619, "qq->ZZ":    6.3951, "DY":    2.5622, "ZG":         0, "RB_I":    0.4266, "RB_II":    0.3482, "Others":   0.45314}
        elif theta_B ==  520:
          rates = {"data":         5, "sig":    17.283, "higgs":         0, "gg->ZZ":   0.40677, "qq->ZZ":    5.6988, "DY":    2.5622, "ZG":         0, "RB_I":   0.23735, "RB_II":   0.41148, "Others":   0.46774}
        elif theta_B ==  540:
          rates = {"data":         6, "sig":     17.19, "higgs":         0, "gg->ZZ":   0.36209, "qq->ZZ":    5.0093, "DY":    2.5622, "ZG":         0, "RB_I":   0.10233, "RB_II":   0.39146, "Others":   0.47683}
        elif theta_B ==  560:
          rates = {"data":         6, "sig":    17.081, "higgs":         0, "gg->ZZ":   0.31861, "qq->ZZ":    4.5361, "DY":    2.5622, "ZG":         0, "RB_I":  0.053893, "RB_II":   0.34669, "Others":   0.47879}
        elif theta_B ==  580:
          rates = {"data":         5, "sig":    16.956, "higgs":         0, "gg->ZZ":   0.29039, "qq->ZZ":    4.1643, "DY":    2.5622, "ZG":   0.43579, "RB_I":         0, "RB_II":    0.3781, "Others":    0.4252}
        elif theta_B ==  600:
          rates = {"data":         5, "sig":    16.818, "higgs":         0, "gg->ZZ":   0.26458, "qq->ZZ":     3.833, "DY":         0, "ZG":   0.43579, "RB_I":  0.024008, "RB_II":   0.21958, "Others":   0.37517}
        elif theta_B ==  640:
          rates = {"data":         4, "sig":    16.502, "higgs":         0, "gg->ZZ":   0.21073, "qq->ZZ":    3.1435, "DY":         0, "ZG":         0, "RB_I":   0.13188, "RB_II":   0.12322, "Others":   0.27755}
        elif theta_B ==  680:
          rates = {"data":         3, "sig":    16.142, "higgs":         0, "gg->ZZ":   0.17431, "qq->ZZ":     2.677, "DY":         0, "ZG":         0, "RB_I":   0.18177, "RB_II":  0.072672, "Others":   0.23957}
        elif theta_B ==  720:
          rates = {"data":         5, "sig":    15.744, "higgs":         0, "gg->ZZ":   0.14551, "qq->ZZ":    2.3323, "DY":         0, "ZG":         0, "RB_I":   0.11996, "RB_II":  0.055687, "Others":   0.20927}
        elif theta_B ==  760:
          rates = {"data":         3, "sig":    15.317, "higgs":         0, "gg->ZZ":   0.12349, "qq->ZZ":    2.0686, "DY":         0, "ZG":         0, "RB_I": 0.0075518, "RB_II":  0.046843, "Others":   0.20015}
        elif theta_B ==  800:
          rates = {"data":         3, "sig":    14.867, "higgs":         0, "gg->ZZ":   0.10439, "qq->ZZ":    1.7374, "DY":         0, "ZG":         0, "RB_I":         0, "RB_II":  0.033523, "Others":   0.15521}
        elif theta_B ==  840:
          rates = {"data":         2, "sig":    14.401, "higgs":         0, "gg->ZZ":  0.087468, "qq->ZZ":    1.5413, "DY":         0, "ZG":         0, "RB_I":   0.06887, "RB_II":  0.025917, "Others":   0.10418}
        elif theta_B ==  880:
          rates = {"data":         0, "sig":    13.923, "higgs":         0, "gg->ZZ":  0.075997, "qq->ZZ":     1.325, "DY":         0, "ZG":         0, "RB_I":   0.14584, "RB_II":  0.026494, "Others":   0.11019}
        elif theta_B ==  920:
          rates = {"data":         0, "sig":    13.439, "higgs":         0, "gg->ZZ":  0.065157, "qq->ZZ":    1.0952, "DY":         0, "ZG":         0, "RB_I":   0.15049, "RB_II":  0.026579, "Others":  0.085466}
        elif theta_B ==  960:
          rates = {"data":         0, "sig":    12.955, "higgs":         0, "gg->ZZ":  0.057184, "qq->ZZ":   0.97347, "DY":         0, "ZG":         0, "RB_I":   0.16386, "RB_II":  0.024526, "Others":   0.08748}
        elif theta_B == 1000:
          rates = {"data":         0, "sig":    12.474, "higgs":         0, "gg->ZZ":  0.051219, "qq->ZZ":    0.8315, "DY":         0, "ZG":         0, "RB_I":  0.083638, "RB_II": 0.0087842, "Others":   0.10013}
        elif theta_B == 1100:
          rates = {"data":         0, "sig":    11.312, "higgs":         0, "gg->ZZ":  0.036651, "qq->ZZ":   0.54082, "DY":         0, "ZG":         0, "RB_I":  0.013294, "RB_II": 0.0059252, "Others":  0.045544}
        elif theta_B == 1200:
          rates = {"data":         0, "sig":    10.244, "higgs":         0, "gg->ZZ":  0.024721, "qq->ZZ":   0.35829, "DY":         0, "ZG":         0, "RB_I":  0.027184, "RB_II":  0.003423, "Others":  0.030879}
        elif theta_B == 1300:
          rates = {"data":         0, "sig":    9.3007, "higgs":         0, "gg->ZZ":  0.018411, "qq->ZZ":   0.31773, "DY":         0, "ZG":         0, "RB_I":   0.03562, "RB_II": 0.0029125, "Others":  0.012652}
        elif theta_B == 1400:
          rates = {"data":         0, "sig":    8.4949, "higgs":         0, "gg->ZZ":  0.014798, "qq->ZZ":   0.35829, "DY":         0, "ZG":         0, "RB_I":         0, "RB_II":0.00034144, "Others":  0.008095}
        elif theta_B == 1500:
          rates = {"data":         0, "sig":    7.8212, "higgs":         0, "gg->ZZ":  0.011471, "qq->ZZ":   0.27717, "DY":         0, "ZG":         0, "RB_I":         0, "RB_II":0.00034144, "Others":  0.008095}
        elif theta_B == 1600:
          rates = {"data":         0, "sig":    7.2561, "higgs":         0, "gg->ZZ": 0.0094064, "qq->ZZ":   0.25689, "DY":         0, "ZG":         0, "RB_I":         0, "RB_II":0.00073388, "Others":  0.008095}
        elif theta_B == 1700:
          rates = {"data":         0, "sig":    6.7581, "higgs":         0, "gg->ZZ": 0.0073416, "qq->ZZ":   0.18929, "DY":         0, "ZG":         0, "RB_I":         0, "RB_II":0.00039244, "Others":  0.008095}
        elif theta_B == 1800:
          rates = {"data":         0, "sig":    6.2678, "higgs":         0, "gg->ZZ": 0.0055636, "qq->ZZ":   0.12844, "DY":         0, "ZG":         0, "RB_I":         0, "RB_II":0.00039244, "Others":  0.008095}
        elif theta_B == 1900:
          rates = {"data":         0, "sig":    5.7075, "higgs":         0, "gg->ZZ": 0.0040149, "qq->ZZ":  0.081122, "DY":         0, "ZG":         0, "RB_I":         0, "RB_II":0.00039244, "Others":  0.008095}
        elif theta_B == 2000:
          rates = {"data":         0, "sig":    4.9816, "higgs":         0, "gg->ZZ": 0.0028678, "qq->ZZ":  0.067602, "DY":         0, "ZG":         0, "RB_I":         0, "RB_II":0.00039244, "Others":  0.008095}




        # lumi = 35900 -> (pb) | 35.9 -> (fb)

        datacard = open("datacard_"+label+".txt", 'w')
        for line in datacard_lines:
          datacard.write(line+"\n")
        datacard.write("observation %(data)s \n" % rates)
        datacard.write("------------\n")
        datacard.write("bin             bin1      bin1      bin1       bin1     bin1    bin1     bin1     bin1      bin1      \n")
        datacard.write("process         sig      higgs      ggZZ       qqZZ      DY      ZG     RB_I     RB_II     others   \n")
        datacard.write("process          0         1         2          3        4       5        6        7         8      \n")
        datacard.write("rate          %(sig)s  %(higgs)s %(gg->ZZ)s %(qq->ZZ)s %(DY)s  %(ZG)s %(RB_I)s %(RB_II)s %(Others)s \n" % rates)
        for line in datacard_lines_unc:
          datacard.write(line+"\n")
        datacard.close()
        print ">>>   datacard_"+label+".txt created"

    return labels


# EXECUTE datacards
def executeDataCards(labels):

    for label in labels:
        file_name = "datacard_"+label+".txt"
        combine_command = "combine -M Asymptotic -m 125 -n %s %s" % (label,file_name)
        print ""
        print ">>> " + combine_command
        p = subprocess.Popen(combine_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            print line.rstrip("\n")
        print ">>>   higgsCombine"+label+".Asymptotic.mH125.root created"
        retval = p.wait()


# GET limits from root file
def getLimits(file_name):

    file = TFile(file_name)
    tree = file.Get("limit")

    limits = [ ]
    for quantile in tree:
        limits.append(tree.limit)
        print ">>>   %.2f" % limits[-1]

    return limits[:6]


# PLOT upper limits
def plotUpperLimits(labels,values):
    # see CMS plot guidelines: https://ghm.web.cern.ch/ghm/plots/
    
    N = len(labels)
    yellow = TGraph(2*N)    # yellow band
    green = TGraph(2*N)     # green band
    median = TGraph(N)      # median line

    obs = TGraph(N)

    up2s = [ ]
    up_obs = []
    for i in range(N):
        file_name = "higgsCombine"+labels[i]+".Asymptotic.mH125.root"
        limit = getLimits(file_name)
        up2s.append(limit[4])
        yellow.SetPoint(    i,    values[i], limit[4] ) # + 2 sigma
        green.SetPoint(     i,    values[i], limit[3] ) # + 1 sigma
        median.SetPoint(    i,    values[i], limit[2] ) # median
        green.SetPoint(  2*N-1-i, values[i], limit[1] ) # - 1 sigma
        yellow.SetPoint( 2*N-1-i, values[i], limit[0] ) # - 2 sigma

        up2s.append(limit[5])
        obs.SetPoint(i, values[i], limit[5] )

    W = 1600 #800
    H = 1200 #600
    T = 0.08*H
    B = 0.12*H
    L = 0.12*W
    R = 0.04*W
    c = TCanvas("c","c",100,100,W,H)
    #c = TCanvas("c","c",100,100,W,H)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin( L/W )
    c.SetRightMargin( R/W )
    c.SetTopMargin( T/H )
    c.SetBottomMargin( B/H )
    c.SetTickx(0)
    c.SetTicky(0)
    c.SetGrid()
    c.cd()
    frame = c.DrawFrame(1.4,0.001, 4.1, 10)
    frame.GetYaxis().CenterTitle()
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetXaxis().SetNdivisions(508)
    frame.GetYaxis().CenterTitle(True)
    frame.GetYaxis().SetTitle("95% upper limit on #sigma (fb)")
    #frame.GetYaxis().SetTitle("95% upper limit on #sigma / #sigma_{SM}")
#    frame.GetYaxis().SetTitle("95% upper limit on #sigma #times BR / (#sigma #times BR)_{SM}")
    frame.GetXaxis().SetTitle("M(4mu) [GeV]")
    frame.SetMinimum(0)

    frame.SetMaximum(max(up2s)*1.05)
    frame.GetXaxis().SetLimits(min(values),max(values))

    yellow.SetFillColor(ROOT.kOrange)
    yellow.SetLineColor(ROOT.kOrange)
    yellow.SetFillStyle(1001)
    yellow.Draw('F')
    
    green.SetFillColor(ROOT.kGreen+1)
    green.SetLineColor(ROOT.kGreen+1)
    green.SetFillStyle(1001)
    green.Draw('Fsame')

    median.SetLineColor(1)
    median.SetLineWidth(2)
    median.SetLineStyle(2)
    median.Draw('Lsame')

    obs.SetLineColor(1)
    obs.SetLineWidth(2)
    obs.SetLineStyle(7)
    obs.Draw('Lsame')


    CMS_lumi.CMS_lumi(c,13,11)
    ROOT.gPad.SetTicks(1,1)
    frame.Draw('sameaxis')

    x1 = 0.7
    x2 = x1 + 0.24
    y1 = 0.70
    y2 = y1 + 0.16
    legend = TLegend(x1,y1,x2,y2)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(0.041)
    legend.SetTextFont(42)
    legend.AddEntry(obs, "observed",'L')
    legend.AddEntry(median, "Asymptotic CL_{s} expected",'L')
    legend.AddEntry(green, "#pm 1 std. deviation",'f')
#    legend.AddEntry(green, "Asymptotic CL_{s} #pm 1 std. deviation",'f')
    legend.AddEntry(yellow,"#pm 2 std. deviation",'f')
#    legend.AddEntry(green, "Asymptotic CL_{s} #pm 2 std. deviation",'f')
    legend.Draw()

    print " "
    c.SaveAs("UpperLimit.png")
    c.Close()


# RANGE of floats
def frange(start, stop, step):
    i = start
    while i <= stop:
        yield i
        i += step


# MAIN
def main():

    labels = [ ]
    values = [ ]

    ll = [200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,420,440,460,480,500,520,540,560,580,600,640,680,720,760,800,840,880,920,960,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000]
    #ll = [200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,360,370,380,390,400,420,440,460,480,500]
    for theta_B in ll:
        values.append(theta_B)
        label = "%d" % (theta_B)
        labels.append(label)

    createDataCardsThetaB(labels,values)
    executeDataCards(labels)
    plotUpperLimits(labels,values)



if __name__ == '__main__':
    main()
    
