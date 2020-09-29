import numpy as np
import uproot
from root_numpy import fill_hist
from ROOT import TH2F, TCanvas, TLegend, TGaxis, gStyle, kRainBow, kRed, kAzure # pylint: disable=import-error,no-name-in-module

gStyle.SetOptStat(0)
gStyle.SetPadLeftMargin(0.12)
gStyle.SetPadBottomMargin(0.12)
gStyle.SetPadRightMargin(0.12)
gStyle.SetPadTopMargin(0.05)
gStyle.SetTitleSize(0.05, 'xy')
gStyle.SetTitleOffset(0.9, 'y')
gStyle.SetLabelSize(0.045, 'xy')
gStyle.SetPalette(kRainBow)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
TGaxis.SetMaxDigits(3)

kineTree = uproot.open('Bevents.root')['treeEvents']
kineDf = kineTree.pandas.df()
dfQuarks = kineDf.query('pdg == 5')
dfAntiQuarks = kineDf.query('pdg == -5')
# do not consider higher mass states
dfHadrons = kineDf.query('pdg in [511, 521, 531, 5122]')
dfAntiHadrons = kineDf.query('pdg in [-511, -521, -531, -5122]')
dfPart = dfQuarks.merge(dfAntiHadrons, left_on='event_number', right_on='event_number', suffixes=('Q', 'H'))
dfAntiPart = dfAntiQuarks.merge(dfHadrons, left_on='event_number', right_on='event_number', suffixes=('Q', 'H'))

hRapidityQuarkVsHadron = TH2F('hRapidityQuarkVsHadron', ';#it{y}(b);#it{y}(H_{b})', 400, -8., 8., 400, -8., 8.)
hRapidityAntiQuarkVsHadron = TH2F('hRapidityAntiQuarkVsHadron', ';#it{y}(b);#it{y}(H_{b})', 400, -8., 8., 400, -8., 8.)
fill_hist(hRapidityQuarkVsHadron, np.array([dfPart['yQ'].values, dfPart['yH'].values]).transpose())
fill_hist(hRapidityAntiQuarkVsHadron, np.array([dfAntiPart['yQ'].values, dfAntiPart['yH'].values]).transpose())
hRapidityQuarkVsHadron.Add(hRapidityAntiQuarkVsHadron)
hRapidityQuark = hRapidityQuarkVsHadron.ProjectionX()
hRapidityQuark.SetLineColor(kAzure+4)
hRapidityQuark.SetLineWidth(2)
hRapidityHadron = hRapidityQuarkVsHadron.ProjectionY()
hRapidityHadron.SetLineColor(kRed+1)
hRapidityHadron.SetLineWidth(2)

yBinMin = hRapidityHadron.GetXaxis().FindBin(-0.5)
yBinMax = hRapidityHadron.GetXaxis().FindBin(0.5)

cCorr = TCanvas('cCorr', '', 500, 500)
cCorr.SetLogz()
hRapidityQuarkVsHadron.Draw('colz')

leg = TLegend(0.6, 0.75, 0.8, 0.9)
leg.SetTextSize(0.045)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.AddEntry(hRapidityQuark, 'b-quark', 'l')
leg.AddEntry(hRapidityHadron, 'b-hadron', 'l')

cY = TCanvas('cY', '', 500, 500)
hRapidityQuark.GetXaxis().SetTitle('#it{y}')
hRapidityQuark.GetYaxis().SetTitle('d#sigma/d#it{y} (a.u)')
hRapidityQuark.Draw()
hRapidityHadron.Draw('same')
leg.Draw()
print(hRapidityHadron.Integral(yBinMin, yBinMax)/hRapidityQuark.Integral(yBinMin, yBinMax))

cCorr.SaveAs('Beauty_Quark_vs_Hadron_rapidity.pdf')
cY.SaveAs('Beauty_Quark_Hadron_rapidity_distribution.pdf')

input('Press enter to exit')
