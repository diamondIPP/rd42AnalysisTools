#!/usr/bin/env python
# from ROOT import TFile, TH2F, TH3F, TH1F, TCanvas, TCutG, kRed, gStyle, TBrowser, Long, TF1
from optparse import OptionParser
# from numpy import array, floor, average, std
import numpy as np
import ROOT as ro
import ipdb  # set_trace, launch_ipdb_on_exception
import progressbar
from copy import deepcopy
from collections import OrderedDict
from pedestalAnalysis import PedestalAnalysis
import os, sys, shutil
from Utils import *

__author__ = 'DA'

dicTypes = {'Char_t': 'int8', 'UChar_t': 'uint8', 'Short_t': 'short', 'UShort_t': 'ushort', 'Int_t': 'int32', 'UInt_t': 'uint32', 'Float_t': 'float32', 'Double_t': 'float64', 'Long64_t': 'int64',
            'ULong64_t': 'uint64', 'Bool_t': 'bool'}

diaChs = 128

fillColor = ro.TColor.GetColor(125, 153, 209)
sigma_axis = {'min': 0, 'max': 35}
adc_axis = {'min': 0, 'max': 2 ** 12 - 1}
ped_axis = {'min': 0, 'max': 2 ** 12 - 1}
cm_axis = {'min': -100, 'max': 100}
color_index = 10000

class CollectionAnalysis:
	def __init__(self, runs, dir='', noise=False, force=False):
		print 'Creating CollectionAnalysis instance for runs:', runs
		self.runs = runs
		self.dir = dir
		self.noise = noise
		self.force = False
		if force:
			ask = raw_input('are you sure you want to re-analyse all the selected runs???')
			if ask.lower() == 'true' or ask.lower() == 'yes' or ask.lower() == 'sure':
				self.force = True
			else:
				print 'The analysis on the selected rusn will not be overwritten'
		self.bar = None
		self.runs_pedAna = None
		self.graphs_pedAna = None
		self.multi_graph_pedAna = None
		self.canvas = None
		self.trash = []
		if self.noise:
			self.CreatePedestalAnalysisInstances()
			self.GetGraphsPedAnaChannels()

	def CreatePedestalAnalysisInstances(self):
		self.runs_pedAna = OrderedDict({runi: PedestalAnalysis(runi, self.dir, force) for runi in self.runs})

	def GetGraphsPedAnaChannels(self):
		branches = self.runs_pedAna[self.runs[0]].listBraNamesChs
		dicBraNames = {bra: name[:-6] for bra, name in self.runs_pedAna[self.runs[0]].dicBraNames.iteritems()}
		self.graphs_pedAna = {}
		for branch in branches:
			self.graphs_pedAna[branch] = OrderedDict()
			for ch in xrange(diaChs + 1):
				self.graphs_pedAna[branch][ch] = ro.TGraphErrors(len(self.runs), np.array(self.runs, dtype='f8'), np.array([runi.dic_bra_mean_ped_chs[branch][ch] for runi in self.runs_pedAna.itervalues()], dtype='f8'), np.zeros(len(self.runs), dtype='f8'), np.array([runi.dic_bra_sigma_ped_chs[branch][ch] for runi in self.runs_pedAna.itervalues()], dtype='f8'))
				self.graphs_pedAna[branch][ch].SetNameTitle(dicBraNames[branch] + '_ch_' + str(ch), dicBraNames[branch] + '_ch_' + str(ch))
				self.graphs_pedAna[branch][ch].GetXaxis().SetTitle('run')
				self.graphs_pedAna[branch][ch].GetYaxis().SetTitle(dicBraNames[branch])
				self.graphs_pedAna[branch][ch].SetMarkerStyle(7)
				self.graphs_pedAna[branch][ch].SetMarkerColor(ro.kBlack)

	def GetPedAnaBranchChannelsGraphs(self, branch='', chi=0, chf=diaChs-1, cmc=False):
		if 'ped' in branch.lower():
			if cmc or 'cmc' in branch.lower() or 'cm' in branch.lower():
				self.Get_Plot_Ped_CMC_Channels(chi, chf)
			else:
				self.Get_Plot_Ped_Channels(chi, chf)
		if 'sig' in branch.lower():
			if cmc or 'cmc' in branch.lower() or 'cm' in branch.lower():
				self.Get_Plot_Sigma_CMC_Channels(chi, chf)
			else:
				self.Get_Plot_Sigma_Channels(chi, chf)

	def Get_Plot_Ped_CMC_Channels(self, chi=0, chf=diaChs - 1):
		self.Create_multi_graph_channels('Ped_CMC_chs:[{i},{f}];run;ped_cmc'.format(i=chi, f=chf), self.graphs_pedAna['diaPedestalMeanCMN'], chi, chf)

	def Get_Plot_Ped_Channels(self, chi=0, chf=diaChs - 1):
		self.Create_multi_graph_channels('Ped_chs:[{i},{f}];run;ped'.format(i=chi, f=chf), self.graphs_pedAna['diaPedestalMean'], chi, chf)

	def Get_Plot_Sigma_CMC_Channels(self, chi=0, chf=diaChs - 1):
		self.Create_multi_graph_channels('Sigma_CMC_chs:[{i},{f}];run;sigma_cmc'.format(i=chi, f=chf), self.graphs_pedAna['diaPedestaSigmaCMN'], chi, chf)

	def Get_Plot_Sigma_Channels(self, chi=0, chf=diaChs - 1):
		self.Create_multi_graph_channels('Sigma_chs:[{i},{f}];run;sigma'.format(i=chi, f=chf), self.graphs_pedAna['diaPedestaSigma'], chi, chf)

	def Create_multi_graph_channels(self, title, graphs, chi, chf):
		if self.multi_graph_pedAna:
			self.multi_graph_pedAna.Delete()
			self.multi_graph_pedAna = None
			self.trash = []
		self.multi_graph_pedAna = ro.TMultiGraph()
		self.multi_graph_pedAna.SetTitle(title)
		graphs_copy = deepcopy(graphs)
		for ch in xrange(chi, chf + 1):
			r, g, b = ReturnRGB(ch, chi, chf)
			self.trash.append(ro.TColor(color_index+ch, r, g, b))
			graphs_copy[ch].SetLineColor(color_index+ch)
			graphs_copy[ch].SetMarkerColor(ro.kBlack)
			self.multi_graph_pedAna.Add(graphs_copy[ch])

	def PlotGraph(self, graph, options='ALP'):
		if self.canvas:
			self.canvas.Clear()
			self.canvas.Close()
			self.canvas.Destructor()
			self.canvas = None
		self.canvas = ro.TCanvas('ca0', 'ca0', 1)
		self.canvas.cd()
		graph.Draw(options)

	def CreateProgressBar(self, maxVal=1):
		widgets = [
			'Processed: ', progressbar.Counter(),
			' out of {mv} '.format(mv=maxVal), progressbar.Percentage(),
			' ', progressbar.Bar(marker='>'),
			' ', progressbar.Timer(),
			' ', progressbar.ETA()
			# ' ', progressbar.AdaptativeETA(),
			#  ' ', progressbar.AdaptativeTransferSpeed()
		]
		self.bar = progressbar.ProgressBar(widgets=widgets, maxval=maxVal)


if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-f', '--first', dest='first', default=22022, type='int', help='Run to be analysed (e.g. 22022)')
	parser.add_option('-l', '--last', dest='last', default=22022, type='int', help='Last run to be analysed e.g. 22023')
	parser.add_option('-d', '--dir', dest='dir', default='.', type='string', help='source folder containing processed data of different runs')
	parser.add_option('-n', '--noise', dest='noise', default=False, action='store_true')
	parser.add_option('-x', '--force', dest='force', default=False, action='store_true')
	# parser.add_option('-o', '--outputDir', dest='output', default='', type='string', help='output folder containing the analysed results')
	# parser.add_option('-c', '--connected', dest='connect', default=1, type='int', help='do connected channels (1) or disconnected (0). Other integer value you would have to specify the range using -l and -H')
	# parser.add_option('-l', '--lowChannel', dest='low', default=0, type='int', help='lower channel to analyse e.g. 1')
	# parser.add_option('-H', '--highChannel', dest='high', default=0, type='int', help='higher channel to analyse e.g. 2')

	(options, args) = parser.parse_args()
	first = int(options.first)
	last = int(options.last)
	dir = str(options.dir)
	force = bool(options.force)
	noise = bool(options.noise)

	if last >= first:
		runs = np.arange(first, last + 1)
		ca = CollectionAnalysis(runs, dir, noise, force)
	else:
		print 'Last run entered "{l}" is not greater or equal than first run "{f}". Try again'.format(l=last, f=first)
		exit()
