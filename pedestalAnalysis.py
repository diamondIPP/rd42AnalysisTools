#!/usr/bin/env python
# from ROOT import TFile, TH2F, TH3F, TH1F, TCanvas, TCutG, kRed, gStyle, TBrowser, Long, TF1
from optparse import OptionParser
# from numpy import array, floor, average, std
import numpy as np
import ROOT as ro
import ipdb  # set_trace, launch_ipdb_on_exception
import progressbar
from copy import deepcopy
from NoiseExtraction import NoiseExtraction
import os, sys
from Utils import *

__author__ = 'DA'

dicTypes = {'Char_t': 'int8', 'UChar_t': 'uint8', 'Short_t': 'short', 'UShort_t': 'ushort', 'Int_t': 'int32', 'UInt_t': 'uint32', 'Float_t': 'float32', 'Double_t': 'float64', 'Long64_t': 'int64',
            'ULong64_t': 'uint64', 'Bool_t': 'bool'}
diaChs = 128


class PedestalAnalysis:
	def __init__(self, run=22011, dir=''):
		print 'Creating NoiseExtraction instance for run:', run
		self.run = run
		self.dir = dir
		self.bar = None
		self.rootFile = self.pedTree = None
		self.adc_vect = self.ped_vect = self.sigma_vect = self.cm_vect = self.ped_cmc_vect = self.sigma_cmc_vect = None
		self.adc_hist, self.ped_hist, self.sigma_hist, self.cm_hist, self.ped_cmc_hist, self.sigma_cmc_hist = ro.TH3D(), ro.TH3D(), ro.TH3D(), ro.TH2D(), ro.TH3D(), ro.TH3D()
		self.dicBraNames = {'rawTree.DiaADC': 'adc_' + str(self.run), 'diaPedestalMean': 'ped_' + str(self.run), 'diaPedestaSigma': 'sigma_' + str(self.run),
		                    'diaPedestalMeanCMN': 'ped_cmc_' + str(self.run), 'diaPedestaSigmaCMN': 'sig_cmc_' + str(self.run), 'commonModeNoise': 'cm_' + str(self.run)}
		self.allBranches = self.dicBraNames.keys()
		self.listBraNames1ch = ['commonModeNoise']
		self.listBraNamesChs = [x for x in self.allBranches if x not in self.listBraNames1ch]
		self.dicBraVectChs = {}
		self.dicBraVect1ch = {}
		self.entries = 0
		self.dictBraHist = {}
		self.ev_axis = {'min': 0, 'max': 0, 'bins': 1000}
		self.ch_axis = {'min': -0.5, 'max': diaChs - 0.5, 'bins': diaChs}
		self.dicBraProf = {}
		self.tcutg = None

		self.dictHasHistos = {bra: False for bra in self.allBranches}
		self.hasHistos = self.CheckHistograms()

		self.dictHasProfiles = {bra: False for bra in self.listBraNamesChs}
		self.hasProfiles = self.CheckProfiles()

		if not self.hasProfiles:
			if not self.hasHistos:
				self.LoadROOTFile()
				self.LoadVectorsFromBranches()
				self.CreateHistogramsForVectors()
				self.FillHistograms()
				self.SaveHistograms()
				self.CreateProgressBar()
				self.SaveProfiles()
			else:
				self.LoadHistograms()
				self.CreateProfiles()
				self.SaveProfiles()
		else:
			self.LoadHistograms()
			self.LoadProfiles()

	def SetTCutG(self, name='def', ch_ini=82, ch_end=82, ev_ini=0, ev_end=0, color=ro.kRed + 3):
		self.tcutg = ro.TCutG(name, 5)
		self.tcutg.SetVarX('xaxis')
		self.tcutg.SetVarY('yaxis')
		ev_fin = self.entries if ev_end == 0 else ev_end
		self.tcutg.SetPoint(0, ev_ini, ch_ini - 0.5)
		self.tcutg.SetPoint(1, ev_ini, ch_end + 0.5)
		self.tcutg.SetPoint(2, ev_fin, ch_end + 0.5)
		self.tcutg.SetPoint(3, ev_fin, ch_ini - 0.5)
		self.tcutg.SetPoint(4, ev_ini, ch_ini - 0.5)
		self.tcutg.SetLineWidth(3)
		self.tcutg.SetLineStyle(7)
		self.tcutg.SetLineColor(color)

	def GetPlots(self):
		cont = False
		while not cont:
			temp = raw_input('Type the channel for the GR (should be between 0 and 127): ')
			if IsInt(temp):
				if 0 <= int(temp) < diaChs:
					gr = int(temp)
					cont = True
		cont = False
		while not cont:
			temp = raw_input('Type the lower channel to analyse (should be between 0 and 127): ')
			if IsInt(temp):
				if 0 <= int(temp) < diaChs:
					lower = int(temp)
					cont = True
		cont = False
		while not cont:
			temp = raw_input('Type the upper channel to analyse (should be between 0 and 127): ')
			if IsInt(temp):
				if 0 <= int(temp) < diaChs:
					upper = int(temp)
					cont = True

		self.SetTCutG('gr_fid', gr, gr, 0, self.entries)


	def ProjectChannel(self, ch):
		for branch, prof in self.dicBraProf.iteritems():
			name = self.dicBraNames[branch]
			temp = prof.ProjectionX(''+name[:-6]+'_ch_'+str(ch), ch, ch)
			

	def SetProfileDefaults(self):
		for branch, prof in self.dicBraProf.iteritems():
			prof.GetXaxis().SetTitle('event')
			prof.GetYaxis().SetTitle('channel')
			prof.GetZaxis().SetTitle(self.dicBraNames[branch][:-6])

	def SetHistogramDefaults(self):
		for branch, hist in self.dictBraHist.iteritems():
			name = self.dicBraNames[branch]
			hist.GetXaxis().SetTitle('event')
			hist.GetYaxis().SetTitle('channel') if branch in self.listBraNamesChs else hist.GetYaxis().SetTitle(name[:-6])
			if branch in self.listBraNamesChs:
				hist.GetZaxis().SetTitle(name[:-6])

	def CheckProfiles(self):
		if not os.path.isdir('{d}/pedestalAnalysis/profiles'.format(d=self.dir)):
			print 'Pedestal analysis directory "profiles" does not exist.'
			return False
		else:
			for branch in self.listBraNamesChs:
				name = self.dicBraNames[branch]
				if not os.path.isfile('{d}/pedestalAnalysis/profiles/{n}_pyx.root'.format(d=self.dir, n=name)):
					self.dictHasProfiles[branch] = False
				else:
					self.dictHasProfiles[branch] = True
			return np.array(self.dictHasProfiles.values(), '?').all()

	def CheckHistograms(self):
		if not os.path.isdir('{d}/pedestalAnalysis/histos'.format(d=self.dir)):
			print 'Pedestal analysis directory "histos" does not exist. All the vectors for the analysis will be created'
			return False
		else:
			for branch in self.allBranches:
				name = self.dicBraNames[branch]
				if not os.path.isfile('{d}/pedestalAnalysis/histos/{n}.root'.format(d=self.dir, n=name)):
					self.dictHasHistos[branch] = False
				else:
					self.dictHasHistos[branch] = True
			return np.array(self.dictHasHistos.values(), '?').all()

	def LoadROOTFile(self):
		print 'Loading ROOT file...',
		sys.stdout.flush()
		self.rootFile = ro.TFile('{d}/pedestalData.{r}.root'.format(d=self.dir, r=self.run), 'READ')
		self.pedTree = self.rootFile.Get('pedestalTree')
		self.entries = self.pedTree.GetEntries()
		# self.entries = int(maxEntries)
		print 'Done'

	def LoadHistograms(self):
		print 'Loading Histograms...',
		sys.stdout.flush()
		for branch in self.allBranches:
			temp = ro.TFile('{d}/pedestalAnalysis/histos/{n}.root'.format(d=self.dir, n=self.dicBraNames[branch]), 'read')
			self.dictBraHist[branch] = deepcopy(temp.Get(self.dicBraNames[branch]))
			temp.Close()
			del temp
		self.entries = int(self.dictBraHist[self.listBraNamesChs[0]].GetXaxis().GetXmax())
		self.SetHistogramDefaults()
		print 'Done'

	def LoadProfiles(self):
		print 'Loading Profiles...',
		sys.stdout.flush()
		for branch in self.listBraNamesChs:
			temp = ro.TFile('{d}/pedestalAnalysis/profiles/{n}_pyx.root'.format(d=self.dir, n=self.dicBraNames[branch]), 'read')
			self.dicBraProf[branch] = deepcopy(temp.Get('{n}_pyx'.format(n=self.dicBraNames[branch])))
			temp.Close()
			del temp
		self.entries = int(self.dicBraProf[self.listBraNamesChs[0]].GetXaxis().GetXmax())
		self.SetProfileDefaults()
		print 'Done'

	def CreateProfiles(self):
		print 'Creating profiles:'
		for branch, histo in self.dictBraHist.iteritems():
			if branch in self.listBraNamesChs:
				self.dicBraProf[branch] = histo.Project3DProfile('yx')

	def LoadVectorsFromBranches(self, first_ev=0):
		print 'Loading vectors from branches...',
		sys.stdout.flush()
		if self.pedTree is None:
			self.LoadROOTFile()

		num_bra_chs = len(self.listBraNamesChs)
		if num_bra_chs < 1:
			print 'The dictionary of branches and vectors is empty! try again'
			return
		channels = self.pedTree.GetLeaf(self.listBraNamesChs[0]).GetLen()
		for branch in self.listBraNamesChs:
			if self.pedTree.GetLeaf(branch).GetLen() != channels:
				print 'The given branches have different sizes! try again'
				return
		leng = self.pedTree.Draw(':'.join(self.listBraNamesChs), '', 'goff para', self.entries, first_ev)
		if leng == -1:
			print 'Error, could not load the branches. try again'
			return
		while leng > self.pedTree.GetEstimate():
			self.pedTree.SetEstimate(leng)
			leng = self.pedTree.Draw(':'.join(self.listBraNamesChs), '', 'goff para', self.entries, first_ev)
		self.entries = leng / channels
		for pos, branch in enumerate(self.listBraNamesChs):
			temp = self.pedTree.GetVal(pos)
			self.dicBraVectChs[branch] = np.array([[temp[ev * channels + ch] for ch in xrange(channels)] for ev in xrange(self.entries)], dtype='f8')
			del temp

		num_bra_1ch = len(self.listBraNames1ch)
		if num_bra_1ch < 1:
			print 'The dictionary of branches and vectors is empty! try again'
			return
		channel = 1
		for branch in self.listBraNames1ch:
			if self.pedTree.GetLeaf(branch).GetLen() != channel:
				print 'The given branches have different sizes different to 1! try again'
				return
		leng = self.pedTree.Draw(':'.join(self.listBraNames1ch), '', 'goff para', self.entries, first_ev)
		if leng == -1:
			print 'Error, could not load the branches. try again'
			return
		while leng > self.pedTree.GetEstimate():
			self.pedTree.SetEstimate(leng)
			leng = self.pedTree.Draw(':'.join(self.listBraNames1ch), '', 'goff para', self.entries, first_ev)
		for pos, branch in enumerate(self.listBraNames1ch):
			temp = self.pedTree.GetVal(pos)
			self.dicBraVect1ch[branch] = np.array([temp[ev] for ev in xrange(self.entries)], dtype='f8')
			del temp

		self.adc_vect, self.ped_vect, self.sigma_vect, self.ped_cmc_vect, self.sigma_cmc_vect = self.dicBraVectChs['rawTree.DiaADC'], self.dicBraVectChs['diaPedestalMean'], self.dicBraVectChs[
			'diaPedestaSigma'], self.dicBraVectChs['diaPedestalMeanCMN'], self.dicBraVectChs['diaPedestaSigmaCMN']
		self.cm_vect = self.dicBraVect1ch['commonModeNoise']
		self.ev_axis['max'] = self.entries
		print 'Done'

	def CreateHistogramsForVectors(self):
		print 'Creating histograms...',
		sys.stdout.flush()
		for branch in self.allBranches:
			name = self.dicBraNames[branch]
			if branch in self.listBraNames1ch:
				ymin, ymax = self.dicBraVect1ch[branch].min(), self.dicBraVect1ch[branch].max()
				ybins = 500
				self.dictBraHist[branch] = ro.TH2D(name, name, self.ev_axis['bins'], self.ev_axis['min'], self.ev_axis['max'], ybins + 1, ymin, ymax + (ymax - ymin) / float(ybins))
				self.dictBraHist[branch].GetXaxis().SetTitle('event')
				self.dictBraHist[branch].GetYaxis().SetTitle(name[:-6])
			elif branch in self.listBraNamesChs:
				zmin, zmax = self.dicBraVectChs[branch].min(), self.dicBraVectChs[branch].max()
				zbins = 500
				self.dictBraHist[branch] = ro.TH3D(name, name, self.ev_axis['bins'], self.ev_axis['min'], self.ev_axis['max'], self.ch_axis['bins'], self.ch_axis['min'], self.ch_axis['max'], zbins + 1, zmin, zmax + (zmax - zmin) / float(zbins))
				self.dictBraHist[branch].GetXaxis().SetTitle('event')
				self.dictBraHist[branch].GetYaxis().SetTitle('channel')
				self.dictBraHist[branch].GetZaxis().SetTitle(name[:-6])
		self.adc_hist, self.ped_hist, self.sigma_hist, self.ped_cmc_hist, self.sigma_cmc_hist, self.cm_hist = self.dictBraHist['rawTree.DiaADC'], self.dictBraHist['diaPedestalMean'], self.dictBraHist['diaPedestaSigma'], self.dictBraHist['diaPedestalMeanCMN'], self.dictBraHist['diaPedestaSigmaCMN'], self.dictBraHist['commonModeNoise']
		print 'Done'

	def FillHistograms(self):
		print 'Filling histograms:'
		self.CreateProgressBar(self.entries)
		if self.bar is not None:
			self.bar.start()
		for ev in xrange(self.entries):
			for ch in xrange(diaChs):
				self.adc_hist.Fill(ev, ch, self.adc_vect[ev, ch])
				self.ped_hist.Fill(ev, ch, self.ped_vect[ev, ch])
				self.sigma_hist.Fill(ev, ch, self.sigma_vect[ev, ch])
				self.ped_cmc_hist.Fill(ev, ch, self.ped_cmc_vect[ev, ch])
				self.sigma_cmc_hist.Fill(ev, ch, self.sigma_cmc_vect[ev, ch])
			self.cm_hist.Fill(ev, self.cm_vect[ev])
			if self.bar is not None:
				self.bar.update(ev + 1)
		if self.bar is not None:
			self.bar.finish()

	def SaveHistograms(self):
		print 'Saving histograms:'
		if not os.path.isdir('{d}/pedestalAnalysis/histos'.format(d=self.dir)):
			os.makedirs('{d}/pedestalAnalysis/histos'.format(d=self.dir))
		for branch, histo in self.dictBraHist.iteritems():
			name = self.dicBraNames[branch]
			histo.SaveAs('{d}/pedestalAnalysis/histos/{n}.root'.format(d=self.dir, n=name))

	def SaveProfiles(self):
		print 'Saving profiles:'
		if not os.path.isdir('{d}/pedestalAnalysis/profiles'.format(d=self.dir)):
			os.makedirs('{d}/pedestalAnalysis/profiles'.format(d=self.dir))
		for branch, prof in self.dicBraProf.iteritems():
			name = self.dicBraNames[branch]
			prof.SaveAs('{d}/pedestalAnalysis/profiles/{n}_pyx.root'.format(d=self.dir, n=name))

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
	parser.add_option('-r', '--run', dest='run', default=22011, type='int', help='Run to be analysed (e.g. 22011)')
	parser.add_option('-d', '--dir', dest='dir', default='.', type='string', help='source folder containing processed data of different runs')
	# parser.add_option('-o', '--outputDir', dest='output', default='', type='string', help='output folder containing the analysed results')
	# parser.add_option('-c', '--connected', dest='connect', default=1, type='int', help='do connected channels (1) or disconnected (0). Other integer value you would have to specify the range using -l and -H')
	# parser.add_option('-l', '--lowChannel', dest='low', default=0, type='int', help='lower channel to analyse e.g. 1')
	# parser.add_option('-H', '--highChannel', dest='high', default=0, type='int', help='higher channel to analyse e.g. 2')

	(options, args) = parser.parse_args()
	run = int(options.run)
	dir = str(options.dir)
	# output = str(options.output)
	# connect = int(options.connect)
	# low = int(options.low)
	# high = int(options.high)

	z = PedestalAnalysis(run, dir)
