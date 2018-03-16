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
from NoiseExtraction import NoiseExtraction
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


class PedestalAnalysis:
	def __init__(self, run=22011, dir='', force=False):
		print 'Creating PedestalAnalysis instance for run:', run
		self.run = run
		self.dir = dir + '/{r}'.format(r=self.run)
		self.force = force
		self.bar = None
		self.rootFile = self.pedTree = None
		self.adc_vect = self.ped_vect = self.sigma_vect = self.cm_vect = self.ped_cmc_vect = self.sigma_cmc_vect = None
		self.adc_hist, self.ped_hist, self.sigma_hist, self.cm_hist, self.ped_cmc_hist, self.sigma_cmc_hist = ro.TH3D(), ro.TH3D(), ro.TH3D(), ro.TH2D(), ro.TH3D(), ro.TH3D()
		self.dicBraNames = {'rawTree.DiaADC': 'adc_' + str(self.run), 'diaPedestalMean': 'ped_' + str(self.run), 'diaPedestaSigma': 'sigma_' + str(self.run),
		                    'diaPedestalMeanCMN': 'ped_cmc_' + str(self.run), 'diaPedestaSigmaCMN': 'sigma_cmc_' + str(self.run), 'commonModeNoise': 'cm_' + str(self.run)}
		self.allBranches = self.dicBraNames.keys()
		self.listBraNames1ch = ['commonModeNoise']
		self.listBraNamesChs = [x for x in self.allBranches if x not in self.listBraNames1ch]
		self.dicBraVectChs = {}
		self.dicBraVect1ch = {}
		self.entries = 0
		self.dicBraHist = {}
		self.ev_axis = {'min': 0, 'max': 0, 'bins': 1000}
		self.ch_axis = {'min': -0.5, 'max': diaChs - 0.5, 'bins': diaChs}
		self.dicBraProf = {}
		self.tcutg_gr = None
		self.tcutg_dia = None
		self.tcutg_chs = {}
		self.gr, self.low_ch, self.up_ch = 0, 0, 0
		self.gr_plots, self.dia_channels_plots, self.nc_channels_plots = {}, {}, {}
		self.event_plots = {}
		self.event_list = None
		self.dia_channel_list = None
		self.nc_channel_list = None
		self.dic_bra_mean_ped_chs = {}
		self.dic_bra_sigma_ped_chs = {}
		self.cm_adc = None

		if self.force:
			self.ClearOldAnalysis()

		self.dicHasHistos = {bra: False for bra in self.allBranches}
		self.hasHistos = self.CheckHistograms()

		self.dicHasProfiles = {bra: False for bra in self.listBraNamesChs}
		self.hasProfiles = self.CheckProfiles()

		if not self.hasProfiles:
			if not self.hasHistos:
				self.LoadROOTFile()
				self.LoadVectorsFromBranches()
				self.CreateHistograms()
				self.FillHistograms()
				self.SaveHistograms()
				self.CreateProfiles()
				self.SaveProfiles()
			else:
				self.LoadHistograms()
				self.CreateProfiles()
				self.SaveProfiles()
		else:
			self.LoadHistograms()
			self.LoadProfiles()
			self.GetMeanSigmaPerChannel()

	def ClearOldAnalysis(self):
		if os.path.isdir('{d}/pedestalAnalysis/histos'.format(d=self.dir)):
			shutil.rmtree('{d}/pedestalAnalysis/histos'.format(d=self.dir))
		if os.path.isdir('{d}/pedestalAnalysis/profiles'.format(d=self.dir)):
			shutil.rmtree('{d}/pedestalAnalysis/profiles'.format(d=self.dir))

	def SetTCutG(self, ev_ini=0, ev_end=0, color=ro.kRed + 3):
		ev_fin = self.entries if ev_end == 0 else ev_end

		self.tcutg_gr = ro.TCutG('TCutG_gr', 5)
		self.tcutg_gr.SetVarX('xaxis')
		self.tcutg_gr.SetVarY('yaxis')
		self.tcutg_gr.SetPoint(0, ev_ini, self.gr - 0.5)
		self.tcutg_gr.SetPoint(1, ev_ini, self.gr + 0.5)
		self.tcutg_gr.SetPoint(2, ev_fin, self.gr + 0.5)
		self.tcutg_gr.SetPoint(3, ev_fin, self.gr - 0.5)
		self.tcutg_gr.SetPoint(4, ev_ini, self.gr - 0.5)
		self.tcutg_gr.SetLineWidth(3)
		self.tcutg_gr.SetLineStyle(7)
		self.tcutg_gr.SetLineColor(color)

		self.tcutg_dia = ro.TCutG('TCutG_dia', 5)
		self.tcutg_dia.SetVarX('xaxis')
		self.tcutg_dia.SetVarY('yaxis')
		self.tcutg_dia.SetPoint(0, ev_ini, self.low_ch - 0.5)
		self.tcutg_dia.SetPoint(1, ev_ini, self.up_ch + 0.5)
		self.tcutg_dia.SetPoint(2, ev_fin, self.up_ch + 0.5)
		self.tcutg_dia.SetPoint(3, ev_fin, self.low_ch - 0.5)
		self.tcutg_dia.SetPoint(4, ev_ini, self.low_ch - 0.5)
		self.tcutg_dia.SetLineWidth(3)
		self.tcutg_dia.SetLineStyle(7)
		self.tcutg_dia.SetLineColor(color)

		for ch in xrange(self.low_ch, self.up_ch + 1):
			self.tcutg_chs[ch] = ro.TCutG('TCutG_ch_'+str(ch), 5)
			self.tcutg_chs[ch].SetVarX('xaxis')
			self.tcutg_chs[ch].SetVarY('yaxis')
			self.tcutg_chs[ch].SetPoint(0, ev_ini, ch - 0.5)
			self.tcutg_chs[ch].SetPoint(1, ev_ini, ch + 0.5)
			self.tcutg_chs[ch].SetPoint(2, ev_fin, ch + 0.5)
			self.tcutg_chs[ch].SetPoint(3, ev_fin, ch - 0.5)
			self.tcutg_chs[ch].SetPoint(4, ev_ini, ch - 0.5)
			self.tcutg_chs[ch].SetLineWidth(3)
			self.tcutg_chs[ch].SetLineStyle(7)
			self.tcutg_chs[ch].SetLineColor(color)

	def SetChannels(self):
		cont = False
		while not cont:
			temp = raw_input('Type the channel for the GR (should be between 0 and 127): ')
			if IsInt(temp):
				if 0 <= int(temp) < diaChs:
					self.gr = int(temp)
					cont = True
		cont = False
		while not cont:
			temp = raw_input('Type the lower channel to analyse (should be between 0 and 127): ')
			if IsInt(temp):
				if 0 <= int(temp) < diaChs:
					self.low_ch = int(temp)
					cont = True
		cont = False
		while not cont:
			temp = raw_input('Type the upper channel to analyse (should be between 0 and 127): ')
			if IsInt(temp):
				if 0 <= int(temp) < diaChs:
					self.up_ch = int(temp)
					cont = True
		self.dia_channel_list = range(self.low_ch, self.up_ch + 1)
		self.nc_channel_list = [ch for ch in xrange(diaChs) if ch not in self.dia_channel_list]
		self.SetTCutG()

	def GetChannelsHistos(self):
		for branch, prof in self.dicBraProf.iteritems():
			name = self.dicBraNames[branch]
			self.gr_plots[branch] = self.ProjectChannel(self.gr, prof, name)
		for ch in self.dia_channel_list:
			# self.dia_channels_plots[ch] = {}
			self.dia_channels_plots[ch] = {branch: self.ProjectChannel(ch, prof, self.dicBraNames[branch]) for branch, prof in self.dicBraProf.iteritems()}
			# for branch, prof in self.dicBraProf.iteritems():
			# 	name = self.dicBraNames[branch]
			# 	self.dia_channels_plots[ch][branch] = self.ProjectChannel(ch, prof, name)
		for ch in self.nc_channel_list:
			# self.nc_channels_plots[ch] = {}
			self.nc_channels_plots[ch] = {branch: self.ProjectChannel(ch, prof, self.dicBraNames[branch]) for branch, prof in self.dicBraProf.iteritems()}
			# for branch, prof in self.dicBraProf.iteritems():
			# 	name = self.dicBraNames[branch]
			# 	self.nc_channels_plots[ch][branch] = self.ProjectChannel(ch, prof, name)

	def GetEventsHistos(self):
		for branch, prof in self.dicBraProf.iteritems():
			name = self.dicBraNames[branch]
			for ev_bin in self.event_list:
				self.event_plots[branch] = self.ProjectEvent(ev_bin, prof, name)
		for branch in self.listBraNames1ch:
			name = self.dicBraNames[branch]
			hist = self.dicBraHist[branch]
			for ev_bin in self.event_list:
				self.event_plots[branch] = self.ProjectEvent(ev_bin, hist, name)

	def SetEventsList(self, bins=11, bin_low=0, bin_up=0):
		bin_l = 1 if bin_low == 0 else bin_low
		bin_u = self.ev_axis['bins'] if bin_up == 0 else bin_up
		n_bins = 1 if bin_l == bin_u else bins
		self.event_list = np.array(np.round(np.linspace(bin_l, bin_u, n_bins)), dtype='uint')

	def GetPlots(self):
		if self.gr == self.low_ch == self.up_ch:
			self.SetChannels()
		if len(self.gr_plots) == 0 or len(self.dia_channels_plots) == 0:
			self.GetChannelsHistos()
		if not self.event_list:
			self.SetEventsList()
		if len(self.event_plots) == 0:
			self.GetEventsHistos()

	def ProjectChannel(self, ch, prof, name):
		temp = prof.ProjectionX('' + name[:-6] + '_ch_' + str(ch), ch + 1, ch + 1)
		temp.SetTitle('' + name[:-6] + '_ch_' + str(ch))
		temp.GetXaxis().SetTitle('event')
		temp.GetYaxis().SetTitle(name[:-6])
		temp.SetFillColor(fillColor)
		temp.GetYaxis().SetRangeUser(prof.GetMinimum(), prof.GetMaximum())
		return temp

	def ProjectEvent(self, ev_bin, prof, name):
		event = int(round(prof.GetXaxis().GetBinCenter(int(ev_bin))))
		temp = prof.ProjectionY('' + name[:-6] + '_ev_' + str(event), int(ev_bin), int(ev_bin))
		temp.SetTitle('' + name[:-6] + '_ev_' + str(event))
		temp.GetXaxis().SetTitle('channel')
		temp.GetYaxis().SetTitle(name[:-6])
		temp.SetFillColor(fillColor)
		temp.GetYaxis().SetRangeUser(prof.GetMinimum(), prof.GetMaximum())
		return temp

	def GetMeanSigmaPerChannel(self):
		for branch, prof in self.dicBraProf.iteritems():
			# ipdb.set_trace()
			temp1 = prof.ProfileY(prof.GetTitle() + '_py')
			self.dic_bra_mean_ped_chs[branch] = np.array([temp1.GetBinContent(ch + 1) for ch in xrange(diaChs + 1)], dtype='f8')
			self.dic_bra_sigma_ped_chs[branch] = np.array([temp1.GetBinError(ch + 1) for ch in xrange(diaChs + 1)], dtype='f8')
			del temp1

	def SetProfileDefaults(self):
		for branch, prof in self.dicBraProf.iteritems():
			prof.GetXaxis().SetTitle('event')
			prof.GetYaxis().SetTitle('channel')
			prof.GetZaxis().SetTitle(self.dicBraNames[branch][:-6])
			prof.GetZaxis().SetRangeUser(self.dicBraHist[branch].GetZaxis().GetXmin(), self.dicBraHist[branch].GetZaxis().GetXmax())
			prof.SetStats(False)

	def SetHistogramDefaults(self):
		for branch, hist in self.dicBraHist.iteritems():
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
					self.dicHasProfiles[branch] = False
				else:
					self.dicHasProfiles[branch] = True
			return np.array(self.dicHasProfiles.values(), '?').all()

	def CheckHistograms(self):
		if not os.path.isdir('{d}/pedestalAnalysis/histos'.format(d=self.dir)):
			print 'Pedestal analysis directory "histos" does not exist. All the vectors for the analysis will be created'
			return False
		else:
			for branch in self.allBranches:
				name = self.dicBraNames[branch]
				if not os.path.isfile('{d}/pedestalAnalysis/histos/{n}.root'.format(d=self.dir, n=name)):
					self.dicHasHistos[branch] = False
				else:
					self.dicHasHistos[branch] = True
			return np.array(self.dicHasHistos.values(), '?').all()

	def LoadROOTFile(self):
		print 'Loading ROOT file...',
		sys.stdout.flush()
		# ipdb.set_trace()
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
			self.dicBraHist[branch] = deepcopy(temp.Get(self.dicBraNames[branch]))
			temp.Close()
			del temp
		self.entries = int(self.dicBraHist[self.listBraNamesChs[0]].GetXaxis().GetXmax())
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
		self.dicBraProf = {branch: histo.Project3DProfile('yx') for branch, histo in self.dicBraHist.iteritems() if branch in self.listBraNamesChs}
		# for branch, histo in self.dicBraHist.iteritems():
		# 	if branch in self.listBraNamesChs:
		# 		self.dicBraProf[branch] = histo.Project3DProfile('yx')
		self.SetProfileDefaults()

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

	def LoadVectorFromBranch(self, branch):
		channels = self.pedTree.GetLeaf(branch).GetLen()
		leng = self.pedTree.Draw(branch, '', 'goff', self.entries, 0)
		if leng == -1:
			print 'Error, could not load the branches. try again'
			return
		while leng > self.pedTree.GetEstimate():
			self.pedTree.SetEstimate(leng)
			leng = self.pedTree.Draw(branch, '', 'goff', self.entries, 0)
		temp = self.pedTree.GetVal(0)
		if branch in self.listBraNamesChs:
			self.dicBraVectChs[branch] = np.array([[temp[ev * channels + ch] for ch in xrange(channels)] for ev in xrange(self.entries)], dtype='f8')
		elif branch in self.listBraNames1ch:
			self.dicBraVect1ch[branch] = np.array([temp[ev] for ev in xrange(self.entries)], dtype='f8')
		del temp

	def SetHistogramLimits(self):
		for branch, vect in self.dicBraVect1ch.iteritems():
			name = self.dicBraNames[branch]
			if name.startswith('cm'):
				ymin, ymax = min(cm_axis['min'], vect.min()), max(cm_axis['max'], vect.max())
				cm_axis['min'], cm_axis['max'] = ymin, ymax
		for branch, vect in self.dicBraVectChs.iteritems():
			name = self.dicBraNames[branch]
			if name.startswith('ped'):
				zmin, zmax = min(ped_axis['min'], vect.min()), max(ped_axis['max'], vect.max())
				ped_axis['min'], ped_axis['max'] = zmin, zmax
			elif name.startswith('sigm'):
				zmin, zmax = min(sigma_axis['min'], vect.min()), max(sigma_axis['max'], vect.max())
				sigma_axis['min'], sigma_axis['max'] = zmin, zmax
			elif name.startswith('adc'):
				zmin, zmax = min(adc_axis['min'], vect.min()), max(adc_axis['max'], vect.max())
				adc_axis['min'], adc_axis['max'] = zmin, zmax

	def CreateHistograms(self):
		print 'Creating histograms...',
		sys.stdout.flush()
		self.SetHistogramLimits()
		for branch in self.allBranches:
			name = self.dicBraNames[branch]
			if branch in self.listBraNames1ch:
				ymin, ymax, ybins = 0, 0, 0
				if name.startswith('cm'):
					ymin, ymax, ybins = cm_axis['min'], cm_axis['max'], int((cm_axis['max'] - cm_axis['min']) / 0.5)
				self.dicBraHist[branch] = ro.TH2D(name, name, self.ev_axis['bins'], self.ev_axis['min'], self.ev_axis['max'], ybins + 1, ymin, ymax + (ymax - ymin) / float(ybins))
				self.dicBraHist[branch].GetXaxis().SetTitle('event')
				self.dicBraHist[branch].GetYaxis().SetTitle(name[:-6])
			elif branch in self.listBraNamesChs:
				zmin, zmax, zbins = 0, 0, 0
				if name.startswith('ped'):
					zmin, zmax, zbins = ped_axis['min'], ped_axis['max'], int(min((ped_axis['max'] - ped_axis['min']) / 10.0, 500))
				elif name.startswith('adc'):
					zmin, zmax, zbins = adc_axis['min'], adc_axis['max'], int(min((adc_axis['max'] - adc_axis['min']) / 10.0, 500))
				elif name.startswith('sigma'):
					zmin, zmax, zbins = sigma_axis['min'], sigma_axis['max'], int(min((sigma_axis['max'] - sigma_axis['min']) / 0.2, 500))
				self.dicBraHist[branch] = ro.TH3D(name, name, self.ev_axis['bins'], self.ev_axis['min'], self.ev_axis['max'], self.ch_axis['bins'], self.ch_axis['min'], self.ch_axis['max'], zbins + 1,
				                                  zmin, zmax + (zmax - zmin) / float(zbins))
				self.dicBraHist[branch].GetXaxis().SetTitle('event')
				self.dicBraHist[branch].GetYaxis().SetTitle('channel')
				self.dicBraHist[branch].GetZaxis().SetTitle(name[:-6])
		self.adc_hist, self.ped_hist, self.sigma_hist, self.ped_cmc_hist, self.sigma_cmc_hist, self.cm_hist = self.dicBraHist['rawTree.DiaADC'], self.dicBraHist['diaPedestalMean'], self.dicBraHist[
			'diaPedestaSigma'], self.dicBraHist['diaPedestalMeanCMN'], self.dicBraHist['diaPedestaSigmaCMN'], self.dicBraHist['commonModeNoise']
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

	def Calculate_CM_ADC(self):
		if len(self.dicBraVectChs) == 0:
			self.LoadROOTFile()
			print 'Loading ADCs...', ; sys.stdout.flush()
			self.LoadVectorFromBranch('rawTree.DiaADC')
			self.adc_vect = self.dicBraVectChs['rawTree.DiaADC']
			print 'Done'
		if len(self.dicBraVect1ch) == 0:
			self.LoadROOTFile()
			print 'Loading CMs...', ; sys.stdout.flush()
			self.LoadVectorFromBranch('commonModeNoise')
			self.cm_vect = self.dicBraVect1ch['commonModeNoise']
			print 'Done'
		temp = []
		self.CreateProgressBar(self.entries)
		if self.bar is not None:
			self.bar.start()
		for ev in xrange(self.entries):
			temp.append(self.adc_vect[ev].mean())
		self.cm_adc = np.array(temp, 'f8')
		del temp

	def SaveHistograms(self):
		print 'Saving histograms:'
		if not os.path.isdir('{d}/pedestalAnalysis/histos'.format(d=self.dir)):
			os.makedirs('{d}/pedestalAnalysis/histos'.format(d=self.dir))
		for branch, histo in self.dicBraHist.iteritems():
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
	parser.add_option('-r', '--run', dest='run', default=22022, type='int', help='Run to be analysed (e.g. 22022)')
	parser.add_option('-d', '--dir', dest='dir', default='.', type='string', help='source folder containing processed data of different runs')
	parser.add_option('-f', '--force', dest='force', default=False, action='store_true')

	(options, args) = parser.parse_args()
	run = int(options.run)
	dir = str(options.dir)
	force = bool(options.force)
	# output = str(options.output)
	# connect = int(options.connect)
	# low = int(options.low)
	# high = int(options.high)

	pedAna = PedestalAnalysis(run, dir, force)
