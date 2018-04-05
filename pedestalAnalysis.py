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
import cPickle as pickle

__author__ = 'DA'

dicTypes = {'Char_t': 'int8', 'UChar_t': 'uint8', 'Short_t': 'short', 'UShort_t': 'ushort', 'Int_t': 'int32', 'UInt_t': 'uint32', 'Float_t': 'float32', 'Double_t': 'float64', 'Long64_t': 'int64',
            'ULong64_t': 'uint64', 'Bool_t': 'bool'}

diaChs = 128

fillColor = ro.TColor.GetColor(125, 153, 209)
sigma_axis = {'min': 0, 'max': 35}
adc_axis = {'min': 0, 'max': 2**12 - 1}
ped_axis = {'min': 0, 'max': 2**12 - 1}
cm_axis = {'min': -100, 'max': 100}
ev_axis = {'min': 0, 'max': 0, 'bins': 1000}
ch_axis = {'min': -0.5, 'max': diaChs - 0.5, 'bins': diaChs}

class PedestalAnalysis:
	def __init__(self, run=22011, dir='', force=False):
		print 'Creating PedestalAnalysis instance for run:', run
		self.run = run
		self.bias = self.GetBiasUser()
		self.dia = self.GetDiamondNameUser()
		self.dir = dir + '/{r}'.format(r=self.run)
		self.force = force
		self.bar = None
		self.canvas = {}
		self.rootFile, self.pedTree = None, None
		self.adc_vect, self.ped_vect, self.sigma_vect, self.cm_vect, self.ped_cmc_vect, self.sigma_cmc_vect = None, None, None, None, None, None
		self.signal_vect, self.signal_cmc_vect = None, None
		self.biggest_adc_vect = None
		self.adc_ev_hist, self.ped_ev_hist, self.sigma_ev_hist, self.cm_ev_hist, self.ped_cmc_ev_hist, self.sigma_cmc_ev_hist = ro.TH3F(), ro.TH3F(), ro.TH3F(), ro.TH2F(), ro.TH3F(), ro.TH3F()
		self.signal_ev_hist, self.signal_cmc_ev_hist, self.biggest_adc_ev_hist = None, None, None
		self.signal_ch_hist, self.signal_cmc_ch_hist, self.biggest_adc_ch_hist = ro.TH2F(), ro.TH2F(), ro.TH2F()
		self.ped_ch_hist, self.ped_cmc_ch_hist, self.sigma_ch_hist, self.sigma_cmc_ch_hist, self.adc_ch_hist = ro.TH2F(), ro.TH2F(), ro.TH2F(), ro.TH2F(), ro.TH2F()
		self.dicBraNames = {'rawTree.DiaADC': 'adc', 'diaPedestalMean': 'ped', 'diaPedestaSigma': 'sigma', 'diaPedestalMeanCMN': 'ped_cmc', 'diaPedestaSigmaCMN': 'sigma_cmc', 'commonModeNoise': 'cm',
		                    'eventNumber': 'event', 'diaChannel': 'dia_ch'}
		self.allBranches = self.dicBraNames.keys()
		self.listBraNames1ch = ['commonModeNoise']
		self.listBraNamesChs = [x for x in self.allBranches if x not in self.listBraNames1ch]
		self.dicNewHistoNames = {'biggestADC': 'biggest_adc_' + str(self.run)}
		self.allNewHistos = self.dicNewHistoNames.keys()

		self.dicBraHist = {}
		self.dicBraProf = {}
		self.dicNewHist = {}

		if self.force:
			self.ClearOldAnalysis()

		self.dicHasNewHistos = {bra: False for bra in self.allNewHistos}
		self.hasNewHistos = self.CheckNewHistograms()

		self.LoadROOTFile()
		self.CheckChannelBranch()



		# TODO REDOING PEDESTAL ANALYSIS SUCH THAT IT UPDATES THE TREE WITH CHANNEL VECTOR AND THEN PLOT DIRECTLY FORM PYTHON WITHOUT SAVING VECTORS, ETC...

		# self.dicBraVectChs = {}
		# self.dicBraVect1ch = {}
		# self.dicNewVectChs = {}
		# self.dicNewVect1ch = {}
		# self.entries = 0
		# self.dicBraHist = {}
		# self.dicNewHist = {}
		# self.dicBraProf = {}
		# self.tcutg_gr = None
		# self.tcutg_dia = None
		# self.tcutg_chs = {}
		# self.gr, self.low_ch, self.up_ch = 0, 0, 0
		# self.gr_plots, self.dia_channels_plots, self.nc_channels_plots = {}, {}, {}
		# self.event_plots = {}
		# self.event_list = None
		# self.dia_channel_list = None
		# self.nc_channel_list = None
		# self.dic_bra_mean_ped_chs = {}
		# self.dic_bra_sigma_ped_chs = {}
		# self.cm_adc = None


		# self.dicHasVectors = {bra: False for bra in self.allBranches}
		# self.hasVectors = self.CheckVectorPickles()

		# self.dicHasHistos = {bra: False for bra in self.allBranches}

		self.GetMeanSigmaPerChannel()

	def GetBiasUser(self):
		bias = 'a'
		while not IsFloat(bias):
			bias = raw_input('Please enter the bias voltage in V used for this run ({r}) (e.g. -234.4): '.format(r=self.run))
		return float(bias)

	def GetDiamondNameUser(self):
		name = raw_input('Please enter the name of the diamond (e.g. II6-99): ')
		return name

	def ClearOldAnalysis(self):
		if os.path.isdir('{d}/pedestalAnalysis/histos'.format(d=self.dir)):
			shutil.rmtree('{d}/pedestalAnalysis/histos'.format(d=self.dir))
		if os.path.isdir('{d}/pedestalAnalysis/profiles'.format(d=self.dir)):
			shutil.rmtree('{d}/pedestalAnalysis/profiles'.format(d=self.dir))

	def CheckNewHistograms(self):
		if not os.path.isdir('{d}/pedestalAnalysis/histos'.format(d=self.dir)):
			print 'Pedestal analysis directory "histos" does not exist. All the histograms for the analysis will be created'
			return False
		else:
			for hist in self.allNewHistos:
				name = self.dicNewHistoNames[hist]
				if os.path.isfile('{d}/pedestalAnalysis/histos/{n}.root'.format(d=self.dir, n=name)):
					self.dicHasNewHistos[hist] = True
					tempf = ro.TFile('{d}/pedestalAnalysis/histos/{n}.root'.format(d=self.dir, n=name), 'READ')
					self.dicNewHist[hist] = deepcopy(tempf.Get(name))
					tempf.Close()
					del tempf
				else:
					self.dicHasNewHistos[hist] = False
			return np.array(self.dicHasNewHistos.values(), '?').all()

	def LoadROOTFile(self, mode='READ'):
		print 'Loading ROOT file...',
		sys.stdout.flush()
		# ipdb.set_trace()
		self.rootFile = ro.TFile('{d}/pedestalData.{r}.root'.format(d=self.dir, r=self.run), mode)
		self.pedTree = self.rootFile.Get('pedestalTree')
		self.entries = self.pedTree.GetEntries()
		# self.entries = int(maxEntries)
		print 'Done'

	def CheckChannelBranch(self):
		if not self.pedTree.GetLeaf('diaChannel'):
			print 'The tree does not have the branch "diaChannel". Creating branch:\n'
			self.rootFile.Close()
			self.LoadROOTFile(mode='UPDATE')
			chs = np.arange(diaChs, dtype='uint8')
			chsBra = self.pedTree.Branch('diaChannel', chs, 'diaChannel[{chs}]/b'.format(chs=diaChs))
			self.CreateProgressBar(self.entries)
			self.bar.start()
			for ev in xrange(self.entries):
				self.pedTree.GetEntry(ev)
				chsBra.Fill()
				self.bar.update(ev + 1)
			self.pedTree.Write()
			self.bar.finish()
			self.pedTree.Delete()
			self.rootFile.Close()
			self.LoadROOTFile(mode='READ')
			del chs, chsBra
			print 'Finished creating updating pedestalTree with branch "diaChannel"'
			return
		else:
			print 'The tree has the branch "diaChannel".\n'

	def Plot2DHisto(self, ybra='', xbra='', option='colz', name=''):
		if not self.pedTree.GetLeaf(xbra) or not self.pedTree.GetLeaf(ybra):
			print 'One of the two branches:', ybra, 'or', xbra, 'does not exist'
			return
		print 'Creating plot...', ; sys.stdout.flush()
		if not self.dicBraHist.has_key(ybra):
			self.dicBraHist[ybra] = {}
		if not self.dicBraHist[ybra].has_key(xbra):
			xbraN, ybraN = self.dicBraNames[xbra], self.dicBraNames[ybra]
			pol = 'Pos' if self.bias >= 0 else 'Neg'
			name0 = '{y}_vs_{x}_{d}_{r}_{p}_{v}V'.format(y=ybraN, x=xbraN, d=self.dia, r=self.run, p=pol, v=self.bias) if name == '' else name
			nameh = 'h_' + name0
			xmin, xmax, xbin = self.GetHistoLimits(xbra)
			ymin, ymax, ybin = self.GetHistoLimits(ybra)
			self.dicBraHist[ybra][xbra] = ro.TH2F(nameh, nameh, xbin, xmin, xmax, ybin, ymin, ymax)
			self.dicBraHist[ybra][xbra].GetXaxis().SetTitle(xbraN)
			self.dicBraHist[ybra][xbra].GetYaxis().SetTitle(ybraN)
			leng = self.pedTree.Draw('{y}:{x}>>{n}'.format(y=ybra, x=xbra, n=nameh), '', 'goff')
			while leng > self.pedTree.GetEstimate():
				self.pedTree.SetEstimate(leng)
				leng = self.pedTree.Draw('{y}:{x}>>{n}'.format(y=ybra, x=xbra, n=nameh), '', 'goff')
		name0 = self.dicBraHist[ybra][xbra].GetName()[2:]
		self.CreateCanvas(name=name0)
		self.dicBraHist[ybra][xbra].Draw(option)
		print 'Done'

	def GetHistoLimits(self, branch):
		if branch.lower().startswith('rawTree.DiaADC'.lower()) or branch.lower().startswith('diaPedestalM'.lower()):
			return adc_axis['min'], adc_axis['max'], int(min((adc_axis['max'] - adc_axis['min']) / 0.5, 10000))
		if branch.lower().startswith('diaPedestaS'.lower()):
			return sigma_axis['min'], sigma_axis['max'], int(min((sigma_axis['max'] - sigma_axis['min']) / 0.01, 10000))
		if branch.lower().startswith('common'.lower()):
			return cm_axis['min'], cm_axis['max'], int(min((cm_axis['max'] - cm_axis['min']) / 0.01, 10000))
		if branch.lower().startswith('diaChann'.lower()):
			return ch_axis['min'], ch_axis['max'], ch_axis['bins']
		if branch.lower().startswith('event'.lower()):
			return ev_axis['min'], ev_axis['max'], int(min((ev_axis['max'] - ev_axis['min']) / 100.0, 10000))

	def CreateCanvas(self, name='0'):
		self.CleanCanvasList()
		if not self.canvas.has_key(name):
			self.canvas[name] = ro.TCanvas('c_'+name, 'c_'+name, 1)
		self.canvas[name].cd()

	def CleanCanvasList(self):
		for key, value in self.canvas.items():
			if not value:
				self.canvas.pop(key)

	def Plot2DProfile(self, zbra='', ybra='', xbra='', option='colz', name=''):
		if not self.pedTree.GetLeaf(xbra) or not self.pedTree.GetLeaf(ybra) or not self.pedTree.GetLeaf(zbra):
			print 'One of the two branches:', zbra, 'or', ybra, 'or', xbra, 'does not exist'
			return
		print 'Creating plot...', ; sys.stdout.flush()
		if not self.dicBraProf.has_key(zbra):
			self.dicBraProf[zbra] = {}
		if not self.dicBraProf[zbra].has_key(ybra):
			self.dicBraProf[zbra][ybra] = {}
		if not self.dicBraProf[zbra][ybra].has_key(xbra):
			xbraN, ybraN, zbraN = self.dicBraNames[xbra], self.dicBraNames[ybra], self.dicBraNames[zbra]
			pol = 'Pos' if self.bias >= 0 else 'Neg'
			name0 = '{z}_{y}_vs_{x}_{d}_{r}_{p}_{v}V'.format(z=zbraN, y=ybraN, x=xbraN, d=self.dia, r=self.run, p=pol, v=self.bias) if name == '' else name
			namep = 'p_' + name0
			xmin, xmax, xbin = self.GetHistoLimits(xbra)
			ymin, ymax, ybin = self.GetHistoLimits(ybra)
			zmin, zmax, zbin = self.GetHistoLimits(zbra)
			self.dicBraProf[zbra][ybra][xbra] = ro.TProfile2D(namep, namep, xbin, xmin, xmax, ybin, ymin, ymax, zmin, zmax)
			self.dicBraProf[zbra][ybra][xbra].GetXaxis().SetTitle(xbraN)
			self.dicBraProf[zbra][ybra][xbra].GetYaxis().SetTitle(ybraN)
			self.dicBraProf[zbra][ybra][xbra].GetZaxis().SetTitle(ybraN)
			leng = self.pedTree.Draw('{z}:{y}:{x}>>{n}'.format(z=zbra, y=ybra, x=xbra, n=namep), '', 'prof goff')
			while leng > self.pedTree.GetEstimate():
				self.pedTree.SetEstimate(leng)
				leng = self.pedTree.Draw('{z}:{y}:{x}>>{n}'.format(z=zbra, y=ybra, x=xbra, n=namep), '', 'prof goff')
		name0 = self.dicBraProf[zbra][ybra][xbra].GetName()[2:]
		self.CreateCanvas(name=name0)
		self.dicBraProf[zbra][ybra][xbra].Draw(option)
		print 'Done'

	def Plot2DHistoBiggestHit(self, xbra='diaChannel', name=''):
		if not self.pedTree.GetLeaf(xbra):
			print 'The branch', xbra, 'does not exist'
			return
		if not self.dicNewHist.has_key('biggest'):
			self.dicNewHist['biggest'] = {}
		if not self.dicNewHist['biggest'].has_key(xbra):
			xbraN = self.dicBraNames[xbra]
			pol = 'Pos' if self.bias >= 0 else 'Neg'
			name0 = 'biggest_vs_{x}_{d}_{r}_{p}_{v}V'.format(x=xbraN, d=self.dia, r=self.run, p=pol, v=self.bias) if name == '' else name
			nameh = 'h_' + name0
			# xmin, xmax, xbin =

	def SaveObject(self, object, location='histos'):
		"""
		Saves a root object in a subfolder specified by location
		:param object: root object. e.g. TH2F, TProfile2D, TCanvas, etc...
		:param location: 'histos' for histograms, 'profiles' for TProfiles...
		:return: Nothing
		"""
		direct = '{d}/pedestalAnalysis/{l}'.format(d=self.dir, l=location)
		self.CheckDirectory(direct)
		name = object.GetName()
		object.SaveAs('{d}/{n}.root'.format(d=direct, n=name))

	def CheckDirectory(self, direct):
		if not os.path.isdir(direct):
			os.makedirs(direct)

	def AssignVectorsFromDictrionaries(self):
		self.adc_vect, self.ped_vect, self.sigma_vect, self.ped_cmc_vect, self.sigma_cmc_vect = self.dicBraVectChs['rawTree.DiaADC'], self.dicBraVectChs['diaPedestalMean'], self.dicBraVectChs[
			'diaPedestaSigma'], self.dicBraVectChs['diaPedestalMeanCMN'], self.dicBraVectChs['diaPedestaSigmaCMN']
		self.cm_vect = self.dicBraVect1ch['commonModeNoise']

	def LoadNewVectorsFromVectors(self):
		self.dicNewVectChs['signal'] = self.adc_vect - self.ped_vect
		self.dicNewVectChs['signalCMN'] = np.subtract(self.adc_vect - self.ped_vect, [[cmi] for cmi in self.cm_vect])
		self.signal_vect = self.dicNewVectChs['signal']
		self.signal_cmc_vect = self.dicNewVectChs['signalCMN']

	def SaveVectorsPickles(self):
		print 'Saving vectors...'
		if not os.path.isdir('{d}/pedestalAnalysis/vectors'.format(d=self.dir)):
			os.makedirs('{d}/pedestalAnalysis/vectors'.format(d=self.dir))
		for branch in self.allBranches:
			print 'Saving', branch, '...', ; sys.stdout.flush()
			name = self.dicBraNames[branch]
			if branch in self.listBraNamesChs:
				with open('{d}/pedestalAnalysis/vectors/{n}.dat'.format(d=self.dir, n=name), 'wb') as fs:
					pickle.dump(self.dicBraVectChs[branch], fs, pickle.HIGHEST_PROTOCOL)
			else:
				with open('{d}/pedestalAnalysis/vectors/{n}.dat'.format(d=self.dir, n=name), 'wb') as fs:
					pickle.dump(self.dicBraVect1ch[branch], fs, pickle.HIGHEST_PROTOCOL)
			print 'Done'

	def GetMeanSigmaPerChannel(self):
		for branch, prof in self.dicBraProf.iteritems():
			# ipdb.set_trace()
			temp1 = prof.ProfileY(prof.GetTitle() + '_py')
			self.dic_bra_mean_ped_chs[branch] = np.array([temp1.GetBinContent(ch + 1) for ch in xrange(diaChs + 1)], dtype='f8')
			self.dic_bra_sigma_ped_chs[branch] = np.array([temp1.GetBinError(ch + 1) for ch in xrange(diaChs + 1)], dtype='f8')
			del temp1

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

	def LoadVectorFromBranch(self, branch):
		if self.dicHasVectors[branch]:
			self.LoadVectorFromPickle(branch)
		else:
			self.LoadVectorFromTree(branch)

	def Calculate_CM_ADC(self):
		self.CheckExistanceVector('rawTree.DiaADC')
		self.adc_vect = self.dicBraVectChs['rawTree.DiaADC']
		self.CheckExistanceVector('commonModeNoise')
		self.cm_vect = self.dicBraVect1ch['commonModeNoise']
		temp = []
		self.CreateProgressBar(self.entries)
		if self.bar is not None:
			self.bar.start()
		for ev in xrange(self.entries):
			temp.append(self.adc_vect[ev].mean())
		self.cm_adc = np.array(temp, 'f8')
		del temp

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
