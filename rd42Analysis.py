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
import os, sys, shutil
from Utils import *
import subprocess as subp

__author__ = 'DA'

dicTypes = {'Char_t': 'int8', 'UChar_t': 'uint8', 'Short_t': 'short', 'UShort_t': 'ushort', 'Int_t': 'int32', 'UInt_t': 'uint32', 'Float_t': 'float32', 'Double_t': 'float64', 'Long64_t': 'int64',
            'ULong64_t': 'uint64', 'Bool_t': 'bool'}

diaChs = 128

fillColor = ro.TColor.GetColor(125, 153, 209)
sigma_axis = {'min': 0, 'max': 35}
adc_axis = {'min': 0, 'max': 2**12 - 1}
ped_axis = {'min': 0, 'max': 2**12 - 1}
cm_axis = {'min': -100, 'max': 100}

class RD42Analysis:
	def __init__(self, run, input, output, numev, runlistdir, settings_dir, pedestal=False, cluster=False, selec=False, autom_all=False):
		print 'Creating NoiseExtraction instance for run:', run
		self.run = run
		self.in_dir = input
		self.out_dir = output
		self.num_ev = numev
		self.runlist_dir = runlistdir
		self.settings_dir = settings_dir
		self.do_pedestal_ana = pedestal
		self.do_cluster_ana = cluster
		self.do_selec_ana = selec
		self.do_autom_all = autom_all
		self.bar = None
		self.process_f = None
		self.process_e = None
		self.process_o = None

	def Create_Run_List(self):
		if not os.path.isdir(self.runlist_dir):
			os.makedirs(self.runlist_dir)
		ped = 1 if self.do_pedestal_ana else 0
		clu = 1 if self.do_cluster_ana else 0
		sele = 1 if self.do_selec_ana else 0
		with open(self.runlist_dir + '/RunList_'+self.run+'.ini', 'w') as rlf:
			rlf.write('{r}\t0\t0\t{n}\t0\t{p}\t{c}\t{s}\t1\t0\t1\n'.format(r=self.run, n=self.num_ev, p=ped, c=clu, s=sele))
		if not os.path.isdir(self.runlist_dir+'/odd'):
			os.makedirs(self.runlist_dir+'/odd')
		with open(self.runlist_dir+'/odd' + '/RunList_'+self.run+'.ini', 'w') as rlf:
			rlf.write('{r}\t0\t0\t{n}\t0\t0\t0\t0\t0\t0\t0\n'.format(r=self.run, n=self.num_ev, p=ped, c=clu, s=sele))
		if not os.path.isdir(self.runlist_dir+'/even'):
			os.makedirs(self.runlist_dir+'/even')
		with open(self.runlist_dir+'/even' + '/RunList_'+self.run+'.ini', 'w') as rlf:
			rlf.write('{r}\t0\t0\t{n}\t0\t0\t0\t0\t0\t0\t0\n'.format(r=self.run, n=self.num_ev, p=ped, c=clu, s=sele))

	def Copy_settings_to_even_odd(self):
		self.Check_settings_full()
		if not os.path.isdir(self.settings_dir + '/even'):
			os.makedirs(self.settings_dir + '/even')
		if not os.path.isdir(self.settings_dir + '/odd'):
			os.makedirs(self.settings_dir + '/odd')
		shutil.copy(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run), self.settings_dir + '/even/')
		shutil.copy(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run), self.settings_dir + '/odd/')

	def Check_settings_full(self):
		if not os.path.isdir(self.settings_dir + '/full'):
			os.makedirs(self.settings_dir + '/full')
		if not os.path.isfile(self.settings_dir + '/full/settings.{r}.ini'.format(r=self.run)):
			CreateDefaultSettingsFile(self.settings_dir + '/full', self.run, self.num_ev)

	def Modify_even_odd(self):
		self.Modify_even()
		self.Modify_odd()

	def Modify_even(self):
		with open(self.settings_dir + '/even/settings.{r}.ini'.format(r=self.run), 'r') as fin:
			with open(self.settings_dir + '/even/settings.new.{r}.ini.'.format(r=self.run), 'w') as fnew:
				for line in fin:
					if not line.startswith('Dia_channel_screen_channels'):
						fnew.write(line)
					else:
						channel_str = line[line.find('{') + 1:line.find('}')].split(',')
						channel_str_new = ''
						for i in xrange(1, len(channel_str)):
							if IsInt(channel_str[i - 1]):
								prev = int(channel_str[i - 1])
								if IsInt(channel_str[i]):
									th = int(channel_str[i])
									if th - prev < 2:
										channel_str_new += str(prev) + ','
									elif prev % 2 == 0:
										for ch in xrange(prev, th + 1):
											channel_str_new += str(ch) + ','
										if th
										channel_str_new += str()
	def Modify_odd(self):

	def First_analysis(self):
		if not os.path.isdir(self.runlist_dir):
			os.makedirs(self.runlist_dir)
		with open(self.runlist_dir + '/RunList_'+self.run+'.ini', 'w') as rlf:
			rlf.write('{r}\t0\t0\t{n}\t0\t0\t0\t0\t0\t0\t0\n'.format(r=self.run, n=self.num_ev))
		CreateDefaultSettingsFile(self.settings_dir + '/no_mask', self.run, self.num_ev)
		self.process_f = subp.Popen(['diamondAnalysis', '-r {d}/RunList_{r}.ini'.format(d=self.runlist_dir, r=self.run), '-s ' + self.settings_dir + '/no_mask', '-o ' + self.out_dir, '-i ' + self.in_dir], bufsize=-1, stdin=subp.PIPE, close_fds=True)
		while self.process_f.poll() is None:
			pass
		CloseSubprocess(self.process_f, stdin=True, stdout=False)

	def ClearOldAnalysis(self):
		if os.path.isdir('{d}/pedestalAnalysis/histos'.format(d=self.dir)):
			shutil.rmtree('{d}/pedestalAnalysis/histos'.format(d=self.dir))
		if os.path.isdir('{d}/pedestalAnalysis/profiles'.format(d=self.dir)):
			shutil.rmtree('{d}/pedestalAnalysis/profiles'.format(d=self.dir))

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
			temp.SetTitle(''+name[:-6]+'_ch_'+str(ch))
			temp.GetXaxis().SetTitle('event')
			temp.GetYaxis().SetTitle(name[:-6])
			temp.SetFilColor(fillColor)
			temp.GetYaxis().SetRangeUser(prof.GetZaxis().GetXmin(), prof.GetZaxis().GetXmax())


	def SetProfileDefaults(self):
		for branch, prof in self.dicBraProf.iteritems():
			prof.GetXaxis().SetTitle('event')
			prof.GetYaxis().SetTitle('channel')
			prof.GetZaxis().SetTitle(self.dicBraNames[branch][:-6])
			prof.GetZaxis().SetRangeUser(self.dicBraHist[branch].GetZaxis().GetXmin(), self.dicBraHist[branch].GetZaxis().GetXmax())

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
		for branch, histo in self.dicBraHist.iteritems():
			if branch in self.listBraNamesChs:
				self.dicBraProf[branch] = histo.Project3DProfile('yx')
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
					ymin, ymax, ybins = cm_axis['min'], cm_axis['max'], int((cm_axis['max'] - cm_axis['min'])/0.5)
				self.dicBraHist[branch] = ro.TH2D(name, name, self.ev_axis['bins'], self.ev_axis['min'], self.ev_axis['max'], ybins + 1, ymin, ymax + (ymax - ymin) / float(ybins))
				self.dicBraHist[branch].GetXaxis().SetTitle('event')
				self.dicBraHist[branch].GetYaxis().SetTitle(name[:-6])
			elif branch in self.listBraNamesChs:
				zmin, zmax, zbins = 0, 0, 0
				if name.startswith('ped'):
					zmin, zmax, zbins = ped_axis['min'], ped_axis['max'], int((ped_axis['max'] - ped_axis['min'])/10.0)
				elif name.startswith('adc'):
					zmin, zmax, zbins = adc_axis['min'], adc_axis['max'], int((adc_axis['max'] - adc_axis['min'])/10.0)
				elif name.startswith('sigma'):
					zmin, zmax, zbins = sigma_axis['min'], sigma_axis['max'], int((sigma_axis['max'] - sigma_axis['min'])/0.1)
				self.dicBraHist[branch] = ro.TH3D(name, name, self.ev_axis['bins'], self.ev_axis['min'], self.ev_axis['max'], self.ch_axis['bins'], self.ch_axis['min'], self.ch_axis['max'], zbins + 1, zmin, zmax + (zmax - zmin) / float(zbins))
				self.dicBraHist[branch].GetXaxis().SetTitle('event')
				self.dicBraHist[branch].GetYaxis().SetTitle('channel')
				self.dicBraHist[branch].GetZaxis().SetTitle(name[:-6])
		self.adc_hist, self.ped_hist, self.sigma_hist, self.ped_cmc_hist, self.sigma_cmc_hist, self.cm_hist = self.dicBraHist['rawTree.DiaADC'], self.dicBraHist['diaPedestalMean'], self.dicBraHist['diaPedestaSigma'], self.dicBraHist['diaPedestalMeanCMN'], self.dicBraHist['diaPedestaSigmaCMN'], self.dicBraHist['commonModeNoise']
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
	parser.add_option('-r', '--run', dest='run', default=22011, type='int', help='Run to be analysed (e.g. 22011)')
	parser.add_option('-i', '--input', dest='input', default='.', type='string', help='source folder containing the run\'s binraries files')
	parser.add_option('-o', '--output', dest='output', default='.', type='string', help='foler where to save the structures (will contain, full, even and odd subfolders)')
	parser.add_option('-l', '--runlistdir', dest='runlist', default='~/RunLists', help='folder which contains the RunLists files')
	parser.add_option('-s', '--settingsdir', dest='sett', default='~/settings', help='folder which contains the settings files')
	parser.add_option('-n', '--numevents', dest='numevents', default=400000, type='int', help='number of events to analyse')
	parser.add_option('-p', '--pedestal', dest='pedestal', default=False, action='store_true', help='enables pedestal analysis')
	parser.add_option('-c', '--cluster', dest='cluster', default=False, action='store_true', help='enables cluster analysis')
	parser.add_option('-e', '--selection', dest='selection', default=False, action='store_true', help='enables selection analysis')
	parser.add_option('-a', '--all', dest='all', default=False, action='store_true', help='enables all types of analysis. Creates odd, even and full')
	parser.add_option('-x', '--first', dest='first', default=False, action='store_true', help='enables first analysis wich has everything un-masked')

	(options, args) = parser.parse_args()
	run = int(options.run)
	input = str(options.input)
	output = bool(options.output)
	numev = int(options.numevents)
	runlist = str(options.runlist)
	settings_dir = str(options.sett)
	pedestal = bool(options.pedestal)
	cluster = bool(options.cluster)
	selec = bool(options.selection)
	autom_all = bool(options.all)
	first_ana = bool(options.first)

	rd42 = RD42Analysis(run=run, input=input, output=output, numev=numev, runlistdir=runlist, settings_dir=settings_dir, pedestal=pedestal, cluster=cluster, selec=selec, autom_all=autom_all)

	if first_ana:
		rd42.First_analysis()
	elif autom_all:
		pass
	else:
		pass # run analysis for full with existing runlist file
	# output = str(options.output)
	# connect = int(options.connect)
	# low = int(options.low)
	# high = int(options.high)

