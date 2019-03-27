#!/usr/bin/env python
import ROOT as ro
import numpy as np
import ipdb, os
from ConfigParser import ConfigParser
from TransparentGrid import TransparentGrid
from ClusterChannelsAnalysis import ClusterChannelsAnalysis
from NoiseAnalysis import NoiseAnalysis
from CenterCellAnalysis import CenterCellAnalysis
from optparse import OptionParser
from Utils import *

color_index = 10000

class TestAreas:
	def __init__(self, configfile='', run=0):
		self.do_fit = False
		self.do_saturation = True
		self.num = 0
		self.clust_size = 2
		self.dir = '.'
		self.run = run
		self.cellsize = 50
		self.cell_resolution = 0
		self.phdelta = 0
		self.phmin = 10000
		self.phmax = -10000
		self.binsperx = 0
		self.binspery = 0
		self.threshold = 800
		self.skip_before_sat = 0
		self.skip_after_sat = 1
		self.do_threshold = False
		self.window_shift = 3
		self.min_snr_neg, self.max_snr_neg, self.delta_snr = -65, 1, 2
		# self.min_snr_neg, self.max_snr_neg, self.delta_snr = -64.25, 0.25, 0.125
		self.min_snr, self.max_snr = -650, 650
		self.min_adc, self.max_adc = -6500, 6500
		self.delta_adc = 20
		self.min_adc_noise, self.max_adc_noise, self.delta_adc_noise = -322.5, 322.5, 0.5
		# self.min_adc_noise, self.max_adc_noise, self.delta_adc_noise = -32.25, 32.25, 0.5
		# self.min_snr_noise, self.max_snr_noise, self.delta_snr_noise = -3.225, 3.225, 0.05
		self.min_snr_noise, self.max_snr_noise, self.delta_snr_noise = -32.25, 32.25, 0.05
		self.minz = {t: {} for t in ['all', 'good', 'bad']}
		self.maxz = {t: {} for t in ['all', 'good', 'bad']}
		self.neg_cut_lines = {}
		self.trash = []
		self.w = 0
		self.num_rows = 0
		self.rows = []
		self.cols = []
		self.cells = []
		self.rcells = []
		self.config_file = configfile
		if self.config_file != '':
			self.config_file = Correct_Path(self.config_file)
			self.ReadConfigFile()
		self.trans_grid = TransparentGrid(self.dir, self.run, self.cellsize)
		self.trans_grid.pkl_sbdir = 'test' + str(self.num)
		self.cluster_ch_ana = None
		self.noise_ana = None
		self.center_cells_ana = None
		self.bias = self.trans_grid.bias
		self.saturated_ADC = self.trans_grid.saturated_ADC
		self.num_strips = self.trans_grid.num_strips if self.trans_grid.num_strips != 0 else 3
		self.cluster_size = self.trans_grid.cluster_size if self.trans_grid.cluster_size != 0 else 3
		self.suffix = {'all': 'all', 'good': 'selection', 'bad': 'not_selection'}
		if self.num_rows != 0:
			self.trans_grid.row_info_diamond['num'] = self.num_rows

		self.noise_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.ph_adc_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.ph_snr_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.ph_adc_h_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.ph_snr_h_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_adc_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_adc_h_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_snr_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_snr_h_cuts = {t: {} for t in ['all', 'good', 'bad']}

		self.noise_varz = {}
		self.ph_adc_h_varz = {}
		self.ph_adc_ch_varz = {}
		self.ph_snr_h_varz = {}
		self.ph_snr_ch_varz = {}

		self.phN_adc_h_varz = {}
		self.phN_adc_ch_varz = {}
		self.phN_snr_h_varz = {}
		self.phN_snr_ch_varz = {}

	def ReadConfigFile(self):
		def unpack_row_col(string):
			elements = string.replace('{', '').replace('}', '')
			elements = elements.split(';')
			elements = [elem.split(',') for elem in elements]
			elements = [[int(elemij) for elemij in elemi] for elemi in elements]
			return elements
		if os.path.isfile(self.config_file):
			pars = ConfigParser()
			pars.read(self.config_file)
			print 'Reading config file for test area...', ; sys.stdout.flush()
			if pars.has_section('SETTINGS'):
				if pars.has_option('SETTINGS', 'run'):
					self.run = pars.getint('SETTINGS', 'run')
				if pars.has_option('SETTINGS', 'dir'):
					self.dir = Correct_Path(pars.get('SETTINGS', 'dir'))
				if pars.has_option('SETTINGS', 'cluster_size'):
					self.cluster_size = pars.getint('SETTINGS', 'cluster_size')
				if pars.has_option('SETTINGS', 'cell_size'):
					self.cellsize = pars.getint('SETTINGS', 'cell_size')
				if pars.has_option('SETTINGS', 'do_fit'):
					self.do_fit = pars.getboolean('SETTINGS', 'do_fit')
				if pars.has_option('SETTINGS', 'do_saturation'):
					self.do_saturation = pars.getboolean('SETTINGS', 'do_saturation')
				if pars.has_option('SETTINGS', 'skip_before_sat'):
					self.skip_before_sat = pars.getint('SETTINGS', 'skip_before_sat')
				if pars.has_option('SETTINGS', 'skip_after_sat'):
					self.skip_after_sat = pars.getint('SETTINGS', 'skip_after_sat')
				if pars.has_option('SETTINGS', 'test_number'):
					self.num = pars.getint('SETTINGS', 'test_number')
				if pars.has_option('SETTINGS', 'threshold'):
					self.threshold = pars.getfloat('SETTINGS', 'threshold')
				if pars.has_option('SETTINGS', 'do_threshold'):
					self.do_threshold = pars.getboolean('SETTINGS', 'do_threshold')
				if pars.has_option('SETTINGS', 'cell_resolution'):
					self.cell_resolution = pars.getfloat('SETTINGS', 'cell_resolution')
				if pars.has_option('SETTINGS', 'phdelta'):
					self.phdelta = pars.getfloat('SETTINGS', 'phdelta')
				if pars.has_option('SETTINGS', 'phmin'):
					self.phmin = pars.getfloat('SETTINGS', 'phmin')
				if pars.has_option('SETTINGS', 'phmax'):
					self.phmax = pars.getfloat('SETTINGS', 'phmax')
				if pars.has_option('SETTINGS', 'binsperx'):
					self.binsperx = pars.getfloat('SETTINGS', 'binsperx')
				if pars.has_option('SETTINGS', 'binspery'):
					self.binspery = pars.getfloat('SETTINGS', 'binspery')
				if pars.has_option('SETTINGS', 'efficiency_subdiv'):
					self.binspery = pars.getfloat('SETTINGS', 'efficiency_subdiv')
			if pars.has_section('ROWS'):
				if pars.has_option('ROWS', 'rows'):
					rows = pars.get('ROWS', 'rows')
					self.rows = unpack_row_col(rows)
				if pars.has_option('ROWS', 'num'):
					self.num_rows = pars.getint('ROWS', 'num')
			if pars.has_section('COLUMNS'):
				if pars.has_option('COLUMNS', 'cols'):
					cols = pars.get('COLUMNS', 'cols')
					self.cols = unpack_row_col(cols)
			if pars.has_section('CELLS'):
				if pars.has_option('CELLS', 'cells'):
					cells = pars.get('CELLS', 'cells')
					self.cells = unpack_row_col(cells)
			if pars.has_section('REMOVECELLS'):
				if pars.has_option('REMOVECELLS', 'rcells'):
					rcells = pars.get('REMOVECELLS', 'rcells')
					self.rcells = unpack_row_col(rcells)
			print 'Done'

	def SetTest(self):
		if self.do_threshold:
			self.trans_grid.ResetAreas()
			print 'Selecting areas with a ph2 greater or equal than', self.threshold, '...', ; sys.stdout.flush()
			self.trans_grid.SelectGoodAndBadByThreshold(self.threshold, 'clusterCharge2')
			print 'Done'
			self.trans_grid.AddRemainingToBadAreas()
			print 'Marked the remaining cells as bad'
			self.trans_grid.gridAreas.SimplifyGoodAndBadAreas()
			self.SetAnalysis()
		elif len(self.rows) + len(self.cols) + len(self.cells) > 0:
			self.trans_grid.ResetAreas()
			for row in self.rows:
				self.trans_grid.AddGoodAreasRow(row[0], row[1], row[2])
				print 'Added row {r} from {coli} to {colj} to selection'.format(r=row[0], coli=row[1], colj= row[2])
			for col in self.cols:
				self.trans_grid.AddGoodAreasCol(col[0], col[1], col[2])
				print 'Added column {r} from {coli} to {colj} to selection'.format(r=col[0], coli=col[1], colj=col[2])
			for cell in self.cells:
				self.trans_grid.AddGoodAreas(cell[0], cell[1])
				print 'Added cell with column {c} and row {r} to selection'.format(c=cell[0], r=cell[1])
			for rcell in self.rcells:
				self.trans_grid.RemoveFromGoodArea(rcell[0], rcell[1])
				print 'Removed cell with column {c} and row {r} from selection'.format(c=rcell[0], r=rcell[1])
			self.trans_grid.AddRemainingToBadAreas()
			print 'Marked the remaining cells as bad'
			self.trans_grid.gridAreas.SimplifyGoodAndBadAreas()
			self.SetAnalysis()
		else:
			print 'Enter a correct settings file for the test area in variable config_file and re run ReadConfigFile before setting the test...'

	def SetAnalysis(self):
		self.trans_grid.SetVarz()
		self.SetCutsInCutManager()
		self.cluster_ch_ana = ClusterChannelsAnalysis(self.trans_grid, self.num_strips, self.cluster_size)
		self.noise_ana = NoiseAnalysis(self.trans_grid, self.num_strips, self.cluster_size)
		self.trans_grid.FindMaxMinVarz()
		self.center_cells_ana = CenterCellAnalysis(self.trans_grid, self.num_strips, self.cluster_size)

	def SetCutsInCutManager(self):
		print 'Setting cuts in cut manager...', ; sys.stdout.flush()
		self.trans_grid.cuts_man.SetCells(selection=self.trans_grid.gridAreas.goodAreasCutNames_simplified_diamond, not_selection=self.trans_grid.gridAreas.badAreasCutNames_simplified_diamond)
		self.trans_grid.cuts_man.SetNoiseCuts()
		self.trans_grid.cuts_man.SetPHCuts()
		print 'Done'

	def PositionCanvas(self, canvas_name):
		if canvas_name in self.trans_grid.canvas.keys():
			self.trans_grid.canvas[canvas_name].SetWindowPosition(self.w, self.w)
			ro.gPad.Update()
			self.w += self.window_shift

	def DoBorderPlots(self):
		self.trans_grid.DrawProfile2DDiamond('PH2_H_map_with_borders', self.phN_adc_h_varz['PH2_H'], draw_top_borders=True)
		self.trans_grid.DrawGoodAreasDiamond('PH2_H_map_with_borders')
		self.trans_grid.DrawBadAreasDiamond('PH2_H_map_with_borders')
		xbinmin, xbinmax = self.trans_grid.profile['PH2_H_map_with_borders'].GetXaxis().FindBin(self.trans_grid.ch_ini - 0.5), self.trans_grid.profile['PH2_H_map_with_borders'].GetXaxis().FindBin(self.trans_grid.ch_ini - 0.5) + self.trans_grid.num_cols * self.trans_grid.bins_per_ch_x - 1
		self.trans_grid.profile['PH2_H_map_with_borders'].GetXaxis().SetRange(xbinmin - 1, xbinmax + 1)
		self.PositionCanvas('PH2_H_map_with_borders')
		self.trans_grid.canvas['PH2_H_map_with_borders_py'] = ro.TCanvas('c_PH2_H_map_with_borders_py', 'c_PH2_H_map_with_borders_py', 1)
		self.trans_grid.canvas['PH2_H_map_with_borders_py'].cd()
		self.trans_grid.histo['PH2_H_map_with_borders_py'] = self.trans_grid.profile['PH2_H_map_with_borders'].ProjectionY('h_PH2_H_map_with_borders_py', xbinmin, xbinmax, 'e hist')
		minbiny, maxbiny = self.trans_grid.histo['PH2_H_map_with_borders_py'].FindFirstBinAbove(), self.trans_grid.histo['PH2_H_map_with_borders_py'].FindLastBinAbove()
		for biny in xrange(maxbiny, int(self.trans_grid.histo['PH2_H_map_with_borders_py'].GetXaxis().GetNbins())):
			if self.trans_grid.histo['PH2_H_map_with_borders_py'].GetBinContent(biny) != 0:
				maxbiny = biny
		miny, maxy = self.trans_grid.histo['PH2_H_map_with_borders_py'].GetXaxis().GetBinLowEdge(minbiny), self.trans_grid.histo['PH2_H_map_with_borders_py'].GetXaxis().GetBinLowEdge(maxbiny + 1)
		self.trans_grid.profile['PH2_H_map_with_borders'].GetYaxis().SetRangeUser(miny, maxy)
		ro.gPad.Update()
		self.trans_grid.histo['PH2_H_map_with_borders_py'].GetXaxis().SetRangeUser(miny, maxy)
		func = ro.TF1('box_fcn', '[0]*(TMath::Erf((x-({l}))/[1])+1)/2-[2]*(TMath::Erf((x-{u})/[3])+1)/2+[4]'.format(l=self.trans_grid.row_info_diamond['0'], u=self.trans_grid.row_info_diamond['up']), miny, maxy)
		func.SetNpx(int(self.trans_grid.row_info_diamond['num'] * self.trans_grid.bins_per_ch_y * 10))
		zmin, zmax = self.trans_grid.histo['PH2_H_map_with_borders_py'].GetMinimum(), self.trans_grid.histo['PH2_H_map_with_borders_py'].GetMaximum()
		y1bin, y2bin = self.trans_grid.histo['PH2_H_map_with_borders_py'].FindFirstBinAbove((zmin + zmax) / 2.0), self.trans_grid.histo['PH2_H_map_with_borders_py'].FindLastBinAbove((zmin + zmax) / 2.0) + 1
		y1, y2 = self.trans_grid.histo['PH2_H_map_with_borders_py'].GetXaxis().GetBinCenter(y1bin), self.trans_grid.histo['PH2_H_map_with_borders_py'].GetXaxis().GetBinCenter(y2bin)
		z0, z1, z2 = self.trans_grid.histo['PH2_H_map_with_borders_py'].GetBinContent(int((minbiny))), self.trans_grid.histo['PH2_H_map_with_borders_py'].GetBinContent(int((y1bin + y2bin) / 2.0)), self.trans_grid.histo['PH2_H_map_with_borders_py'].GetBinContent(int((maxbiny)))
		func.SetParLimits(0, abs(z1 - z0) / 10.0, 2.0 * abs(z1 - z0))
		func.SetParLimits(1, 0.1, 50)
		func.SetParLimits(2, abs(z1 - z2) / 10.0, 2.0 * abs(z1 - z2))
		func.SetParLimits(3, 0.1, 50)
		func.SetParLimits(4, -2.0 * abs(z0), 10 * abs(z0))
		params = np.array((abs(z1 - z0), 20, abs(z1 - z2), 20, z0), 'float64')
		func.SetParameters(params)
		self.trans_grid.fits['PH2_H_map_with_borders_py'] = self.trans_grid.histo['PH2_H_map_with_borders_py'].Fit('box_fcn', 'QIEBMS', 'goff', self.trans_grid.histo['PH2_H_map_with_borders_py'].GetBinLowEdge(int((minbiny))), self.trans_grid.histo['PH2_H_map_with_borders_py'].GetBinLowEdge(int((maxbiny))))
		self.PositionCanvas('PH2_H_map_with_borders_py')

	def PlotNoiseNotInCluster(self, cells='all'):
		temp_cut_noise = self.noise_cuts[cells]
		temph = ro.TH1F('temph0', 'temph0', int(RoundInt((self.max_adc_noise - self.min_adc_noise) / float(self.delta_adc_noise))), self.min_adc_noise, self.max_adc_noise)
		self.trans_grid.trans_tree.Draw('diaChSignal>>temph0', temp_cut_noise, 'goff')
		mean, sigma = temph.GetMean(), temph.GetRMS()
		temph.Delete()
		self.min_snr_noise, self.max_snr_noise, self.delta_snr_noise = (ni / float(sigma) for ni in [self.min_adc_noise, self.max_adc_noise, self.delta_adc_noise])
		suffix = self.suffix[cells]
		self.trans_grid.DrawHisto1D('signal_noise_{c}_snr'.format(c=suffix), self.min_snr_noise, self.max_snr_noise, self.delta_snr_noise, self.noise_varz['snr'], varname='Signal not in cluster (SNR)', cuts=temp_cut_noise, option='e hist')
		self.trans_grid.FitGaus('signal_noise_{c}_snr'.format(c=suffix))
		self.trans_grid.histo['signal_noise_{c}_snr'.format(c=suffix)].GetXaxis().SetRangeUser(-3.2, 3.2)
		self.PositionCanvas('signal_noise_{c}_snr'.format(c=suffix))
		self.trans_grid.DrawHisto1D('signal_noise_{c}_adc'.format(c=suffix), self.min_adc_noise, self.max_adc_noise, self.delta_adc_noise, self.noise_varz['adc'], varname='Signal not in cluster (ADC)', cuts=temp_cut_noise, option='e hist')
		self.trans_grid.FitGaus('signal_noise_{c}_adc'.format(c=suffix))
		self.trans_grid.histo['signal_noise_{c}_adc'.format(c=suffix)].GetXaxis().SetRangeUser(-32, 32)
		self.PositionCanvas('signal_noise_{c}_adc'.format(c=suffix))

	def DoSaturationStudies(self, cells='all'):
		pass

	def OverlayNoiseDistribution(self, histo, cells='all'):
		suffix = self.suffix[cells]
		hname = histo.GetName().split('h_')[1]
		typ = 'adc' if 'adc' in hname.lower() else 'snr'
		noise_name0 = 'signal_noise_{s}_{t}'.format(s=suffix, t=typ)
		if not noise_name0 in self.trans_grid.histo.keys():
			self.PlotNoiseNotInCluster(cells)
		elif not self.trans_grid.histo[noise_name0]:
			del self.trans_grid.histo[noise_name0]
			self.PlotNoiseNotInCluster(cells)

		noise_name_new = noise_name0 + '_' + hname
		nbins = histo.GetNbinsX()
		xbins = np.zeros(nbins, 'float64')
		histo.GetXaxis().GetLowEdge(xbins)
		xbins = np.append(xbins, 2 * xbins[-1] - xbins[-2])
		self.trans_grid.histo[noise_name_new] = self.trans_grid.histo[noise_name0].Rebin(nbins, 'h_' + noise_name_new, xbins)
		if self.trans_grid.histo[noise_name_new]:
			self.trans_grid.histo[noise_name_new].SetTitle(noise_name0 + '(scaled)')
			scale = histo.GetMaximum() / self.trans_grid.histo[noise_name_new].GetMaximum()
			self.trans_grid.histo[noise_name_new].Scale(scale)
			self.trans_grid.histo[noise_name_new].SetLineColor(ro.kGray + 1)
			self.trans_grid.histo[noise_name_new].SetStats(0)
			self.trans_grid.canvas[hname].cd()
			if self.trans_grid.histo[noise_name_new].GetFunction('f_gaus_' + noise_name0):
				self.trans_grid.histo[noise_name_new].GetFunction('f_gaus_' + noise_name0).SetBit(ro.TF1.kNotDraw)
			self.trans_grid.histo[noise_name_new].Draw('same')
			ro.gPad.Update()

	def DoNoiseStudies(self, cells='all'):
		self.noise_ana.w, self.noise_ana.window_shift = self.w, self.window_shift
		self.noise_ana.DoNoiseAnalysis(cells)
		self.w, self.window_shift = self.noise_ana.w, self.noise_ana.window_shift

	def DoClusterStudies(self, cells='all'):
		self.cluster_ch_ana.w, self.cluster_ch_ana.window_shift = self.w, self.window_shift
		self.cluster_ch_ana.noise_ana = self.noise_ana
		self.cluster_ch_ana.DoClusterStudies(cells)
		self.w, self.window_shift = self.cluster_ch_ana.w, self.cluster_ch_ana.window_shift

	def DoCenterCellStudies(self, cells='all'):
		self.center_cells_ana.w, self.center_cells_ana.window_shift = self.w, self.window_shift
		dists = np.arange(0.2, 1, 0.05)
		percents = np.floor(np.power(dists, 2) * 100 + 0.5).astype('int32')
		for percent in percents:
			self.trans_grid.CreateTCutGSymmetricRectangle(percent)
		self.center_cells_ana.DoCenterRegionStudies(percents, cells)
		self.w, self.window_shift = self.center_cells_ana.w, self.center_cells_ana.window_shift

	def DoSaturationPlots(self, cells='all'):
		self.trans_grid.AddFriendWithSaturationRegions(self.skip_after_sat, self.skip_before_sat)

	def PlotSaturation(self):
		for c in xrange(1, self.cluster_size):
			self.trans_grid.DrawPHGoodAreas('ph{c}_saturated'.format(c=c), 'clusterCharge{c}'.format(c=c), '((transparentEvent)&&((diaChADC[clusterChannel0]=={s})||(diaChADC[clusterChannel1]=={s})||(diaChADC[clusterChannel2]=={s})))'.format(s=self.saturated_ADC))
			self.trans_grid.canvas['ph{c}_saturated'.format(c=c)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.DrawPHGoodAreas('ph{c}_not_saturated'.format(c=c), 'clusterCharge{c}'.format(c=c), '((transparentEvent)&&((diaChADC[clusterChannel0]<{s})&&(diaChADC[clusterChannel1]<{s})&&(diaChADC[clusterChannel2]<{s})))'.format(s=self.saturated_ADC))
			self.trans_grid.canvas['ph{c}_not_saturated'.format(c=c)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
		self.trans_grid.DrawPHGoodAreas('ph{c}_saturated'.format(c='N'), 'clusterCharge{c}'.format(c='N'), '((transparentEvent)&&((diaChADC[clusterChannel0]=={s})||(diaChADC[clusterChannel1]=={s})||(diaChADC[clusterChannel2]=={s})))'.format(s=self.saturated_ADC))
		self.trans_grid.canvas['ph{c}_saturated'.format(c='N')].SetWindowPosition(self.w, self.w)
		self.w += self.window_shift
		self.trans_grid.DrawPHGoodAreas('ph{c}_not_saturated'.format(c='N'), 'clusterCharge{c}'.format(c='N'), '((transparentEvent)&&((diaChADC[clusterChannel0]<{s})&&(diaChADC[clusterChannel1]<{s})&&(diaChADC[clusterChannel2]<{s})))'.format(s=self.saturated_ADC))
		self.trans_grid.canvas['ph{c}_not_saturated'.format(c='N')].SetWindowPosition(self.w, self.w)
		self.w += self.window_shift

	def PlotTestForNegative(self, cells='all'):
		y0, rowpitch, numrows, xoff, yoff, colpitch, numcols, yup = self.trans_grid.row_info_diamond['0'], self.trans_grid.row_info_diamond['pitch'], self.trans_grid.row_info_diamond['num'], self.trans_grid.row_info_diamond['x_off'], self.trans_grid.row_info_diamond['y_off'], self.trans_grid.col_pitch, self.trans_grid.num_cols, self.trans_grid.row_info_diamond['up']
		list_cuts_clusters = {i: ['(({y0}<diaChYPred)&&(diaChYPred<{yup})&&(diaChSignal/(diaChPedSigmaCmc+1e-12)<{m})&&(diaChPedSigmaCmc>0)&&(diaChannels==clusterChannel{n}))'.format(y0=y0, yup=yup, m=self.max_snr_neg, n=i)] for i in xrange(self.cluster_size)}
		list_cuts_noise = ['(({y0}<diaChYPred)&&(diaChYPred<{yup})&&(diaChHits==0)&&(diaChSignal/(diaChPedSigmaCmc+1e-12)<{m})&&(diaChPedSigmaCmc>0))'.format(y0=y0, yup=yup, m=self.max_snr_neg)]
		if cells == 'good':
			list_cuts_noise.append(self.trans_grid.gridAreas.goodAreasCutNames_simplified_diamond)
		elif cells == 'bad':
			list_cuts_noise.append(self.trans_grid.gridAreas.badAreasCutNames_simplified_diamond)
		temp_cut_clusters = {i: '&&'.join(list_cut) for i, list_cut in list_cuts_clusters.iteritems()}
		temp_cut_noise = '&&'.join(list_cuts_noise)
		lastbin = RoundInt((self.max_snr_neg - self.min_snr_neg) / float(self.delta_snr))
		tempmin, tempmax, tempbins = self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins
		suffix = 'not_selection' if cells == 'bad' else 'selection' if cells == 'good' else 'all'
		if cells == 'all':
			temp = [self.trans_grid.DrawHisto1D('ph_neg_c{n}_{c}'.format(n=i, c=suffix), self.min_snr_neg, self.max_snr_neg, self.delta_snr, var='diaChSignal/(diaChPedSigmaCmc+1e-12)', varname='PH closestStrip{n} (SNR)'.format(n=i), cuts=temp_cut_clusters[i]) for i in xrange(self.cluster_size)]
		elif cells == 'good':
			self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins = self.min_snr_neg, self.max_snr_neg, RoundInt((self.max_snr_neg - self.min_snr_neg) / float(self.delta_snr))
			temp = [self.trans_grid.DrawPHGoodAreas('ph_neg_c{n}_{c}'.format(n=i, c=suffix), 'diaChSignal/(diaChPedSigmaCmc+1e-12)', temp_cut_clusters[i], 'diamond', varname='PH closestStrip{n} (SNR)'.format(n=i)) for i in xrange(self.cluster_size)]
		elif cells == 'bad':
			self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins = self.min_snr_neg, self.max_snr_neg, RoundInt((self.max_snr_neg - self.min_snr_neg) / float(self.delta_snr))
			temp = [self.trans_grid.DrawPHBadAreas('ph_neg_c{n}_{c}'.format(n=i, c=suffix), 'diaChSignal/(diaChPedSigmaCmc+1e-12)', temp_cut_clusters[i], 'diamond', varname='PH closestStrip{n} (SNR)'.format(n=i)) for i in xrange(self.cluster_size)]
		self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins = tempmin, tempmax, tempbins
		temp = [self.trans_grid.DrawHisto1D('signal_noise_c{n}_{c}'.format(n=i, c=suffix), self.min_snr_neg, self.max_snr_neg, self.delta_snr, 'diaChSignal/(diaChPedSigmaCmc+1e-12)', varname='Signal not in cluster (SNR)', cuts=temp_cut_noise, option='goff') for i in xrange(self.cluster_size)]
		temp = [self.trans_grid.histo['signal_noise_c{n}_{c}'.format(n=i, c=suffix)].SetLineColor(ro.kGray + 1) for i in xrange(self.cluster_size)]
		temp = [self.trans_grid.histo['signal_noise_c{n}_{c}'.format(n=i, c=suffix)].Scale(self.trans_grid.histo['ph_neg_c{n}_{c}'.format(n=i, c=suffix)].GetBinContent(lastbin)/self.trans_grid.histo['signal_noise_c{n}_{c}'.format(n=i, c=suffix)].GetBinContent(lastbin)) for i in xrange(self.cluster_size)]
		temp = [self.trans_grid.histo['signal_noise_c{n}_{c}'.format(n=i, c=suffix)].SetStats(0) for i in xrange(self.cluster_size)]
		temp = [self.trans_grid.histo['signal_noise_c{n}_{c}'.format(n=i, c=suffix)].SetTitle('signal_not_in_cluster(scaled)') for i in xrange(self.cluster_size)]
		for c in xrange(self.cluster_size):
			ro.gPad.Update()
			self.trans_grid.canvas['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].SetWindowPosition(self.w, self.w)
			self.trans_grid.canvas['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].cd()
			self.trans_grid.histo['signal_noise_c{n}_{c}'.format(n=c, c=suffix)].Draw('same')
			ro.gPad.Update()
			if self.trans_grid.histo['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].FindObject('stats'):
				self.trans_grid.histo['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].FindObject('stats').SetX1NDC(0.15)
				self.trans_grid.histo['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].FindObject('stats').SetX2NDC(0.45)
				ro.gPad.Update()
			hbla = ro.TH1F('bla{n}_{c}'.format(n=c, c=suffix), 'Negative cut', 2, self.max_snr_neg + 1, self.max_snr_neg + 2)
			hbla.SetLineStyle(7)
			hbla.SetLineColor(ro.kRed)
			hbla.SetLineWidth(2)
			hbla.Draw('same')
			self.trash.append(hbla)
			ro.gPad.Update()
			legend = self.trans_grid.canvas['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].BuildLegend()
			ro.gPad.Update()
			legend.SetX1NDC(0.15)
			legend.SetX2NDC(0.45)
			legend.SetY1NDC(0.5)
			legend.SetY2NDC(0.6)
			ro.gPad.Update()
			self.neg_cut_lines['ph_neg_c{n}_{c}'.format(n=c, c=suffix)] = ro.TLine(-self.trans_grid.neg_cut, 0, -self.trans_grid.neg_cut, self.trans_grid.histo['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].GetMaximum())
			self.neg_cut_lines['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].SetLineColor(ro.kRed)
			self.neg_cut_lines['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].SetLineStyle(7)
			self.neg_cut_lines['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].SetLineWidth(2)
			self.neg_cut_lines['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].Draw('same')
			self.trans_grid.canvas['ph_neg_c{n}_{c}'.format(n=c, c=suffix)].SetLogy()
			ro.gPad.Update()
			self.w += self.window_shift

		self.trans_grid.canvas['ph_neg_{c}'.format(c=suffix)] = ro.TCanvas('ph_neg_{c}'.format(c=suffix), 'ph_neg_{c}'.format(c=suffix), 1)
		for c in xrange(self.cluster_size):
			self.trans_grid.histo['ph_neg_c{n}_{c}'.format(n=self.cluster_size - c - 1, c=suffix)].SetStats(0)
			self.trans_grid.histo['ph_neg_c{n}_{c}'.format(n=self.cluster_size - c - 1, c=suffix)].Draw('same')
			r, g, b = ReturnRGB(c, 0, self.cluster_size)
			self.trash.append(ro.TColor(color_index + c, r, g, b))
			self.trans_grid.histo['ph_neg_c{n}_{c}'.format(n=self.cluster_size - c - 1, c=suffix)].SetLineColor(color_index + c)
			self.trans_grid.histo['ph_neg_c{n}_{c}'.format(n=self.cluster_size - c - 1, c=suffix)].SetMarkerColor(ro.kBlack)

		hbla = ro.TH1F('bla000', 'Negative cut', 2, self.max_snr_neg + 1, self.max_snr_neg + 2)
		hbla.SetLineStyle(7)
		hbla.SetLineColor(ro.kRed)
		hbla.SetLineWidth(2)
		hbla.Draw('same')
		self.trash.append(hbla)
		self.trans_grid.canvas['ph_neg_{c}'.format(c=suffix)].SetLogy()
		legend = self.trans_grid.canvas['ph_neg_{c}'.format(c=suffix)].BuildLegend()
		ro.gPad.Update()
		legend.SetX1NDC(0.15)
		legend.SetX2NDC(0.45)
		legend.SetY1NDC(0.7)
		legend.SetY2NDC(0.9)
		ro.gPad.Update()
		self.neg_cut_lines['ph_neg_c{n}_{c}'.format(n=self.cluster_size - 1, c=suffix)].Draw('same')
		ro.gPad.Update()
		self.trans_grid.canvas['ph_neg_{c}'.format(c=suffix)].SetWindowPosition(self.w, self.w)
		self.w += self.window_shift

		for c in xrange(self.cluster_size):
			self.trans_grid.DrawProfile2DDiamond('ph_neg_map_c{n}'.format(n=c, c=suffix), 'diaChSignal[clusterChannel{i}]'.format(i=c), '(diaChSignal[clusterChannel{i}]/diaChPedSigmaCmc[clusterChannel{i}]<=-{c})'.format(i=c, c=self.trans_grid.neg_cut))
			# self.trans_grid.DrawProfile2DDiamond('ph_neg_map_c{n}'.format(n=c, c=suffix), 'diaChSignal', '(diaChannels==clusterChannel{i})'.format(i=c))
			self.trans_grid.profile['ph_neg_map_c{n}'.format(n=c, c=suffix)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile['ph_neg_map_c{n}'.format(n=c, c=suffix)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
			self.trans_grid.canvas['ph_neg_map_c{n}'.format(n=c, c=suffix)].SetWindowPosition(self.w, self.w)
			self.trans_grid.DrawTCutGs('ph_neg_map_c{n}'.format(n=c, c=suffix), 'diamond')
			self.trans_grid.DrawGoodAreasDiamondCenters('ph_neg_map_c{n}'.format(n=c, c=suffix))
			self.w += self.window_shift
			self.trans_grid.DrawProfile2DNoTopBottomBorders('signal_neg_map_c{n}'.format(n=c, c=suffix), -0.5, 127.5, 1, 'dia X ch', self.trans_grid.row_info_diamond['0']-RoundInt(self.trans_grid.row_info_diamond['0'] / self.trans_grid.row_info_diamond['pitch'], 'f8') * self.trans_grid.row_info_diamond['pitch'],
			                                                self.trans_grid.row_info_diamond['0'] + (256 - RoundInt(self.trans_grid.row_info_diamond['0'] / self.trans_grid.row_info_diamond['pitch'], 'f8')) * self.trans_grid.row_info_diamond['pitch'], float(self.trans_grid.row_info_diamond['pitch'])/self.trans_grid.bins_per_ch_y,
			                                                'sil pred Y [#mum]', 'diaChannels', 'diaChYPred', 'diaChSignal', 'Ch Signal [ADC]', '(diaChannels==clusterChannel{i})&&(diaChSignal/diaChPedSigmaCmc<=-{c})'.format(i=c, c=self.trans_grid.neg_cut))
			self.trans_grid.canvas['signal_neg_map_c{n}'.format(n=c, c=suffix)].SetWindowPosition(self.w, self.w)
			self.trans_grid.profile['signal_neg_map_c{n}'.format(n=c, c=suffix)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile['signal_neg_map_c{n}'.format(n=c, c=suffix)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
			self.trans_grid.DrawTCutGs('signal_neg_map_c{n}'.format(n=c, c=suffix), 'diamond')
			self.trans_grid.DrawGoodAreasDiamondCenters('signal_neg_map_c{n}'.format(n=c, c=suffix))
			self.w += self.window_shift
			self.trans_grid.DrawProfile2DDiamondCellOverlay('ph_neg_cells_c{n}_{c}'.format(n=c, c=suffix), 'diaChSignal', cells, '(diaChannels==clusterChannel{n})&&(diaChSignal/diaChPedSigmaCmc<=-{c})'.format(n=c, c=self.trans_grid.neg_cut))
			self.trans_grid.canvas['ph_neg_cells_c{n}_{c}'.format(n=c, c=suffix)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift


	def PlotTest(self):
		num = self.num
		y0, rowpitch, numrows, xoff, yoff, colpitch, numcols, yup = self.trans_grid.row_info_diamond['0'], self.trans_grid.row_info_diamond['pitch'], self.trans_grid.row_info_diamond['num'], self.trans_grid.row_info_diamond['x_off'], self.trans_grid.row_info_diamond['y_off'], self.trans_grid.col_pitch, self.trans_grid.num_cols, self.trans_grid.row_info_diamond['up']
		tempn = self.num_strips if self.num_strips == 1 else 2
		for ch in xrange(1, tempn + 1):
			#  plot map
			self.trans_grid.DrawProfile2DDiamond('ph{c}_map_test{n}'.format(c=ch, n=num), varz='clusterCharge' + str(ch), cuts='({l}<=diaChYPred)&&(diaChYPred<={h})'.format(l=y0, h=yup))
			# ro.gPad.Update()
			self.trans_grid.canvas['ph{c}_map_test{n}'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.gridTextDiamond.Draw('same TEXT0')
			self.trans_grid.profile['ph{c}_map_test'.format(c=ch) + str(num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile['ph{c}_map_test'.format(c=ch) + str(num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
			#  draw selected areas
			self.trans_grid.DrawGoodAreasDiamond('ph{c}_map_test{n}'.format(c=ch, n=num))
			self.trans_grid.DrawGoodAreasDiamondCenters('ph{c}_map_test{n}'.format(c=ch, n=num))
			#  plot cell overlay
			self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{c}_cells_test{n}'.format(c=ch, n=num), var='clusterCharge' + str(ch), cells='good')
			self.trans_grid.canvas['ph{c}_cells_test{n}'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.trans_grid.DrawTCutCentersInCellOverlay('ph{c}_cells_test{n}'.format(c=ch, n=num))
			self.w += self.window_shift
			self.trans_grid.GetOccupancyFromProfile('ph{c}_cells_test{n}'.format(c=ch, n=num), 'goff')
			#  show center region in cell
			#  Efficiency plots
			self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num), var='clusterCharge' + str(ch), cells='good', cuts='(clusterCharge{c}>=5*diaChPedSigmaCmc[clusterChannel0])'.format(c=ch), plot_option='prof colz goff')
			self.trans_grid.GetOccupancyFromProfile('ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num))
			self.trans_grid.canvas['hit_map_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.histo['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)] = self.trans_grid.histo['hit_map_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)].Clone('h_Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num))
			self.trans_grid.histo['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)].SetTitle('h_Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num))
			self.trans_grid.histo['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)].Sumw2()
			self.trans_grid.histo['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)].Divide(self.trans_grid.histo['hit_map_ph{c}_cells_test{n}'.format(c=ch, n=num)])
			self.trans_grid.canvas['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)] = ro.TCanvas('c_Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num), 'c_Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num), 1)
			self.trans_grid.histo['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)].Draw('colz')
			SetDefault2DStats(self.trans_grid.histo['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)])
			self.trans_grid.histo['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)].GetZaxis().SetTitle('Efficiency')
			self.trans_grid.histo['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)].GetZaxis().SetRangeUser(0.9, 1)
			self.trans_grid.canvas['Eff_ph{c}_cells_5sigma_cut_test{n}'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.DrawEfficiencyADCCut('Eff_ph{c}VsADC_test{n}'.format(c=ch, n=num), 'clusterCharge' + str(ch), cells='good', cut='transparentEvent', ymin_plot=0.95)
			self.trans_grid.canvas['Eff_ph{c}VsADC_test{n}'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			#  draw ph of selected areas
			self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}'.format(c=ch, n=num), var='clusterCharge' + str(ch))
			self.trans_grid.canvas['ph{c}_test{n}'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			if len(self.trans_grid.gridAreas.goodAreas_diamond_centers) < 900:
				self.trans_grid.DrawPHCentralRegion('ph{c}_test{n}_centers'.format(c=ch, n=num), cells='good', var='clusterCharge' + str(ch))
				self.trans_grid.canvas['ph{c}_test{n}_centers'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				#  fit distribution for central region
				if self.do_fit: self.trans_grid.FitLanGaus('ph{c}_test{n}_centers'.format(c=ch, n=num), color=ro.kRed)
				#  get difference between cell and center
				self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}_periphery'.format(c=ch, n=num))
				self.trans_grid.histo['ph{c}_test{n}_periphery'.format(c=ch, n=num)].Reset('ICES')
				self.trans_grid.histo['ph{c}_test{n}_periphery'.format(c=ch, n=num)].Add(self.trans_grid.histo['ph{c}_test{n}'.format(c=ch, n=num)], self.trans_grid.histo['ph{c}_test{n}_centers'.format(c=ch, n=num)], 1, -1)
				self.trans_grid.canvas['ph{c}_test{n}_periphery'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				if self.do_fit:
					self.trans_grid.FitLanGaus('ph{c}_test{n}_periphery'.format(c=ch, n=num), color=ro.kBlue)
					self.trans_grid.canvas['ph{c}_test{n}'.format(c=ch, n=num)].cd()
					self.trans_grid.langaus['ph{c}_test{n}_centers'.format(c=ch, n=num)].fit.Draw('same')
					self.trans_grid.langaus['ph{c}_test{n}_periphery'.format(c=ch, n=num)].fit.Draw('same')
					self.trans_grid.DrawDoubleLangaus('ph{c}_test{n}'.format(c=ch, n=num), 'ph{c}_test{n}_centers'.format(c=ch, n=num), 'ph{c}_test{n}_periphery'.format(c=ch, n=num), color=ro.kBlack)
			# ro.gPad.Update()
			#  position of negative clusters
			cut_no_neg = '(Sum$((diaChHits)&&(diaChSignal>-{c}*diaChPedSigmaCmc))=={n})'.format(c=self.trans_grid.neg_cut, n=self.cluster_size)
			cut_any_neg = '(Sum$((diaChHits)&&(diaChSignal<-{c}*diaChPedSigmaCmc))>0)'.format(c=self.trans_grid.neg_cut)
			self.trans_grid.DrawProfile2DDiamond('ph{c}_map_test{n}_negative'.format(c=ch, n=num), varz='clusterCharge' + str(ch), cuts=cut_any_neg)
			self.trans_grid.DrawGoodAreasDiamondCenters('ph{c}_map_test{n}_negative'.format(c=ch, n=num))
			self.trans_grid.profile['ph{c}_map_test{n}_negative'.format(c=ch, n=num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile['ph{c}_map_test{n}_negative'.format(c=ch, n=num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
			self.trans_grid.canvas['ph{c}_map_test{n}_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.GetOccupancyFromProfile('ph{c}_map_test{n}_negative'.format(c=ch, n=num))
			self.trans_grid.histo['hit_map_ph{c}_map_test{n}_negative'.format(c=ch, n=num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.histo['hit_map_ph{c}_map_test{n}_negative'.format(c=ch, n=num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
			self.trans_grid.canvas['hit_map_ph{c}_map_test{n}_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.trans_grid.DrawGoodAreasDiamond('hit_map_ph{c}_map_test{n}_negative'.format(c=ch, n=num))
			self.trans_grid.DrawBadAreasDiamond('hit_map_ph{c}_map_test{n}_negative'.format(c=ch, n=num))
			self.trans_grid.DrawGoodAreasDiamondCenters('hit_map_ph{c}_map_test{n}_negative'.format(c=ch, n=num))
			self.w += self.window_shift
			self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{c}_cells_test{n}_negative'.format(c=ch, n=num), var='clusterCharge' + str(ch), cells='good', cuts=cut_any_neg)
			self.trans_grid.canvas['ph{c}_cells_test{n}_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.GetOccupancyFromProfile('ph{c}_cells_test{n}_negative'.format(c=ch, n=num))
			self.trans_grid.canvas['hit_map_ph{c}_cells_test{n}_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}_negative'.format(c=ch, n=num), var='clusterCharge' + str(ch), cuts=cut_any_neg)
			self.trans_grid.canvas['ph{c}_test{n}_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			if self.do_fit: self.trans_grid.FitLanGaus('ph{c}_test{n}_negative'.format(c=ch, n=num), color=ro.kBlue)
			self.w += self.window_shift

			self.trans_grid.DrawProfile2DDiamond('ph{c}_map_test{n}_no_negative'.format(c=ch, n=num), varz='clusterCharge' + str(ch), cuts=cut_no_neg)
			self.trans_grid.DrawGoodAreasDiamondCenters('ph{c}_map_test{n}_no_negative'.format(c=ch, n=num))
			self.trans_grid.profile['ph{c}_map_test{n}_no_negative'.format(c=ch, n=num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile['ph{c}_map_test{n}_no_negative'.format(c=ch, n=num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
			self.trans_grid.canvas['ph{c}_map_test{n}_no_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.GetOccupancyFromProfile('ph{c}_map_test{n}_no_negative'.format(c=ch, n=num))
			self.trans_grid.histo['hit_map_ph{c}_map_test{n}_no_negative'.format(c=ch, n=num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.histo['hit_map_ph{c}_map_test{n}_no_negative'.format(c=ch, n=num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
			self.trans_grid.canvas['hit_map_ph{c}_map_test{n}_no_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.trans_grid.DrawGoodAreasDiamond('hit_map_ph{c}_map_test{n}_no_negative'.format(c=ch, n=num))
			self.trans_grid.DrawBadAreasDiamond('hit_map_ph{c}_map_test{n}_no_negative'.format(c=ch, n=num))
			self.trans_grid.DrawGoodAreasDiamondCenters('hit_map_ph{c}_map_test{n}_no_negative'.format(c=ch, n=num))
			self.w += self.window_shift
			self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{c}_cells_test{n}_no_negative'.format(c=ch, n=num), var='clusterCharge' + str(ch), cells='good', cuts=cut_no_neg)
			self.trans_grid.canvas['ph{c}_cells_test{n}_no_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.GetOccupancyFromProfile('ph{c}_cells_test{n}_no_negative'.format(c=ch, n=num))
			self.trans_grid.canvas['hit_map_ph{c}_cells_test{n}_no_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}_no_negative'.format(c=ch, n=num), var='clusterCharge' + str(ch), cuts=cut_no_neg)
			self.trans_grid.canvas['ph{c}_test{n}_no_negative'.format(c=ch, n=num)].SetWindowPosition(self.w, self.w)
			if self.do_fit: self.trans_grid.FitLanGaus('ph{c}_test{n}_no_negative'.format(c=ch, n=num), color=ro.kBlue)
			self.w += self.window_shift

		if self.clust_size >= 2:
			self.trans_grid.DrawProfile2DDiamondCellOverlay('ph2_minus_ph1_map_test{n}'.format(n=num), 'clusterCharge2-clusterCharge1', cells='good')
			self.trans_grid.canvas['ph2_minus_ph1_map_test{n}'.format(n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.DrawTCutCentersInCellOverlay('ph2_minus_ph1_map_test{n}'.format(n=num))
			tempmin, tempmax, tempbins = self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins
			self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins = self.trans_grid.phmin_neg, self.trans_grid.phmax_neg, self.trans_grid.phbins_neg
			self.trans_grid.DrawPHGoodAreas('ph2_minus_ph1_test{n}'.format(n=num), 'clusterCharge2-clusterCharge1')
			self.trans_grid.canvas['ph2_minus_ph1_test{n}'.format(n=num)].SetWindowPosition(self.w, self.w)
			self.w += self.window_shift
			self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins = tempmin, tempmax, tempbins

			if self.num_strips != 2:
				self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{s}_minus_ph2_map_test{n}'.format(s=self.num_strips, n=num), 'clusterChargeN-clusterCharge2', cells='good')
				self.trans_grid.canvas['ph{s}_minus_ph2_map_test{n}'.format(s=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				self.trans_grid.DrawTCutCentersInCellOverlay('ph{s}_minus_ph2_map_test{n}'.format(s=self.num_strips, n=num))
				tempmin, tempmax, tempbins = self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins
				self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins = self.trans_grid.phmin_neg, self.trans_grid.phmax_neg, self.trans_grid.phbins_neg
				self.trans_grid.DrawPHGoodAreas('ph{s}_minus_ph2_test{n}'.format(s=self.num_strips, n=num), 'clusterChargeN-clusterCharge2')
				self.trans_grid.canvas['ph{s}_minus_ph2_test{n}'.format(s=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins = tempmin, tempmax, tempbins

				self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{c}_cells_test{n}'.format(c=self.num_strips, n=num), var='clusterChargeN', cells='good')
				self.trans_grid.canvas['ph{c}_cells_test{n}'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				#  show center region in cell
				self.trans_grid.DrawTCutCentersInCellOverlay('ph{c}_cells_test{n}'.format(c=self.num_strips, n=num))
				#  draw ph of selected areas
				self.trans_grid.DrawEfficiencyADCCut('Eff_ph{c}VsADC_test{n}'.format(c=self.num_strips, n=num), 'clusterChargeN', cells='good', cut='transparentEvent', ymin_plot=0.95)
				self.trans_grid.canvas['Eff_ph{c}VsADC_test{n}'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}'.format(c=self.num_strips, n=num), var='clusterChargeN')
				self.trans_grid.canvas['ph{c}_test{n}'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				if len(t.trans_grid.gridAreas.goodAreas_diamond_centers) < 900:
					self.trans_grid.DrawPHCentralRegion('ph{c}_test{n}_centers'.format(c=self.num_strips, n=num), cells='good', var='clusterChargeN')
					self.trans_grid.canvas['ph{c}_test{n}_centers'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
					self.w += self.window_shift
					#  fit distribution for central region
					if self.do_fit: self.trans_grid.FitLanGaus('ph{c}_test{n}_centers'.format(c=self.num_strips, n=num), color=ro.kRed)
					#  get difference between cell and center
					self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}_periphery'.format(c=self.num_strips, n=num))
					self.trans_grid.histo['ph{c}_test{n}_periphery'.format(c=self.num_strips, n=num)].Reset('ICES')
					self.trans_grid.histo['ph{c}_test{n}_periphery'.format(c=self.num_strips, n=num)].Add(self.trans_grid.histo['ph{c}_test{n}'.format(c=self.num_strips, n=num)], self.trans_grid.histo['ph{c}_test{n}_centers'.format(c=self.num_strips, n=num)], 1, -1)
					self.trans_grid.canvas['ph{c}_test{n}_periphery'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
					self.w += self.window_shift
					if self.do_fit:
						self.trans_grid.FitLanGaus('ph{c}_test{n}_periphery'.format(c=self.num_strips, n=num), color=ro.kBlue)
						self.trans_grid.canvas['ph{c}_test{n}'.format(c=self.num_strips, n=num)].cd()
						self.trans_grid.langaus['ph{c}_test{n}_centers'.format(c=self.num_strips, n=num)].fit.Draw('same')
						self.trans_grid.langaus['ph{c}_test{n}_periphery'.format(c=self.num_strips, n=num)].fit.Draw('same')
						self.trans_grid.DrawDoubleLangaus('ph{c}_test{n}'.format(c=self.num_strips, n=num), 'ph{c}_test{n}_centers'.format(c=self.num_strips, n=num), 'ph{c}_test{n}_periphery'.format(c=self.num_strips, n=num), color=ro.kBlack)
				# ro.gPad.Update()
				#  position of negative clusters
				cut_no_neg = '(Sum$((diaChHits)&&(diaChSignal>-{c}*diaChPedSigmaCmc))=={n})'.format(c=self.trans_grid.neg_cut, n=self.cluster_size)
				cut_any_neg = '(Sum$((diaChHits)&&(diaChSignal<-{c}*diaChPedSigmaCmc))>0)'.format(c=self.trans_grid.neg_cut)
				self.trans_grid.DrawProfile2DDiamond('ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num), varz='clusterChargeN', cuts=cut_any_neg)
				self.trans_grid.DrawGoodAreasDiamondCenters('ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.profile['ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
				self.trans_grid.profile['ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
				self.trans_grid.canvas['ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				self.trans_grid.GetOccupancyFromProfile('ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.histo['hit_map_ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
				self.trans_grid.histo['hit_map_ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
				self.trans_grid.canvas['hit_map_ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.trans_grid.DrawGoodAreasDiamond('hit_map_ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.DrawBadAreasDiamond('hit_map_ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.DrawGoodAreasDiamondCenters('hit_map_ph{c}_map_test{n}_negative'.format(c=self.num_strips, n=num))
				self.w += self.window_shift
				self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{c}_cells_test{n}_negative'.format(c=self.num_strips, n=num), var='clusterChargeN', cells='good', cuts=cut_any_neg)
				self.trans_grid.canvas['ph{c}_cells_test{n}_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				self.trans_grid.GetOccupancyFromProfile('ph{c}_cells_test{n}_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.canvas['hit_map_ph{c}_cells_test{n}_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}_negative'.format(c=self.num_strips, n=num), var='clusterChargeN', cuts=cut_any_neg)
				self.trans_grid.canvas['ph{c}_test{n}_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				if self.do_fit: self.trans_grid.FitLanGaus('ph{c}_test{n}_negative'.format(c=self.num_strips, n=num), color=ro.kBlue)
				self.w += self.window_shift

				self.trans_grid.DrawProfile2DDiamond('ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num), varz='clusterChargeN', cuts=cut_no_neg)
				self.trans_grid.DrawGoodAreasDiamondCenters('ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.profile['ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
				self.trans_grid.profile['ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
				self.trans_grid.canvas['ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				self.trans_grid.GetOccupancyFromProfile('ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.histo['hit_map_ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
				self.trans_grid.histo['hit_map_ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
				self.trans_grid.canvas['hit_map_ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.trans_grid.DrawGoodAreasDiamond('hit_map_ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.DrawBadAreasDiamond('hit_map_ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.DrawGoodAreasDiamondCenters('hit_map_ph{c}_map_test{n}_no_negative'.format(c=self.num_strips, n=num))
				self.w += self.window_shift
				self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{c}_cells_test{n}_no_negative'.format(c=self.num_strips, n=num), var='clusterChargeN', cells='good', cuts=cut_no_neg)
				self.trans_grid.canvas['ph{c}_cells_test{n}_no_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				self.trans_grid.GetOccupancyFromProfile('ph{c}_cells_test{n}_no_negative'.format(c=self.num_strips, n=num))
				self.trans_grid.canvas['hit_map_ph{c}_cells_test{n}_no_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				self.w += self.window_shift
				self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}_no_negative'.format(c=self.num_strips, n=num), var='clusterChargeN', cuts=cut_no_neg)
				self.trans_grid.canvas['ph{c}_test{n}_no_negative'.format(c=self.num_strips, n=num)].SetWindowPosition(self.w, self.w)
				if self.do_fit: self.trans_grid.FitLanGaus('ph{c}_test{n}_no_negative'.format(c=self.num_strips, n=num), color=ro.kBlue)
				self.w += self.window_shift

	def SaveCanvas(self):
		self.trans_grid.SaveCanvasInlist(self.trans_grid.canvas.keys())

	def DoAutomatic(self, cells='good', do_save=True):
		self.DoBorderPlots()
		self.DoNoiseStudies(cells)
		self.DoClusterStudies(cells)
		# self.PlotTestClusterStudies(cells)
		# self.PlotTestForNegative(cells)
		# self.PlotTest()
		# if self.do_saturation:
		# 	self.PlotSaturation()
		if do_save:
			self.SaveCanvas()

	def SetTransparentGrid(self):
		self.trans_grid.SetLines()
		self.trans_grid.CreateTCutGs()
		if self.trans_grid.saturated_ADC != 0:
			self.saturated_ADC = self.trans_grid.saturated_ADC
		else:
			self.saturated_ADC = Get_From_User_Value('saturation_ADC for run ' + str(self.run), 'int', self.trans_grid.saturated_ADC, True)
			self.trans_grid.saturated_ADC = self.saturated_ADC
			self.trans_grid.SavePickle()
		if self.trans_grid.bias != 0:
			self.bias = self.trans_grid.bias
		else:
			self.bias = Get_From_User_Value('bias for run ' + str(self.run), 'float', self.trans_grid.bias, update=True)
			self.trans_grid.bias = self.bias
			self.trans_grid.SavePickle()
		if self.trans_grid.neg_cut_adc == 4100:
			self.trans_grid.DrawHisto1D('NoiseAllCells', -32.25, 32.25, 0.5, 'diaChSignal', 'signal not in cluster cmc [ADC]', self.trans_grid.cuts_man.not_in_transp_cluster, option='goff')
			self.trans_grid.neg_cut_adc = self.trans_grid.neg_cut * self.trans_grid.histo['NoiseAllCells'].GetRMS()
			print 'Setting neg_cut_adc to', self.trans_grid.neg_cut_adc, '. If you wish to change it, do it and save the pickle again in transparent grid object'
			self.trans_grid.SavePickle()

if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-r', '--run', dest='run', type='int', default=0, help='run number to be analysed (e.g. 25209)')
	parser.add_option('-a', '--auto', dest='auto', default=False, action='store_true', help='Sets up test, creates plots and saves them automatically if toggled')
	parser.add_option('-c', '--config', dest='config', default='', type='string', help='gives the path to a config file for the test area')

	(options, args) = parser.parse_args()
	run = int(options.run)
	autom = bool(options.auto)
	config = str(options.config)

	t = TestAreas(config, run)
	t.SetTransparentGrid()
	t.SetTest()
	if t.trans_grid.loaded_pickle:
		t.trans_grid.LoadPickle()
	else:
		t.trans_grid.FindPickleValues()
		ExitMessage('Run it again to load the generated pickle :)', os.EX_OK)
	if t.trans_grid.loaded_default_pickle:
		t.trans_grid.FindBinningAndResolution()
		t.trans_grid.SavePickle()
	if autom:
		t.DoAutomatic('good')

