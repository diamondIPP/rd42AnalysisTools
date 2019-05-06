#!/usr/bin/env python
import ROOT as ro
import numpy as np
import ipdb, os
from ConfigParser import ConfigParser
from TransparentGrid import TransparentGrid
from ClusterChannelsAnalysis import ClusterChannelsAnalysis
from NoiseAnalysis import NoiseAnalysis
from CenterCellAnalysis import CenterCellAnalysis
from NegativeChargesAnalysis import NegativeChargesAnalysis
from SaturationAnalysis import SaturationAnalysis
from FinalAnalysis import FinalAnalysis
from optparse import OptionParser
from Utils import *
import time

color_index = 10000

class TestAreas:
	def __init__(self, configfile='', run=0):
		self.do_fit = False
		self.do_saturation = True
		self.num = 0
		self.clust_size = 2
		self.dir = '.'
		self.run = run
		self.col_pitch = 50
		self.cell_resolution = 0
		self.phmin = 10000
		self.phmax = -10000
		self.binsperx = 0
		self.binspery = 0
		self.threshold = 0
		self.skip_before_sat = 0
		self.skip_after_sat = 1
		self.do_threshold = False
		self.neg_cut_adc = 4100
		self.neg_cut_snr = 410
		self.efficiency_subdiv = 1
		self.window_shift = 3
		self.min_snr_neg, self.max_snr_neg, self.delta_snr = -65, 1, 2
		self.min_snr, self.max_snr = -650, 650
		self.min_adc, self.max_adc = -6500, 6500
		self.delta_adc = 0
		self.min_adc_noise, self.max_adc_noise, self.delta_adc_noise = -322.5, 322.5, 0.5
		self.min_snr_noise, self.max_snr_noise, self.delta_snr_noise = -32.25, 32.25, 0.05
		self.minz = {t: {} for t in ['all', 'good', 'bad']}
		self.maxz = {t: {} for t in ['all', 'good', 'bad']}
		self.neg_cut_lines = {}
		self.trash = []
		self.w = 0
		self.num_rows_even = 0
		self.num_rows_odd = 0
		self.rows_pitch = 0
		self.cells_width = 0
		self.saturated_ADC = 0
		self.num_strips = 0
		self.cluster_size = 0
		self.conv_steps = 0
		self.sigma_conv = 0
		self.mpshift = 0
		self.num_parallel = 0
		self.hit_factor = 0
		self.seed_factor = 0
		self.num_sides = 0
		self.rows = []
		self.cols = []
		self.cells = []
		self.rcells = []
		self.config_file = configfile
		if self.config_file != '':
			self.config_file = Correct_Path(self.config_file)
			self.ReadConfigFile()
		self.trans_grid = TransparentGrid(self.dir, self.run, self.col_pitch)
		self.cluster_ch_ana = None
		self.noise_ana = None
		self.center_cells_ana = None
		self.neg_ana = None
		self.sat_ana = None
		self.final_ana = None
		self.bias = 0
		self.suffix = {'all': 'all', 'good': 'selection', 'bad': 'not_selection'}

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
				if pars.has_option('SETTINGS', 'col_pitch'):
					self.col_pitch = pars.getint('SETTINGS', 'col_pitch')
				if pars.has_option('SETTINGS', 'cluster_size'):
					self.cluster_size = pars.getint('SETTINGS', 'cluster_size')
				if pars.has_option('SETTINGS', 'num_strips'):
					self.num_strips = pars.getint('SETTINGS', 'num_strips')
				if pars.has_option('SETTINGS', 'test_number'):
					self.num = pars.getint('SETTINGS', 'test_number')
				if pars.has_option('SETTINGS', 'threshold'):
					self.threshold = pars.getfloat('SETTINGS', 'threshold')
				if pars.has_option('SETTINGS', 'do_threshold'):
					self.do_threshold = pars.getboolean('SETTINGS', 'do_threshold')
				if pars.has_option('SETTINGS', 'saturated_ADC'):
					self.saturated_ADC = pars.getint('SETTINGS', 'saturated_ADC')
				if pars.has_option('SETTINGS', 'num_sides'):
					self.num_sides = pars.getint('SETTINGS', 'num_sides')
				if pars.has_option('SETTINGS', 'neg_cut_adc'):
					self.neg_cut_adc = pars.getint('SETTINGS', 'neg_cut_adc')
				if pars.has_option('SETTINGS', 'neg_cut_snr'):
					self.neg_cut_snr = pars.getint('SETTINGS', 'neg_cut_snr')
				if pars.has_option('SETTINGS', 'skip_before_sat'):
					self.skip_before_sat = pars.getint('SETTINGS', 'skip_before_sat')
				if pars.has_option('SETTINGS', 'skip_after_sat'):
					self.skip_after_sat = pars.getint('SETTINGS', 'skip_after_sat')
				if pars.has_option('SETTINGS', 'hit_factor'):
					self.hit_factor = pars.getint('SETTINGS', 'hit_factor')
				if pars.has_option('SETTINGS', 'seed_factor'):
					self.seed_factor = pars.getint('SETTINGS', 'seed_factor')
				if pars.has_option('SETTINGS', 'cell_resolution'):
					self.cell_resolution = pars.getfloat('SETTINGS', 'cell_resolution')
				if pars.has_option('SETTINGS', 'delta_adc'):
					self.delta_adc = pars.getfloat('SETTINGS', 'delta_adc')
				if pars.has_option('SETTINGS', 'phmin'):
					self.phmin = pars.getfloat('SETTINGS', 'phmin')
				if pars.has_option('SETTINGS', 'phmax'):
					self.phmax = pars.getfloat('SETTINGS', 'phmax')
				if pars.has_option('SETTINGS', 'binsperx'):
					self.binsperx = pars.getfloat('SETTINGS', 'binsperx')
				if pars.has_option('SETTINGS', 'binspery'):
					self.binspery = pars.getfloat('SETTINGS', 'binspery')
				if pars.has_option('SETTINGS', 'efficiency_subdiv'):
					self.efficiency_subdiv = pars.getfloat('SETTINGS', 'efficiency_subdiv')
				if pars.has_option('SETTINGS', 'ch_ini'):
					self.ch_ini = pars.getint('SETTINGS', 'ch_ini')
				if pars.has_option('SETTINGS', 'ch_end'):
					self.ch_end = pars.getint('SETTINGS', 'ch_end')
				if pars.has_option('SETTINGS', 'num_parallel'):
					self.num_parallel = pars.getint('SETTINGS', 'num_parallel')
				if pars.has_option('SETTINGS', 'conv_steps'):
					self.conv_steps = pars.getint('SETTINGS', 'conv_steps')
				if pars.has_option('SETTINGS', 'sigma_conv'):
					self.sigma_conv = pars.getfloat('SETTINGS', 'sigma_conv')
				if pars.has_option('SETTINGS', 'mpshift'):
					self.mpshift = pars.getfloat('SETTINGS', 'mpshift')

			if pars.has_section('ROWS'):
				if pars.has_option('ROWS', 'rows'):
					rows = pars.get('ROWS', 'rows')
					self.rows = unpack_row_col(rows)
				if pars.has_option('ROWS', 'num'):
					self.num_rows_even = pars.getint('ROWS', 'num')
					self.num_rows_odd = pars.getint('ROWS', 'num')
				if pars.has_option('ROWS', 'num_even'):
					self.num_rows_even = pars.getint('ROWS', 'num_even')
				if pars.has_option('ROWS', 'num_odd'):
					self.num_rows_odd = pars.getint('ROWS', 'num_odd')
				if pars.has_option('ROWS', 'height'):
					self.rows_pitch = pars.getfloat('ROWS', 'height')
				if pars.has_option('ROWS', 'width'):
					self.cells_width = pars.getfloat('ROWS', 'width')
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
		"""
		This method selects the cells to take into account from the region from the selection analysis.
		This is done either by giving a detailed list of coloumns, rows or cells in the config file,
		or by giving a threshold and setting do_threshold as True to automatically select cells above the threshold
		:return: Nothing
		"""
		if self.do_threshold:
			if self.threshold == 0:
				self.trans_grid.FindThresholdCutFromCells('clusterChargeN', 'adc', self.phmin, self.phmax, self.delta_adc / 2.)
				self.threshold = self.trans_grid.threshold
			self.trans_grid.ResetAreas()
			print 'Selecting areas with a ph2 greater or equal than', self.threshold, '...', ; sys.stdout.flush()
			self.trans_grid.SelectGoodAndBadByThreshold(self.threshold, 'clusterChargeN')
			print 'Done'
			self.trans_grid.AddRemainingToBadAreas()
			print 'Marked the remaining cells as bad'
			if len(self.trans_grid.gridAreas.goodAreas_diamond) < 2:
				print 'There is only', len(self.trans_grid.gridAreas.goodAreas_diamond), 'cell in the selection. Check the thresholds or the area config file'
				if self.do_threshold:
					self.threshold = int(RoundInt(self.threshold * 0.75))
					print 'Trying with a new threshold of', self.threshold
					self.SetTest()
				else:
					return
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
		self.SetCutsInCutManager()
		# Update cell resolution value for analysis
		if self.cell_resolution == 0:
			if self.trans_grid.cell_resolution == 0:
				self.trans_grid.FindBinningAndResolution()
				self.cell_resolution = self.trans_grid.FindBinningAndResolution()
				self.trans_grid.SavePickle()
			else:
				self.cell_resolution = self.trans_grid.cell_resolution
		else:
			self.trans_grid.cell_resolution = self.cell_resolution
			self.trans_grid.SavePickle()
		self.noise_ana = NoiseAnalysis(self.trans_grid, self.num_strips, self.cluster_size)
		self.cluster_ch_ana = ClusterChannelsAnalysis(self.trans_grid, self.num_strips, self.cluster_size, self.noise_ana)
		self.trans_grid.FindMaxMinVarz()
		self.neg_ana = NegativeChargesAnalysis(self.trans_grid, self.num_strips, self.cluster_size, self.noise_ana)
		self.sat_ana = SaturationAnalysis(self.trans_grid, self.num_strips, self.cluster_size, self.noise_ana)
		self.center_cells_ana = CenterCellAnalysis(self.trans_grid, self.num_strips, self.cluster_size)
		self.final_ana = FinalAnalysis(self.trans_grid, self.num_strips, self.cluster_size, self.noise_ana, self.center_cells_ana)

	def SetCutsInCutManager(self):
		print 'Setting cuts in cut manager...', ; sys.stdout.flush()
		self.trans_grid.cuts_man.SetChs(self.trans_grid.gridAreas.good_channels, self.trans_grid.gridAreas.bad_channels)
		self.trans_grid.cuts_man.SetCells(selection=self.trans_grid.gridAreas.goodAreasCutNames_simplified_diamond, not_selection=self.trans_grid.gridAreas.badAreasCutNames_simplified_diamond)
		self.trans_grid.cuts_man.SetNoiseCuts()
		# self.trans_grid.cuts_man.SetPHCuts()
		print 'Done'

	def PositionCanvas(self, canvas_name):
		if canvas_name in self.trans_grid.canvas.keys():
			self.trans_grid.canvas[canvas_name].SetWindowPosition(self.w, self.w)
			ro.gPad.Update()
			self.w += self.window_shift

	def DoBorderPlots(self):
		"""
		This Method creates a 2D Profile Map of all the transparent events showing the cells divistion and which cells will be used for the analysis.
		Then it makes a projection in the Y axis and shows the fit that it was used to determine the edges Y limits.
		:return: Nothing
		"""
		self.trans_grid.DrawProfile2DDiamond('PH2_H_map_with_borders', self.trans_grid.GetPHNChsVar(2, 'H', False))
		self.trans_grid.DrawTCutGs('PH2_H_map_with_borders', 'diamond')
		self.trans_grid.DrawGoodAreasDiamondCenters('PH2_H_map_with_borders')
		# self.trans_grid.DrawGoodAreasDiamond('PH2_H_map_with_borders')
		# self.trans_grid.DrawBadAreasDiamond('PH2_H_map_with_borders')
		xbinmin, xbinmax = int(self.trans_grid.profile['PH2_H_map_with_borders'].GetXaxis().FindBin(self.trans_grid.ch_ini - 0.5)), int(self.trans_grid.profile['PH2_H_map_with_borders'].GetXaxis().FindBin(self.trans_grid.ch_ini - 0.5) + self.trans_grid.num_cols * self.trans_grid.bins_per_ch_x - 1)
		self.trans_grid.profile['PH2_H_map_with_borders'].GetXaxis().SetRange(xbinmin - 1, xbinmax + 1)
		self.trans_grid.canvas['PH2_H_map_with_borders_py'] = ro.TCanvas('c_PH2_H_map_with_borders_py', 'c_PH2_H_map_with_borders_py', 1)
		self.trans_grid.canvas['PH2_H_map_with_borders_py'].cd()
		self.trans_grid.histo['PH2_H_map_with_borders_py'] = self.trans_grid.profile['PH2_H_map_with_borders'].ProjectionY('h_PH2_H_map_with_borders_py', xbinmin, xbinmax, 'e hist')
		minbiny, maxbiny = self.trans_grid.histo['PH2_H_map_with_borders_py'].FindFirstBinAbove(), self.trans_grid.histo['PH2_H_map_with_borders_py'].FindLastBinAbove()
		for biny in xrange(maxbiny, int(self.trans_grid.histo['PH2_H_map_with_borders_py'].GetXaxis().GetNbins())):
			if self.trans_grid.histo['PH2_H_map_with_borders_py'].GetBinContent(biny) != 0:
				maxbiny = biny
		miny, maxy = self.trans_grid.histo['PH2_H_map_with_borders_py'].GetXaxis().GetBinLowEdge(minbiny), self.trans_grid.histo['PH2_H_map_with_borders_py'].GetXaxis().GetBinLowEdge(maxbiny + 1)
		self.trans_grid.canvas['PH2_H_map_with_borders'].cd()
		self.trans_grid.profile['PH2_H_map_with_borders'].GetYaxis().SetRangeUser(miny, maxy)
		ro.gPad.Update()
		self.PositionCanvas('PH2_H_map_with_borders')
		self.trans_grid.canvas['PH2_H_map_with_borders_py'].cd()
		self.trans_grid.histo['PH2_H_map_with_borders_py'].GetXaxis().SetRangeUser(miny, maxy)
		func = ro.TF1('box_fcn', '[0]*(TMath::Erf((x-({l}))/[1])+1)/2-[2]*(TMath::Erf((x-{u})/[3])+1)/2+[4]'.format(l=np.array([self.trans_grid.row_cell_info_diamond['0_even'], self.trans_grid.row_cell_info_diamond['0_odd']]).mean(), u=np.array([self.trans_grid.row_cell_info_diamond['up_even'], self.trans_grid.row_cell_info_diamond['up_even']]).mean()), miny, maxy)
		func.SetNpx(10000)
		zmin, zmax = self.trans_grid.histo['PH2_H_map_with_borders_py'].GetMinimum(), self.trans_grid.histo['PH2_H_map_with_borders_py'].GetMaximum()
		y1bin, y2bin = self.trans_grid.histo['PH2_H_map_with_borders_py'].FindFirstBinAbove((zmin + zmax) / 2.0), self.trans_grid.histo['PH2_H_map_with_borders_py'].FindLastBinAbove((zmin + zmax) / 2.0) + 1
		y1, y2 = self.trans_grid.histo['PH2_H_map_with_borders_py'].GetXaxis().GetBinCenter(y1bin), self.trans_grid.histo['PH2_H_map_with_borders_py'].GetXaxis().GetBinCenter(y2bin)
		z0, z1, z2 = self.trans_grid.histo['PH2_H_map_with_borders_py'].GetBinContent(int((minbiny))), self.trans_grid.histo['PH2_H_map_with_borders_py'].GetBinContent(int((y1bin + y2bin) / 2.0)), self.trans_grid.histo['PH2_H_map_with_borders_py'].GetBinContent(int((maxbiny)))
		func.SetParLimits(0, abs(z1 - z0) / 10.0, 10.0 * abs(z1 - z0))
		func.SetParLimits(1, 0.1, 50)
		func.SetParLimits(2, abs(z1 - z2) / 10.0, 10.0 * abs(z1 - z2))
		func.SetParLimits(3, 0.1, 50)
		func.SetParLimits(4, -100.0 * abs(z0), 100 * abs(z0))
		params = np.array((abs(z1 - z0), 20, abs(z1 - z2), 20, z0), 'float64')
		func.SetParameters(params)
		self.trans_grid.fits['PH2_H_map_with_borders_py'] = self.trans_grid.histo['PH2_H_map_with_borders_py'].Fit('box_fcn', 'QIEBMSL', 'goff', self.trans_grid.histo['PH2_H_map_with_borders_py'].GetBinLowEdge(int((minbiny))), self.trans_grid.histo['PH2_H_map_with_borders_py'].GetBinLowEdge(int((maxbiny))))
		params = np.zeros(5, 'f8')
		func.GetParameters(params)
		func.SetParameters(params)
		self.trans_grid.fits['PH2_H_map_with_borders_py'] = self.trans_grid.histo['PH2_H_map_with_borders_py'].Fit('box_fcn', 'QIEBMSL', 'goff', self.trans_grid.histo['PH2_H_map_with_borders_py'].GetBinLowEdge(int((minbiny))), self.trans_grid.histo['PH2_H_map_with_borders_py'].GetBinLowEdge(int((maxbiny))))
		self.PositionCanvas('PH2_H_map_with_borders_py')
		ro.gPad.Update()

	def DoCellHistograms(self, typ='adc'):
		self.trans_grid.mean_ph_cell_dic[typ] = {}
		fact = 1. if typ == 'adc' else 10.
		minx, maxx, deltax = self.phmin / fact, self.phmax / fact, self.delta_adc / (2. * fact)
		self.trans_grid.FindThresholdCutFromCells(self.trans_grid.GetPHNChsVar(self.num_strips, 'H', typ == 'snr'), typ, minx, maxx, deltax)
		self.PositionCanvas('mean_ph_per_cell_{t}'.format(t=typ))
		# self.trans_grid.FindThresholdCutFromCells(self.trans_grid.GetPHNChsVar(self.num_strips, 'H', False), typ, self.phmin, self.phmax, self.delta_adc / 2.)
		# self.PositionCanvas('mean_ph_per_cell_adc')

	def DoNoiseStudiesDifferentBuffers(self, cells='good', typ='adc'):
		suffix = self.suffix[cells]
		divs, minbuffer, maxbuff = 3., 10, 2000
		tempbuff = np.unique(np.floor(np.power(10, np.divide(np.arange(divs), float(divs), dtype='float64'), dtype='float64') + 0.5).astype('int32'))
		tempbuff = np.append(np.append(np.append(np.append(np.append(tempbuff, 10*tempbuff), 100*tempbuff), 1000*tempbuff), 10000*tempbuff), 100000*tempbuff)
		buffers = np.extract(np.bitwise_and(minbuffer <= tempbuff, tempbuff <= maxbuff), tempbuff)
		cond_do_buffs = np.bitwise_not([self.trans_grid.CheckIfPedTreeFriendExists(buff) for buff in buffers])
		buffers_to_create = np.extract(cond_do_buffs, buffers)
		self.trans_grid.CreatePedTreeFriendsForStudy(buffers_to_create)
		for buff in buffers:
			self.trans_grid.AddFriendWithNewPedestalBuffer(buff)
			self.DoNoiseFriendStudies(cells, typ)
			self.trans_grid.UnfriendTree(self.trans_grid.trans_tree.GetFriend('pedTree'))

		temp0 = np.array([buffers[1] - buffers[0]] + [buffers[i] - buffers[i-1] for i in xrange(1, buffers.size)])
		temp0fact = (buffers[-1] - buffers[0]) / (2 * temp0.sum(dtype='float64') - temp0[0] - temp0[-1])
		xarrayerrs = np.multiply(temp0, temp0fact, dtype='float64')

		ymeanarray = np.array([self.trans_grid.histo['signal_noise_buffer_{b}_{s}_adc'.format(s=suffix, b=buff)].GetRMS() for buff in buffers], dtype='float64')
		ymeanarrayerrs = np.array([self.trans_grid.histo['signal_noise_buffer_{b}_{s}_adc'.format(s=suffix, b=buff)].GetRMSError() for buff in buffers], dtype='float64')
		# maxy = yarray.max() + yarrayerrs.max()

		tgraph = ro.TGraphErrors(int(buffers.size), buffers.astype('float64'), ymeanarray, xarrayerrs, ymeanarrayerrs)
		graphname = 'signal_noise_Vs_buffer_sizes_in_adc_{s}'.format(s=suffix)
		tgraph.SetNameTitle('g_' + graphname, 'g_' + graphname)
		tgraph.GetXaxis().SetTitle('buffer size')
		tgraph.GetYaxis().SetTitle('signal noise [ADC]')
		# tgraph.GetYaxis().SetRangeUser(0, maxy)
		tgraph.SetMarkerStyle(8)
		tgraph.SetMarkerColor(ro.kBlack)
		tgraph.SetLineColor(ro.kBlack)
		self.trans_grid.graph[graphname] = tgraph
		self.trans_grid.canvas[graphname] = ro.TCanvas('c_' + graphname, 'c_' + graphname, 1)
		self.trans_grid.graph[graphname].Draw('ALP')
		SetDefault1DCanvasSettings(self.trans_grid.canvas[graphname])
		self.trans_grid.canvas[graphname].SetLogx()
		self.PositionCanvas(graphname)

		ymeanncarray = np.array([self.trans_grid.histo['signal_noise_NC_chs_buffer_{b}_{s}_adc'.format(b=buff, s='all')].GetRMS() for buff in buffers], dtype='float64')
		ymeanncarrayerrs = np.array([self.trans_grid.histo['signal_noise_NC_chs_buffer_{b}_{s}_adc'.format(b=buff, s='all')].GetRMSError() for buff in buffers], dtype='float64')

		tgraphnc = ro.TGraphErrors(int(buffers.size), buffers.astype('float64'), ymeanncarray, xarrayerrs, ymeanncarrayerrs)
		graphname = 'signal_noise_NC_chs_Vs_buffer_sizes_in_adc_{s}'.format(s=suffix)
		tgraphnc.SetNameTitle('g_' + graphname, 'g_' + graphname)
		tgraphnc.GetXaxis().SetTitle('buffer size')
		tgraphnc.GetYaxis().SetTitle('signal noise [ADC]')
		# tgraphnc.GetYaxis().SetRangeUser(0, maxy)
		tgraphnc.SetMarkerStyle(8)
		tgraphnc.SetMarkerColor(ro.kBlack)
		tgraphnc.SetLineColor(ro.kBlack)
		self.trans_grid.graph[graphname] = tgraphnc
		self.trans_grid.canvas[graphname] = ro.TCanvas('c_' + graphname, 'c_' + graphname, 1)
		self.trans_grid.graph[graphname].Draw('ALP')
		SetDefault1DCanvasSettings(self.trans_grid.canvas[graphname])
		self.trans_grid.canvas[graphname].SetLogx()
		self.PositionCanvas(graphname)

	def DoNoiseStudies(self, cells='all', typ='adc', isFriend=False):
		self.noise_ana.w, self.noise_ana.window_shift = self.w, self.window_shift
		self.noise_ana.DoNoiseAnalysis(cells, typ, isFriend)
		self.w, self.window_shift = self.noise_ana.w, self.noise_ana.window_shift

	def DoNoiseFriendStudies(self, cells='all', typ='adc', isFriend=False):
		self.noise_ana.w, self.noise_ana.window_shift = self.w, self.window_shift
		self.noise_ana.DoFriendNoiseAnalysis(cells, typ, isFriend)
		self.w, self.window_shift = self.noise_ana.w, self.noise_ana.window_shift

	def DoClusterStudies(self, cells='all', typ='adc', isFriend=False):
		self.cluster_ch_ana.w, self.cluster_ch_ana.window_shift = self.w, self.window_shift
		self.cluster_ch_ana.DoClusterStudies(cells, typ, isFriend)
		self.w, self.window_shift = self.cluster_ch_ana.w, self.cluster_ch_ana.window_shift

	def DoNegativeEventsStudies(self, cells='all', typ='adc', isFriend=False):
		self.neg_ana.w, self.neg_ana.window_shift = self.w, self.window_shift
		self.neg_ana.DoNegativeAnalysis(cells, typ, isFriend)
		self.w, self.window_shift = self.neg_ana.w, self.neg_ana.window_shift

	def DoSaturationStudies(self, cells='all', typ='adc', isFriend=False):
		self.sat_ana.w, self.sat_ana.window_shift = self.w, self.window_shift
		self.sat_ana.DoSaturationAnalysis(cells, self.skip_before_sat, self.skip_after_sat, typ=typ, isFriend=isFriend)
		self.w, self.window_shift = self.sat_ana.w, self.sat_ana.window_shift

	def DoFinalStudies(self, typ='adc', cummulative_chs=[2], isFriend=False):
		self.final_ana.w, self.final_ana.window_shift = self.w, self.window_shift
		self.final_ana.DoFinalAnalysis(typ=typ, cummulative_chs=cummulative_chs, isFriend=isFriend)
		self.w, self.window_shift = self.final_ana.w, self.final_ana.window_shift

	def DoCenterCellStudies(self, cells='all'):
		suffix = self.suffix[cells]
		self.center_cells_ana.w, self.center_cells_ana.window_shift = self.w, self.window_shift
		dists = np.arange(0.1, 1, 0.025)
		self.center_cells_ana.DoCenterRegionStudies(dists, cells, suffix)
		self.w, self.window_shift = self.center_cells_ana.w, self.center_cells_ana.window_shift

	def DoCenterCellSaturationStudies(self, cells='all'):
		self.trans_grid.AddFriendWithSaturationRegions(skipAfter=1, skipBefore=0)
		suffix = self.suffix[cells]
		dists = np.arange(0.1, 1, 0.025)
		percents = np.unique(np.floor(np.power(dists, 2) * 100 + 0.5).astype('int32'))
		for percent in percents:
			self.trans_grid.CreateTCutGSymmetricRectangle(percent)
		self.center_cells_ana.GetCutsFromCutManager(cells)
		self.center_cells_ana.GetVarzFromTranspGrid()

		dists = np.arange(0.1, 1, 0.025)
		percents = np.unique(np.floor(np.power(dists, 2) * 100 + 0.5).astype('int32'))

		xarray = percents.astype('float64')
		temp0 = np.array([percents[1] - percents[0]] + [percents[i] - percents[i-1] for i in xrange(1, percents.size)])
		temp0fact = (percents[-1] - percents[0]) / (2 * temp0.sum(dtype='float64') - temp0[0] - temp0[-1])
		xarrayerrs = np.multiply(temp0, temp0fact, dtype='float64')
		# ipdb.set_trace()
		ysatin = [self.trans_grid.trans_tree.Draw('event', self.center_cells_ana.sat_adc_inside_cut[p], 'goff') for p in percents]
		ysatout = [self.trans_grid.trans_tree.Draw('event', self.center_cells_ana.sat_adc_outside_cut[p], 'goff') for p in percents]
		ynosatin = [self.trans_grid.trans_tree.Draw('event', self.center_cells_ana.nosat_adc_inside_cut[p], 'goff') for p in percents]
		ynosatout = [self.trans_grid.trans_tree.Draw('event', self.center_cells_ana.nosat_adc_outside_cut[p], 'goff') for p in percents]
		ysat_nosat_ratio_in = np.divide(ysatin, np.add(ynosatin, ysatin), dtype='float64')
		ysat_nosat_ratio_out = np.divide(ysatout, np.add(ysatout, ynosatout), dtype='float64')
		ysat_nosat_ratio_in_errs = np.divide(np.sqrt(np.add(np.power(np.multiply(ysatin, np.sqrt(ynosatin)), 2), np.power(np.multiply(ynosatin, np.sqrt(ysatin)), 2))), np.power(np.add(ysatin, ynosatin), 2), dtype='float64')
		ysat_nosat_ratio_out_errs = np.divide(np.sqrt(np.add(np.power(np.multiply(ysatout, np.sqrt(ynosatout)), 2), np.power(np.multiply(ynosatout, np.sqrt(ysatout)), 2))), np.power(np.add(ysatout, ynosatout), 2), dtype='float64')

		tgraphe_sat_ratio_in = ro.TGraphErrors(int(xarray.size), xarray, ysat_nosat_ratio_in, xarrayerrs, ysat_nosat_ratio_in_errs)
		ingraphrationame = 'Sat_events_ratio_in_rect_Vs_percent_area_in'
		tgraphe_sat_ratio_in.SetNameTitle('g_' + ingraphrationame, 'g_' + ingraphrationame)
		tgraphe_sat_ratio_in.GetXaxis().SetTitle('percentage of rectangular area inside')
		tgraphe_sat_ratio_in.GetYaxis().SetTitle('sat_tracks/all_tracks')
		tgraphe_sat_ratio_in.GetYaxis().SetRangeUser(0, 1)
		tgraphe_sat_ratio_in.SetMarkerStyle(8)
		tgraphe_sat_ratio_in.SetMarkerColor(ro.kRed)
		tgraphe_sat_ratio_in.SetLineColor(ro.kRed)
		self.trans_grid.graph[ingraphrationame] = tgraphe_sat_ratio_in
		self.trans_grid.canvas[ingraphrationame] = ro.TCanvas('c_' + ingraphrationame, 'c_' + ingraphrationame, 1)
		self.trans_grid.graph[ingraphrationame].Draw('ALP')
		SetDefault1DCanvasSettings(self.trans_grid.canvas[ingraphrationame])
		self.PositionCanvas(ingraphrationame)

		tgraphe_sat_ratio_out = ro.TGraphErrors(int(xarray.size), xarray, ysat_nosat_ratio_out, xarrayerrs, ysat_nosat_ratio_out_errs)
		outgraphrationame = 'Sat_events_ratio_out_rect_Vs_percent_area_in'
		tgraphe_sat_ratio_out.SetNameTitle('g_' + outgraphrationame, 'g_' + outgraphrationame)
		tgraphe_sat_ratio_out.GetXaxis().SetTitle('percentage of rectangular area inside')
		tgraphe_sat_ratio_out.GetYaxis().SetTitle('sat_tracks/all_tracks')
		tgraphe_sat_ratio_out.GetYaxis().SetRangeUser(0, 1)
		tgraphe_sat_ratio_out.SetMarkerStyle(8)
		tgraphe_sat_ratio_out.SetMarkerColor(ro.kBlue)
		tgraphe_sat_ratio_out.SetLineColor(ro.kBlue)
		self.trans_grid.graph[outgraphrationame] = tgraphe_sat_ratio_out
		self.trans_grid.canvas[outgraphrationame] = ro.TCanvas('c_' + outgraphrationame, 'c_' + outgraphrationame, 1)
		self.trans_grid.graph[outgraphrationame].Draw('ALP')
		SetDefault1DCanvasSettings(self.trans_grid.canvas[outgraphrationame])
		self.PositionCanvas(outgraphrationame)

	def SaveCanvas(self):
		self.trans_grid.SaveCanvasInlist(self.trans_grid.canvas.keys())

	def DoAutomatic(self, cells='good', types=['adc'], isFriend=False, SaveAllPlots=False, isFirst=False):
		"""
		Makes a complete analysis
		:param cells: Specify which cells to analyse. Options are: "good": selected cells, "bad": not selected cells, "all": all the cells in the transparent grid
		:param types: Is a list that contains the type of variable to analyse. It can be ['adc'] (only adc), ['snr'] (only snr) or ['adc', 'snr'] (both adc and snr)
		:param isFriend: This flag is to do the analysis with a different pedestal calculation whose data is stored in a pedTree which should be a friend of the transparentTree. The method AddFriendWithNewPedestalBuffer should have been used before.
		:param SaveAllPlots: If true, all the plots created and which have a corresponding Canvas, will be stored as png and root files
		:param isFirst: If it is true, then the analysis will stop after DoClusterStudies, so that the user can select accurate values for neg_cut_adc and neg_cut_snr and then save the pickle.
		:return: Nothing
		"""

		self.window_shift = 1
		self.DoBorderPlots()
		for typ in ['adc', 'snr']:
			if typ in types:
				self.DoCellHistograms(typ)
				self.DoNoiseStudies(cells, typ, isFriend)
				self.DoClusterStudies(cells, typ, isFriend)
				if not isFirst:
					self.DoNegativeEventsStudies(cells, typ, isFriend)
					self.DoSaturationStudies(cells, typ, isFriend)
					self.DoFinalStudies(typ, cummulative_chs=[self.num_strips], isFriend=isFriend)
				# self.DoCenterCellStudies(cells)
				# self.DoCenterCellSaturationStudies(cells)
		# self.PlotTestClusterStudies(cells)
		# self.PlotTestForNegative(cells)
		# self.PlotTest()
		# if self.do_saturation:
		# 	self.PlotSaturation()
		if SaveAllPlots and not isFirst:
			self.SaveCanvas()

	def SetTransparentGrid(self, isFirst=False):
		"""
		This method tries to load a saved pickle from the test subdirectory. If it has not been created,
		it will create it with info from the config file and also by asking the user. Then it will update certain
		parameters that are set by the user in the config file. At the end it save the pickle file to take into account
		all the modifications
		:return: Nothing
		"""
		self.trans_grid.pkl_sbdir = 'test' + str(self.num)
		self.trans_grid.LoadPickle()

		is_first_time = False
		if not self.trans_grid.loaded_pickle or isFirst:
			is_first_time = True
			print 'It is first time the analysis runs on this data set test' + str(self.num) + '. Creating pickle file:'
			self.trans_grid.phmax = self.trans_grid.phmax if self.phmax == -10000 else self.phmax
			self.trans_grid.phmin = self.trans_grid.phmin if self.phmin == 10000 else self.phmin
			if self.delta_adc != 0:
				self.trans_grid.phbins = RoundInt(float(self.phmax - self.phmin) / self.delta_adc)
			self.trans_grid.row_cell_info_diamond['num_even'] = self.trans_grid.row_cell_info_diamond['num_even'] if self.num_rows_even == 0 else self.num_rows_even
			self.trans_grid.row_cell_info_diamond['num_odd'] = self.trans_grid.row_cell_info_diamond['num_odd'] if self.num_rows_odd == 0 else self.num_rows_odd
			self.trans_grid.row_cell_info_diamond['height'] = self.trans_grid.row_cell_info_diamond['height'] if self.rows_pitch == 0 else self.rows_pitch
			self.trans_grid.row_cell_info_diamond['width'] = self.trans_grid.row_cell_info_diamond['width'] if self.cells_width == 0 else self.cells_width

			self.trans_grid.SetupCutManager()
			if self.trans_grid.row_cell_info_diamond['0_even'] == 0 and self.trans_grid.row_cell_info_diamond['0_odd'] == 0:
				self.trans_grid.FindPickleValues(True, True)

			self.trans_grid.cluster_size = Get_From_User_Value('cluster size for run ' + str(self.run), 'int', self.trans_grid.cluster_size, True) if self.cluster_size == 0 else self.cluster_size
			self.trans_grid.num_strips = Get_From_User_Value('number of strips for run ' + str(self.run), 'int', self.trans_grid.num_strips, True) if self.num_strips == 0 else self.num_strips
			self.trans_grid.threshold = self.trans_grid.threshold if self.threshold == 0 else self.threshold
			self.trans_grid.cell_resolution = self.trans_grid.cell_resolution if self.cell_resolution == 0 else self.cell_resolution
			self.trans_grid.bins_per_ch_x = self.trans_grid.bins_per_ch_x if self.binsperx == 0 else self.binsperx
			self.trans_grid.bins_per_ch_y = self.trans_grid.bins_per_ch_y if self.binspery == 0 else self.binspery
			self.trans_grid.efficiency_subdiv = self.trans_grid.efficiency_subdiv if self.efficiency_subdiv == 1 else self.efficiency_subdiv
			self.trans_grid.saturated_ADC = Get_From_User_Value('saturated_ADC for run ' + str(self.run), 'int', self.trans_grid.saturated_ADC, True) if self.saturated_ADC == 0 else self.saturated_ADC
			self.trans_grid.bias = Get_From_User_Value('bias for run ' + str(self.run), 'float', self.trans_grid.bias, update=True)
			self.trans_grid.neg_cut_adc = self.neg_cut_adc if self.neg_cut_adc != 4100 else self.trans_grid.neg_cut_adc
			self.trans_grid.neg_cut_snr = self.neg_cut_snr if self.neg_cut_snr != 410 else self.trans_grid.neg_cut_snr
			self.trans_grid.conv_steps = self.conv_steps if self.conv_steps != 0 else self.trans_grid.conv_steps
			self.trans_grid.sigma_conv = self.trans_grid.sigma_conv if self.sigma_conv == 0 else self.sigma_conv
			self.trans_grid.mpshift = self.trans_grid.mpshift if self.mpshift == 0 else self.mpshift
			self.trans_grid.num_parallel = self.trans_grid.num_parallel if self.num_parallel == 0 else self.num_parallel
			self.trans_grid.hit_factor = self.trans_grid.hit_factor if self.hit_factor == 0 else self.hit_factor
			self.trans_grid.seed_factor = self.trans_grid.seed_factor if self.seed_factor == 0 else self.seed_factor
			self.trans_grid.num_sides = self.trans_grid.num_sides if self.num_sides == 0 else self.num_sides
			self.trans_grid.SavePickle()

			del self.trans_grid.cuts_man
			self.trans_grid.cuts_man = None

			print 'It is suggested to run the FindPickleValues method in transparent_grid to find the remaining parameters'

		print 'Updating variables not saved in the pickle...', ; sys.stdout.flush()
		self.cluster_size = self.trans_grid.cluster_size if self.cluster_size == 0 else self.cluster_size
		self.threshold = self.trans_grid.threshold if self.threshold == 0 else self.threshold
		self.num_strips = self.trans_grid.num_strips if self.num_strips == 0 else self.num_strips
		self.phmax = self.trans_grid.phmax if self.phmax == -10000 else self.phmax
		self.phmin = self.trans_grid.phmin if self.phmin == 10000 else self.phmin
		self.delta_adc = float(self.phmax - self.phmin) / float(self.trans_grid.phbins) if self.delta_adc == 0 else self.delta_adc
		self.binsperx = self.trans_grid.bins_per_ch_x if self.binsperx == 0 else self.binsperx
		self.binxpery = self.trans_grid.bins_per_ch_y if self.binspery == 0 else self.binspery
		self.efficiency_subdiv = self.trans_grid.efficiency_subdiv if self.efficiency_subdiv == 1 else self.efficiency_subdiv
		self.saturated_ADC = self.trans_grid.saturated_ADC if self.saturated_ADC == 0 else self.saturated_ADC
		self.num_rows_even = self.trans_grid.row_cell_info_diamond['num_even'] if self.num_rows_even == 0 else self.num_rows_even
		self.num_rows_odd = self.trans_grid.row_cell_info_diamond['num_odd'] if self.num_rows_odd == 0 else self.num_rows_odd
		self.rows_pitch = self.trans_grid.row_cell_info_diamond['height'] if self.rows_pitch == 0 else self.rows_pitch
		self.cells_width = self.trans_grid.row_cell_info_diamond['width'] if self.cells_width == 0 else self.cells_width
		self.bias = self.trans_grid.bias if self.bias == 0 else self.bias
		self.neg_cut_adc = self.trans_grid.neg_cut_adc if self.neg_cut_adc == 4100 else self.neg_cut_adc
		self.neg_cut_snr = self.trans_grid.neg_cut_snr if self.neg_cut_snr == 410 else self.neg_cut_snr
		self.conv_steps = self.trans_grid.conv_steps if self.conv_steps == 0 else self.conv_steps
		self.sigma_conv = self.trans_grid.sigma_conv if self.sigma_conv == 0 else self.sigma_conv
		self.mpshift = self.trans_grid.mpshift if self.mpshift == 0 else self.mpshift
		self.num_parallel = self.trans_grid.num_parallel if self.num_parallel == 0 else self.num_parallel
		self.hit_factor = self.trans_grid.hit_factor if self.hit_factor == 0 else self.hit_factor
		self.seed_factor = self.trans_grid.seed_factor if self.seed_factor == 0 else self.seed_factor
		self.num_sides = self.trans_grid.num_sides if self.num_sides == 0 else self.num_sides
		print 'Done'

		self.trans_grid.cluster_size = self.cluster_size
		self.trans_grid.num_strips = self.num_strips
		self.trans_grid.threshold = self.threshold
		self.trans_grid.cell_resolution = self.cell_resolution
		self.trans_grid.phmax = self.phmax
		self.trans_grid.phmin = self.phmin
		self.trans_grid.phbins = RoundInt(float(self.phmax - self.phmin) / self.delta_adc)
		self.trans_grid.bins_per_ch_x = self.binsperx
		self.trans_grid.bins_per_ch_y = self.binspery
		self.trans_grid.efficiency_subdiv = self.efficiency_subdiv
		self.trans_grid.saturated_ADC = self.saturated_ADC
		self.trans_grid.row_cell_info_diamond['num_even'] = self.num_rows_even
		self.trans_grid.row_cell_info_diamond['num_odd'] = self.num_rows_odd
		self.trans_grid.row_cell_info_diamond['height'] = self.rows_pitch
		self.trans_grid.row_cell_info_diamond['width'] = self.cells_width
		self.trans_grid.bias = self.bias
		self.trans_grid.neg_cut_adc = self.neg_cut_adc
		self.trans_grid.neg_cut_snr = self.neg_cut_snr
		self.trans_grid.conv_steps = self.conv_steps
		self.trans_grid.sigma_conv = self.sigma_conv
		self.trans_grid.mpshift = self.mpshift
		self.trans_grid.num_parallel = self.num_parallel
		self.trans_grid.hit_factor = self.hit_factor
		self.trans_grid.seed_factor = self.seed_factor
		self.trans_grid.num_sides = self.num_sides
		self.trans_grid.SavePickle()

		if self.trans_grid.gridAreas:
			self.trans_grid.ResetAreas()
		self.trans_grid.CreateGridAreas()
		self.trans_grid.CreateTCutGs()
		self.trans_grid.AddFriendWithCells(self.trans_grid.row_cell_info_diamond['x_off'], self.trans_grid.row_cell_info_diamond['y_off'])

if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-r', '--run', dest='run', type='int', default=0, help='run number to be analysed (e.g. 25209)')
	parser.add_option('-a', '--auto', dest='auto', default=False, action='store_true', help='Sets up test, creates plots and saves them automatically if toggled')
	parser.add_option('-c', '--config', dest='config', default='', type='string', help='gives the path to a config file for the test area')
	parser.add_option('-t', '--type', dest='typ', default='[adc]', type='string', help='selects a type of analysis to run automatically. Options are [adc], [adc,snr], or [snr] order does not matter.')
	parser.add_option('-f', '--first', dest='first', default=False, action='store_true', help='Use when negative cuts are unknown. Will do analysis until cluster channel part.')

	(options, args) = parser.parse_args()
	run = int(options.run)
	autom = bool(options.auto)
	config = str(options.config)
	typ = str(options.typ).strip('[').strip(']').split(',')
	first = bool(options.first)

	t = TestAreas(config, run)
	t.SetTransparentGrid(first)
	t.SetTest()
	if autom:
		t.DoAutomatic('good', types=typ if not first else ['adc', 'snr'], SaveAllPlots=False, isFirst=first)

