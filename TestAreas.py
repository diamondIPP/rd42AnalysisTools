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
		self.cellsize = 50
		self.cell_resolution = 0
		self.phmin = 10000
		self.phmax = -10000
		self.binsperx = 0
		self.binspery = 0
		self.threshold = 800
		self.skip_before_sat = 5
		self.skip_after_sat = 95
		self.do_threshold = False
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
		self.num_rows = 0
		self.saturated_ADC = 0
		self.num_strips = 0
		self.cluster_size = 0
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
		self.neg_ana = None
		self.sat_ana = None
		self.final_ana = None
		self.bias = self.trans_grid.bias
		self.suffix = {'all': 'all', 'good': 'selection', 'bad': 'not_selection'}
		if self.num_rows != 0:
			self.trans_grid.row_info_diamond['num'] = self.num_rows

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
				if pars.has_option('SETTINGS', 'num_strips'):
					self.num_strips = pars.getint('SETTINGS', 'num_strips')
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
				if pars.has_option('SETTINGS', 'saturated_ADC'):
					self.saturated_ADC = pars.getint('SETTINGS', 'saturated_ADC')
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
			if self.threshold == 0:
				self.trans_grid.FindThresholdCutFromCells('clusterCharge2', 'adc', self.phmin, self.phmax, self.delta_adc / 2.)
				self.threshold = self.trans_grid.threshold
			self.trans_grid.ResetAreas()
			print 'Selecting areas with a ph2 greater or equal than', self.threshold, '...', ; sys.stdout.flush()
			self.trans_grid.SelectGoodAndBadByThreshold(self.threshold, 'clusterCharge2')
			print 'Done'
			self.trans_grid.AddRemainingToBadAreas()
			print 'Marked the remaining cells as bad'
			if len(self.trans_grid.gridAreas.goodAreas_diamond) < 2:
				print 'There is only', len(self.trans_grid.gridAreas.goodAreas_diamond), 'cell in the selection. Check the thresholds or the area config file. Break'
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
		func = ro.TF1('box_fcn', '[0]*(TMath::Erf((x-({l}))/[1])+1)/2-[2]*(TMath::Erf((x-{u})/[3])+1)/2+[4]'.format(l=self.trans_grid.row_info_diamond['0'], u=self.trans_grid.row_info_diamond['up']), miny, maxy)
		func.SetNpx(int(self.trans_grid.row_info_diamond['num'] * self.trans_grid.bins_per_ch_y * 100))
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

	def DoCellHistograms(self):
		self.trans_grid.mean_ph_cell_dic = {'adc': {}, 'snr': {}}
		self.trans_grid.FindThresholdCutFromCells(self.trans_grid.GetPHNChsVar(2, 'H', True), 'snr', self.phmin / 10., self.phmax / 10., self.delta_adc / 20.)
		self.PositionCanvas('mean_ph_per_cell_snr')
		self.trans_grid.FindThresholdCutFromCells(self.trans_grid.GetPHNChsVar(2, 'H', False), 'adc', self.phmin, self.phmax, self.delta_adc / 2.)
		self.PositionCanvas('mean_ph_per_cell_adc')

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
		self.window_shift = 1
		self.DoBorderPlots()
		self.DoCellHistograms()
		for typ in ['adc', 'snr']:
			if typ in types:
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

	def SetTransparentGrid(self):
		self.trans_grid.SetLines()
		self.trans_grid.CreateTCutGs()
		if self.saturated_ADC == 0:
			if self.trans_grid.saturated_ADC != 0:
				self.saturated_ADC = self.trans_grid.saturated_ADC
			else:
				self.saturated_ADC = Get_From_User_Value('saturated_ADC for run ' + str(self.run), 'int', self.trans_grid.saturated_ADC, True)
				self.trans_grid.saturated_ADC = self.saturated_ADC
				self.trans_grid.SavePickle()
		else:
			self.trans_grid.saturated_ADC = self.saturated_ADC
			self.trans_grid.SavePickle()
		if self.trans_grid.bias != 0:
			self.bias = self.trans_grid.bias
		else:
			self.bias = Get_From_User_Value('bias for run ' + str(self.run), 'float', self.trans_grid.bias, update=True)
			self.trans_grid.bias = self.bias
			self.trans_grid.SavePickle()
		if self.trans_grid.neg_cut_adc == 4100:
			if self.trans_grid.neg_cut_snr != 410:
				self.trans_grid.DrawHisto1D('NoiseAllCells', -32.25, 32.25, 0.5, 'diaChSignal', 'signal not in cluster cmc [ADC]', self.trans_grid.cuts_man.not_in_transp_cluster, option='goff')
				self.trans_grid.neg_cut_adc = self.trans_grid.neg_cut_snr * self.trans_grid.histo['NoiseAllCells'].GetRMS()
				print 'Setting neg_cut_adc to', self.trans_grid.neg_cut_adc, '. If you wish to change it, do it and save the pickle again in the transparent grid object'
				self.trans_grid.SavePickle()
			else:
				print 'Run the noise and the cluster analysis to determine the neg_cut_snr and neg_cut_adc values'

		self.trans_grid.FindMaxMinVarz()
		# Update cluster size:
		if self.cluster_size == 0:
			self.cluster_size = self.trans_grid.cluster_size
		else:
			self.trans_grid.cluster_size = self.cluster_size
		# Update num_strips
		if self.num_strips == 0:
			self.num_strips = self.trans_grid.num_strips
		else:
			self.trans_grid.num_strips = self.num_strips
		# Update phmax for value analysis
		if self.phmax != -10000:
			self.trans_grid.phmax = self.phmax
		else:
			self.phmax = self.trans_grid.phmax
		# Update phmin for value analysis
		if self.phmin != 10000:
			self.trans_grid.phmin = self.phmin
		else:
			self.phmin = self.trans_grid.phmin
		# update ph binning:
		if self.delta_adc != 0:
			self.trans_grid.phbins = RoundInt(float(self.phmax - self.phmin) / self.delta_adc)
		else:
			self.delta_adc = float(self.phmax - self.phmin) / float(self.trans_grid.phbins)
		# ph binx x update:
		if self.binsperx != 0:
			self.trans_grid.bins_per_ch_x = self.binsperx
		# ph bins y update:
		if self.binspery != 0:
			self.trans_grid.bins_per_ch_y = self.binspery
		# efficiency subdiv
		if self.efficiency_subdiv != 1:
			self.trans_grid.efficiency_subdiv = self.efficiency_subdiv



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
	t.SetTransparentGrid()
	t.SetTest()
	if t.trans_grid.loaded_pickle:
		t.trans_grid.LoadPickle()
	else:
		t.trans_grid.FindPickleValues(False)
		ExitMessage('Run it again to load the generated pickle :)', os.EX_OK)
	if t.trans_grid.loaded_default_pickle:
		t.trans_grid.FindBinningAndResolution()
		t.trans_grid.SavePickle()
	if autom:
		t.DoAutomatic('good', types=typ if not first else ['adc', 'snr'], SaveAllPlots=True, isFirst=first)

