#!/usr/bin/env python
import ROOT as ro
import numpy as np
import ipdb, os
from ConfigParser import ConfigParser
from TransparentGrid import TransparentGrid
from optparse import OptionParser
from NoiseAnalysis import NoiseAnalysis
from CenterCellAnalysis import CenterCellAnalysis
from Utils import *


color_index = 10000

class FinalAnalysis:
	def __init__(self, trans_grid, numstrips, clustersize, noise_ana=None, center_reg_ana=None):
		self.window_shift = 2
		self.min_snr_neg, self.max_snr_neg = -64, 0
		self.min_adc_neg, self.max_adc_neg = -650, 0
		self.delta_ev = 100
		self.min_snr, self.max_snr = -650, 650
		self.min_adc, self.max_adc = -6500, 6500
		self.delta_adc, self.delta_snr = 20, 2
		self.min_adc_noise, self.max_adc_noise, self.delta_adc_noise = -322.5, 322.5, 0.5
		self.min_snr_noise, self.max_snr_noise, self.delta_snr_noise = -32.25, 32.25, 0.05
		self.trash = []
		self.w = 0
		self.trans_grid = trans_grid
		self.num_strips = numstrips
		self.cluster_size = clustersize
		self.noise_ana = noise_ana
		self.center_reg_ana = center_reg_ana
		self.cell_dists = np.arange(0.1, 1, 0.025)

		self.suffix = {'all': 'all', 'good': 'selection', 'bad': 'not_selection'}

		self.noise_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.noise_nc_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.noise_friend_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.noise_nc_friend_cuts = {t: '' for t in ['all', 'good', 'bad']}

		self.ph_adc_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.ph_snr_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.ph_adc_h_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.ph_snr_h_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_adc_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_snr_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_adc_h_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_snr_h_cuts = {t: {} for t in ['all', 'good', 'bad']}

		self.minz = {t: {} for t in ['all', 'good', 'bad']}
		self.maxz = {t: {} for t in ['all', 'good', 'bad']}

		self.sat_adc_ch_cut = {}
		self.sat_adc_h_cut = {}
		self.not_sat_adc_ch_cut = {}
		self.not_sat_adc_h_cut = {}
		self.sat_adc_N_ch_cut = {}
		self.sat_adc_N_h_cut = {}
		self.not_sat_adc_N_ch_cut = {}
		self.not_sat_adc_N_h_cut = {}

		self.sat_evts_region = ''
		self.not_sat_evts_region = ''

		self.in_transp_cluster = ''

		self.noise_varz = {}
		self.noise_friend_varz = {}
		
		self.neg_snr_ph_ch = {}
		self.neg_snr_ph_h = {}
		self.not_neg_snr_ph_ch = {}
		self.not_neg_snr_ph_h = {}

		self.neg_adc_ph_ch = {}
		self.neg_adc_ph_h = {}
		self.not_neg_adc_ph_ch = {}
		self.not_neg_adc_ph_h = {}

		self.neg_snr_phN_ch = {}
		self.neg_snr_phN_h = {}
		self.not_neg_snr_phN_ch = {}
		self.not_neg_snr_phN_h = {}

		self.neg_adc_phN_ch = {}
		self.neg_adc_phN_h = {}
		self.not_neg_adc_phN_ch = {}
		self.not_neg_adc_phN_h = {}

		self.noise_varz = {}
		self.ph_adc_h_varz = {}
		self.ph_adc_ch_varz = {}
		self.ph_snr_h_varz = {}
		self.ph_snr_ch_varz = {}

		self.phN_adc_h_varz = {}
		self.phN_adc_ch_varz = {}
		self.phN_snr_h_varz = {}
		self.phN_snr_ch_varz = {}

		self.analysis_cummulative_ch = np.arange(1, self.cluster_size + 1)

	def PosCanvas(self, canvas_name):
		self.w = PositionCanvas(self.trans_grid, canvas_name, self.w, self.window_shift)

	def DoDeviceMaps(self, cells, cuts='', suffix='no_cuts', typ='adc'):
		def DrawProfile2D(name, varz, varzname, cut, getOccupancy=False):
			self.trans_grid.DrawProfile2DDiamondMap(name, varz, varzname, cells, cut, 'prof colz')
			self.trans_grid.profile[name].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile[name].GetYaxis().SetRangeUser(self.trans_grid.row_info_diamond['0'] - int(self.trans_grid.row_info_diamond['pitch'] / self.trans_grid.bins_per_ch_y), self.trans_grid.row_info_diamond['up'] + int(self.trans_grid.row_info_diamond['pitch'] / self.trans_grid.bins_per_ch_y))
			self.trans_grid.DrawTCutGs(name, 'diamond')
			self.trans_grid.DrawGoodAreasDiamondCenters(name)
			self.PosCanvas(name)
			if getOccupancy:
				self.trans_grid.GetOccupancyFromProfile(name)
				self.trans_grid.histo['hit_map_' + name].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
				self.trans_grid.histo['hit_map_' + name].GetYaxis().SetRangeUser(self.trans_grid.row_info_diamond['0'] - int(self.trans_grid.row_info_diamond['pitch'] / self.trans_grid.bins_per_ch_y), self.trans_grid.row_info_diamond['up'] + int(self.trans_grid.row_info_diamond['pitch'] / self.trans_grid.bins_per_ch_y))
				self.trans_grid.DrawTCutGs('hit_map_' + name, 'diamond')
				self.trans_grid.DrawGoodAreasDiamondCenters('hit_map_' + name)
				self.PosCanvas('hit_map_' + name)

		tempcsnr = self.not_neg_snr_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)] if cuts == 'no_neg' else self.trans_grid.cuts_man.ConcatenateCuts(self.not_neg_snr_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)], self.not_sat_evts_region) if cuts == 'no_neg_no_sat' else ''
		tempcadc = self.not_neg_adc_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)] if cuts == 'no_neg' else self.trans_grid.cuts_man.ConcatenateCuts(self.not_neg_adc_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)], self.not_sat_evts_region) if cuts == 'no_neg_no_sat' else ''
		for ch in self.analysis_cummulative_ch:
			if typ == 'adc':
				DrawProfile2D('PH{c}_Ch_map_adc_{s}'.format(c=ch, s=suffix), self.phN_adc_ch_varz['PH{c}_Ch'.format(c=ch)], 'PH{c} cluster chs [ADC]'.format(c=ch), tempcadc, False)
				if ch != self.cluster_size: DrawProfile2D('PH{c}_H_map_adc_{s}'.format(c=ch, s=suffix), self.phN_adc_h_varz['PH{c}_H'.format(c=ch)], 'PH{c} highest chs [ADC]'.format(c=ch), tempcadc, ch == self.cluster_size)
			else:
				DrawProfile2D('PH{c}_Ch_map_snr_{s}'.format(c=ch, s=suffix), self.phN_snr_ch_varz['PH{c}_Ch'.format(c=ch)], 'PH{c} cluster chs [SNR]'.format(c=ch), tempcsnr, False)
				if ch != self.cluster_size: DrawProfile2D('PH{c}_H_map_snr_{s}'.format(c=ch, s=suffix), self.phN_snr_h_varz['PH{c}_H'.format(c=ch)], 'PH{c} highest chs [SNR]'.format(c=ch), tempcsnr, ch == self.cluster_size)

	def DoCellMaps(self, cells, cuts='', suffix='no_cuts', typ='adc'):
		def PlotCellsProfiles(name, varz, zmin, zmax, varname, cut, doOccupancy=False):
			self.trans_grid.DrawProfile2DDiamondCellOverlay(name, varz, cells, cut, varname=varname)
			self.trans_grid.profile[name].SetMinimum(min(0, zmin))
			self.trans_grid.profile[name].SetMaximum(zmax)
			self.PosCanvas(name)
			if doOccupancy:
				self.trans_grid.GetOccupancyFromProfile(name)
				self.PosCanvas('hit_map_' + name)

		tempcsnr = self.not_neg_snr_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)] if cuts == 'no_neg' else self.trans_grid.cuts_man.ConcatenateCuts(self.not_neg_snr_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)], self.not_sat_evts_region) if cuts == 'no_neg_no_sat' else '(1)'
		tempcadc = self.not_neg_adc_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)] if cuts == 'no_neg' else self.trans_grid.cuts_man.ConcatenateCuts(self.not_neg_adc_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)], self.not_sat_evts_region) if cuts == 'no_neg_no_sat' else '(1)'
		for ch in self.analysis_cummulative_ch:
			if typ == 'adc':
				minz, maxz = self.trans_grid.minz[cells]['PH{c}_Ch_adc'.format(c=ch)], self.trans_grid.minz[cells]['PH{c}_Ch_adc'.format(c=ch)]
				PlotCellsProfiles('PH{c}_Ch_cell_map_adc_{s}'.format(c=ch, s=suffix), self.phN_adc_ch_varz['PH{c}_Ch'.format(c=ch)], minz, max(maxz, 4800), 'PH{c} cluster chs [ADC]', tempcadc)
				minz, maxz = self.trans_grid.minz[cells]['PH{c}_H_adc'.format(c=ch)], self.trans_grid.minz[cells]['PH{c}_H_adc'.format(c=ch)]
				if ch != self.cluster_size: PlotCellsProfiles('PH{c}_H_cell_map_adc_{s}'.format(c=ch, s=suffix), self.phN_adc_h_varz['PH{c}_H'.format(c=ch)], minz, max(maxz, 4800), 'PH{c} highest chs [ADC]', tempcadc, ch == self.cluster_size)
			else:
				minz, maxz = self.trans_grid.minz[cells]['PH{c}_Ch_snr'.format(c=ch)], self.trans_grid.minz[cells]['PH{c}_Ch_snr'.format(c=ch)]
				PlotCellsProfiles('PH{c}_Ch_cell_map_snr_{s}'.format(c=ch, s=suffix), self.phN_snr_ch_varz['PH{c}_Ch'.format(c=ch)], minz, max(maxz, 480), 'PH{c} cluster chs [SNR]', tempcsnr)
				minz, maxz = self.trans_grid.minz[cells]['PH{c}_H_snr'.format(c=ch)], self.trans_grid.minz[cells]['PH{c}_H_snr'.format(c=ch)]
				if ch != self.cluster_size: PlotCellsProfiles('PH{c}_H_cell_map_snr_{s}'.format(c=ch, s=suffix), self.phN_snr_h_varz['PH{c}_H'.format(c=ch)], minz, max(maxz, 480), 'PH{c} highest chs [SNR]', tempcsnr, ch == self.cluster_size)

	def DoStripHistograms(self, cells, cuts='', suffix='no_cuts', typ='adc'):
		minx, maxx, deltax, xname, xvar = -0.5, 0.5, self.trans_grid.cell_resolution / float(self.trans_grid.row_info_diamond['pitch']), 'dia pred. strip hit pos', 'diaChXPred-TMath::Floor(diaChXPred+0.5)'
		def Draw2DHistogram(name, zmin, zmax, yname, yvar, cuts, typ='adc'):
			deltay = 4 * self.delta_adc if typ == 'adc' else 4 * self.delta_snr
			histo_limits = Get1DLimits(zmin, zmax, deltay)
			self.trans_grid.DrawHisto2D(name, minx, maxx, deltax, xname, min(0, histo_limits['min']), histo_limits['max'], deltay, yname, xvar, yvar, cuts)
			self.PosCanvas(name)

		def Draw1DHistogram(name, cuts):
			self.trans_grid.DrawHisto1D(name, minx, maxx, deltax, xvar, xname, cuts)
			maxbin = self.trans_grid.histo[name].GetMaximumBin()
			self.trans_grid.histo[name].GetYaxis().SetRangeUser(0, self.trans_grid.histo[name].GetBinContent(maxbin) + self.trans_grid.histo[name].GetBinError(maxbin))
			self.PosCanvas(name)

		tempcsnr = self.not_neg_snr_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)] if cuts == 'no_neg' else self.trans_grid.cuts_man.ConcatenateCuts(self.not_neg_snr_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)], self.not_sat_evts_region) if cuts == 'no_neg_no_sat' else '(1)'
		tempcadc = self.not_neg_adc_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)] if cuts == 'no_neg' else self.trans_grid.cuts_man.ConcatenateCuts(self.not_neg_adc_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)], self.not_sat_evts_region) if cuts == 'no_neg_no_sat' else '(1)'
		tempcsnr = self.trans_grid.cuts_man.ConcatenateCutWithCells(tempcsnr, cells)
		tempcadc = self.trans_grid.cuts_man.ConcatenateCutWithCells(tempcadc, cells)

		for ch in self.analysis_cummulative_ch:
			if typ == 'adc':
				minzadc, maxzadc = min(0, self.trans_grid.minz[cells]['PH{c}_Ch_adc'.format(c=ch)]), self.trans_grid.maxz[cells]['PH{c}_Ch_adc'.format(c=ch)]
				Draw2DHistogram('PH{c}_Ch_Vs_strip_location_adc_{s}'.format(c=ch, s=suffix), minzadc, maxzadc, 'PH{c} cluster chs [ADC]', self.phN_adc_ch_varz['PH{c}_Ch'.format(c=ch)], tempcadc, 'adc')
				minzadc, maxzadc = min(0, self.trans_grid.minz[cells]['PH{c}_H_adc'.format(c=ch)]), self.trans_grid.maxz[cells]['PH{c}_H_adc'.format(c=ch)]
				if ch != self.cluster_size: Draw2DHistogram('PH{c}_H_Vs_strip_location_adc_{s}'.format(c=ch, s=suffix), minzadc, maxzadc, 'PH{c} highest chs [ADC]', self.phN_adc_h_varz['PH{c}_H'.format(c=ch)], tempcadc, 'adc')
			else:
				minzsnr, maxzsnr = min(0, self.trans_grid.minz[cells]['PH{c}_Ch_snr'.format(c=ch)]), self.trans_grid.maxz[cells]['PH{c}_Ch_snr'.format(c=ch)]
				Draw2DHistogram('PH{c}_Ch_Vs_strip_location_snr_{s}'.format(c=ch, s=suffix), minzsnr, maxzsnr, 'PH{c} cluster chs [SNR]', self.phN_snr_ch_varz['PH{c}_Ch'.format(c=ch)], tempcsnr, 'snr')
				minzsnr, maxzsnr = min(0, self.trans_grid.minz[cells]['PH{c}_H_snr'.format(c=ch)]), self.trans_grid.maxz[cells]['PH{c}_H_snr'.format(c=ch)]
				if ch != self.cluster_size: Draw2DHistogram('PH{c}_H_Vs_strip_location_snr_{s}'.format(c=ch, s=suffix), minzsnr, maxzsnr, 'PH{c} highest chs [SNR]', self.phN_snr_h_varz['PH{c}_H'.format(c=ch)], tempcsnr, 'snr')

		if typ == 'adc':
			Draw1DHistogram('strip_location_adc_{s}'.format(s=suffix), tempcadc)
		else:
			Draw1DHistogram('strip_location_snr_{s}'.format(s=suffix), tempcsnr)

	def DoEfficiencyPlots(self, cells, cuts='', suffix='no_cuts', typ='adc'):

		def DrawEfficiencyGraphs(name, var, cells, cuts, typ='adc', show_only_95=True):
			xmin, xmax, deltax = (0, 4800, 50) if typ == 'adc' else (0, 480, 5)
			ymin = 0.95 if show_only_95 else 0
			self.trans_grid.DrawEfficiencyGraph(name, var, cells, cuts, xmin, xmax, deltax, typ, ymin)
			self.PosCanvas(name)

		def DrawEfficiencyCellMaps(name, var, cells, cuts, ncuts=[1, 2, 5, 10, 20], typ='adc'):
			self.trans_grid.DrawProfile2DDiamondCellOverlay(name + '_h0_', var=var, cells=cells, cuts=cuts, plot_option='prof colz goff')
			self.trans_grid.GetOccupancyFromProfile(name + '_h0_', 'colz goff')
			typ_cut = 'sigma' if typ == 'snr' else 'adc'
			for ns in ncuts:
				tempcut = '({v}>={th})'.format(v=var, th=ns)
				cut = tempcut if cuts == '' else self.trans_grid.cuts_man.ConcatenateCuts(cuts, tempcut)
				self.trans_grid.DrawProfile2DDiamondCellOverlay(name + '_{n}{t}_cut'.format(n=ns, t=typ_cut), var=var, cells=cells, cuts=cut, plot_option='prof colz goff')
				self.trans_grid.GetOccupancyFromProfile(name + '_{n}{t}_cut'.format(n=ns, t=typ_cut), 'colz goff')
				self.trans_grid.histo[name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)] = self.trans_grid.histo['hit_map_' + name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)].Clone('h_' + name + '_{n}{t}_cut'.format(n=ns, t=typ_cut))
				self.trans_grid.histo[name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)].SetTitle('h_' + name + '_{n}{t}_cut'.format(n=ns, t=typ_cut))
				self.trans_grid.histo[name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)].Sumw2()
				self.trans_grid.histo[name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)].Divide(self.trans_grid.histo['hit_map_' + name + '_h0_'])
				self.trans_grid.canvas[name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)] = ro.TCanvas('c_' + name + '_{n}{t}_cut'.format(n=ns, t=typ_cut), 'c_' + name + '_{n}{t}_cut'.format(n=ns, t=typ_cut), 1)
				self.trans_grid.histo[name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)].GetZaxis().SetTitle('Efficiency')
				self.trans_grid.histo[name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)].GetZaxis().SetRangeUser(0.9, 1)
				self.trans_grid.histo[name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)].Draw('colz')
				SetDefault2DStats(self.trans_grid.histo[name + '_{n}{t}_cut'.format(n=ns, t=typ_cut)])
				self.PosCanvas(name + '_{n}{t}_cut'.format(n=ns, t=typ_cut))

		tempcsnr = self.not_neg_snr_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)] if cuts == 'no_neg' else self.trans_grid.cuts_man.ConcatenateCuts(self.not_neg_snr_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)], self.not_sat_evts_region) if cuts == 'no_neg_no_sat' else '(1)'
		tempcadc = self.not_neg_adc_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)] if cuts == 'no_neg' else self.trans_grid.cuts_man.ConcatenateCuts(self.not_neg_adc_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)], self.not_sat_evts_region) if cuts == 'no_neg_no_sat' else '(1)'
		for ch in self.analysis_cummulative_ch:
			if typ == 'adc':
				DrawEfficiencyCellMaps('PH{c}_H_Efficiency_adc_{s}'.format(c=ch, s=suffix), self.trans_grid.phN_adc_h_varz['PH{c}_H'.format(c=ch)], cells, tempcadc, ncuts=[10, 50, 100, 200, 500], typ=typ)
				DrawEfficiencyGraphs('Eff_PH{c}_H_Vs_Threshold_adc_{s}'.format(c=ch, s=suffix), self.trans_grid.phN_adc_h_varz['PH{c}_H'.format(c=ch)], cells, tempcadc, typ)
				# DrawEfficiencyGraphs('Eff_PH{c}_Ch_Vs_Threshold_adc_{s}'.format(c=ch, s=suffix), self.trans_grid.phN_adc_ch_varz['PH{c}_Ch'.format(c=ch)], cells, tempcadc, typ)
			else:
				DrawEfficiencyCellMaps('PH{c}_H_Efficiency_snr_{s}'.format(c=ch, s=suffix), self.trans_grid.phN_snr_h_varz['PH{c}_H'.format(c=ch)], cells, tempcsnr, ncuts=[1, 5, 10, 20, 50], typ=typ)
				DrawEfficiencyGraphs('Eff_PH{c}_H_Vs_Threshold_snr_{s}'.format(c=ch, s=suffix), self.trans_grid.phN_snr_h_varz['PH{c}_H'.format(c=ch)], cells, tempcsnr, typ)
				# DrawEfficiencyGraphs('Eff_PH{c}_Ch_Vs_Threshold_snr_{s}'.format(c=ch, s=suffix), self.trans_grid.phN_snr_ch_varz['PH{c}_Ch'.format(c=ch)], cells, tempcsnr, typ)

	def DoPHHistograms(self, cells, cuts='', suffix='no_cuts', typ='adc'):
		def DrawHisto(name, xmin, xmax, deltax, varz, varname, cuts):
			self.trans_grid.DrawPHInArea(name, varz, cells, cuts, varname=varname, xmin=xmin, xmax=xmax, deltax=deltax)
			self.PosCanvas(name)

		tempcsnr = self.not_neg_snr_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)] if cuts == 'no_neg' else self.trans_grid.cuts_man.ConcatenateCuts(self.not_neg_snr_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)], self.not_sat_evts_region) if cuts == 'no_neg_no_sat' else '(1)'
		tempcadc = self.not_neg_adc_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)] if cuts == 'no_neg' else self.trans_grid.cuts_man.ConcatenateCuts(self.not_neg_adc_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)], self.not_sat_evts_region) if cuts == 'no_neg_no_sat' else '(1)'

		for ch in self.analysis_cummulative_ch:
			if typ == 'adc':
				DrawHisto('PH{c}_Ch_adc_{s}'.format(c=ch, s=suffix), 0, 4800, 4800. / self.trans_grid.phbins, self.trans_grid.phN_adc_ch_varz['PH{c}_Ch'.format(c=ch)], 'PH{c} cluster chs [ADC]'.format(c=ch), tempcadc)
				self.trans_grid.FitLanGaus('PH{c}_Ch_adc_{s}'.format(c=ch, s=suffix))
				if ch != self.cluster_size:
					DrawHisto('PH{c}_H_adc_{s}'.format(c=ch, s=suffix), 0, 4800, 4800. / self.trans_grid.phbins, self.trans_grid.phN_adc_h_varz['PH{c}_H'.format(c=ch)], 'PH{c} highest chs [ADC]'.format(c=ch), tempcadc)
					self.trans_grid.FitLanGaus('PH{c}_H_adc_{s}'.format(c=ch, s=suffix))
			else:
				DrawHisto('PH{c}_Ch_snr_{s}'.format(c=ch, s=suffix), 0, 480, 480. / self.trans_grid.phbins / 10.0, self.trans_grid.phN_snr_ch_varz['PH{c}_Ch'.format(c=ch)], 'PH{c} cluster chs [SNR]'.format(c=ch), tempcsnr)
				self.trans_grid.FitLanGaus('PH{c}_Ch_snr_{s}'.format(c=ch, s=suffix))
				if ch != self.cluster_size:
					DrawHisto('PH{c}_H_snr_{s}'.format(c=ch, s=suffix), 0, 480, 480. / self.trans_grid.phbins / 10.0, self.trans_grid.phN_snr_h_varz['PH{c}_H'.format(c=ch)], 'PH{c} highest chs [SNR]'.format(c=ch), tempcsnr)
					self.trans_grid.FitLanGaus('PH{c}_H_snr_{s}'.format(c=ch, s=suffix))

	def DoPH2DHistograms(self, cells, cuts='', suffix='no_cuts', typ='adc'):
		tempcsnr = self.not_neg_snr_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)] if cuts == 'no_neg' else self.trans_grid.cuts_man.ConcatenateCuts(self.not_neg_snr_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)], self.not_sat_evts_region) if cuts == 'no_neg_no_sat' else '(1)'
		tempcadc = self.not_neg_adc_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)] if cuts == 'no_neg' else self.trans_grid.cuts_man.ConcatenateCuts(self.not_neg_adc_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)], self.not_sat_evts_region) if cuts == 'no_neg_no_sat' else '(1)'
		tempcsnr = self.trans_grid.cuts_man.ConcatenateCutWithCells(tempcsnr, cells)
		tempcadc = self.trans_grid.cuts_man.ConcatenateCutWithCells(tempcadc, cells)
		for ch in self.analysis_cummulative_ch:
			nameh = 'PH{c}_H_Vs_hit_channel_adc_{s}'.format(c=ch, s=suffix) if typ == 'adc' else 'PH{c}_H_Vs_hit_channel_snr_{s}'.format(c=ch, s=suffix)
			xmin, xmax, deltax = self.trans_grid.ch_ini, self.trans_grid.ch_end, 1
			ymin, ymax, deltay = (0, 4800, 4800. / self.trans_grid.phbins) if typ == 'adc' else (0, 4800, 4800. / self.trans_grid.phbins / 10.)
			vary = self.phN_adc_h_varz['PH{c}_H'.format(c=ch)] if typ == 'adc' else self.phN_snr_h_varz['PH{c}_H'.format(c=ch)]
			tempc = tempcadc if typ == 'adc' else tempcsnr
			self.trans_grid.DrawHisto2D(nameh, xmin, xmax, deltax, 'dia pred hit ch', ymin, ymax, deltay, 'PH{c} highest chs'.format(c=ch), 'diaChannels[int(TMath::Floor(diaChXPred+0.5))]', vary, tempc)
			self.PosCanvas(nameh)
			nameh = 'PH{c}_H_Vs_event_adc_{s}'.format(c=ch, s=suffix) if typ == 'adc' else 'PH{c}_H_Vs_event_snr_{s}'.format(c=ch, s=suffix)
			xmin, xmax, deltax = self.trans_grid.trans_tree.GetMinimum('event'), self.trans_grid.trans_tree.GetMaximum('event'), 100 * self.delta_ev
			self.trans_grid.DrawHisto2D(nameh, xmin, xmax, deltax, 'event', ymin, ymax, deltay, 'PH{c} highest chs'.format(c=ch), 'event', vary, tempc)
			self.PosCanvas(nameh)

			namehp = 'PH{c}_H_mean_Vs_event_adc_{s}'.format(c=ch, s=suffix) if typ == 'adc' else 'PH{c}_H_mean_Vs_event_snr_{s}'.format(c=ch, s=suffix)
			self.trans_grid.histo[namehp] = self.trans_grid.histo[nameh].ProfileX('h_' + namehp)
			self.trans_grid.histo[namehp].SetTitle('h_' + namehp)
			self.trans_grid.canvas[namehp] = ro.TCanvas('c_' + namehp, 'c_' + namehp, 1)
			self.trans_grid.histo[namehp].Draw('e hist')
			SetDefault1DStats(self.trans_grid.histo[namehp], y1=0.15, y2=0.45)
			SetDefault1DCanvasSettings(self.trans_grid.canvas[namehp])
			ro.gPad.Update()
			self.PosCanvas(namehp)

	def DoCenterPHStudies(self, cells, cuts='', suffix='no_cuts', typ='adc', do_sat=True, do_all_plots=False):
		self.center_reg_ana.w, self.center_reg_ana.window_shift = self.w, self.window_shift
		tempcsnr = self.not_neg_snr_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)] if cuts == 'no_neg' else self.trans_grid.cuts_man.ConcatenateCuts(self.not_neg_snr_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)], self.not_sat_evts_region) if cuts == 'no_neg_no_sat' else '(1)'
		tempcadc = self.not_neg_adc_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)] if cuts == 'no_neg' else self.trans_grid.cuts_man.ConcatenateCuts(self.not_neg_adc_phN_ch['PH{c}_Ch'.format(c=self.cluster_size)], self.not_sat_evts_region) if cuts == 'no_neg_no_sat' else '(1)'
		tempc = tempcsnr if typ == 'snr' else tempcadc
		self.center_reg_ana.DoCenterRegionStudies(self.cell_dists, cells, do_all_plots, suffix, tempc, typ, do_sat)
		self.w, self.window_shift = self.center_reg_ana.w, self.center_reg_ana.window_shift

	def DoFinalAnalysis(self, typ='adc', cummulative_chs=None):
		if cummulative_chs:
			self.analysis_cummulative_ch = cummulative_chs
		self.SetFinalAnalysis()
		self.DefineSatRegion(before=0, after=1)
		self.DoDeviceMaps('all', '', 'no_cuts', typ=typ)
		self.DoStripHistograms('all', '', 'no_cuts', typ=typ)
		self.DoPH2DHistograms('all', '', 'no_cuts', typ=typ)
		self.DoEfficiencyPlots('all', '', 'no_cuts', typ=typ)
		self.DoPHHistograms('all', '', 'no_cuts', typ=typ)
		self.DoCenterPHStudies('all', '', 'no_cuts', typ=typ)
		list_cuts = ['', 'no_neg', 'no_neg_no_sat']
		cells = 'good'
		for cut in list_cuts:
			self.DoDeviceMaps(cells, cut, '{s}_{c}'.format(s=self.suffix[cells], c=cut), typ=typ)
			self.DoCellMaps(cells, cut, '{s}_{c}'.format(s=self.suffix[cells], c=cut), typ=typ)
			self.DoStripHistograms(cells, cut, '{s}_{c}'.format(s=self.suffix[cells], c=cut), typ=typ)
			self.DoPH2DHistograms(cells, cut, '{s}_{c}'.format(s=self.suffix[cells], c=cut), typ=typ)
			self.DoEfficiencyPlots(cells, cut, '{s}_{c}'.format(s=self.suffix[cells], c=cut), typ=typ)
			self.DoPHHistograms(cells, cut, '{s}_{c}'.format(s=self.suffix[cells], c=cut), typ=typ)
			self.DoCenterPHStudies(cells, cut, '{s}_{c}'.format(s=self.suffix[cells], c=cut), typ=typ, do_sat=(cut != 'no_neg_no_sat'))

	def SetFinalAnalysis(self):
		self.minz = self.trans_grid.minz
		self.maxz = self.trans_grid.maxz
		self.GetCutsFromCutManager('all')
		self.GetCutsFromCutManager('good')
		self.GetVarzFromTranspGrid()

	def GetCutsFromCutManager(self, cells):
		self.noise_cuts[cells] = self.trans_grid.cuts_man.noise_cuts[cells]
		self.noise_friend_cuts[cells] = self.trans_grid.cuts_man.noise_friend_cuts[cells]
		self.in_transp_cluster = self.trans_grid.cuts_man.ConcatenateCuts(cut2=self.trans_grid.cuts_man.in_transp_cluster, cut1=self.trans_grid.cuts_man.transp_ev)
		self.noise_nc_cuts[cells] = self.trans_grid.cuts_man.noise_nc_cuts[cells]
		self.noise_nc_friend_cuts[cells] = self.trans_grid.cuts_man.noise_nc_friend_cuts[cells]

		self.ph_adc_ch_cuts[cells] = self.trans_grid.cuts_man.ph_adc_ch_cuts[cells]
		self.ph_snr_ch_cuts[cells] = self.trans_grid.cuts_man.ph_snr_ch_cuts[cells]
		self.ph_adc_h_cuts[cells] = self.trans_grid.cuts_man.ph_adc_h_cuts[cells]
		self.ph_snr_h_cuts[cells] = self.trans_grid.cuts_man.ph_snr_h_cuts[cells]
		self.phN_adc_ch_cuts[cells] = self.trans_grid.cuts_man.phN_adc_ch_cuts[cells]
		self.phN_snr_ch_cuts[cells] = self.trans_grid.cuts_man.phN_snr_ch_cuts[cells]
		self.phN_adc_h_cuts[cells] = self.trans_grid.cuts_man.phN_adc_h_cuts[cells]
		self.phN_snr_h_cuts[cells] = self.trans_grid.cuts_man.phN_snr_h_cuts[cells]

		self.sat_adc_N_ch_cut = self.trans_grid.cuts_man.sat_adc_N_ch
		self.sat_adc_N_h_cut = self.trans_grid.cuts_man.sat_adc_N_h
		self.not_sat_adc_N_ch_cut = self.trans_grid.cuts_man.not_sat_adc_N_ch
		self.not_sat_adc_N_h_cut = self.trans_grid.cuts_man.not_sat_adc_N_h

		self.sat_evts_region = self.trans_grid.cuts_man.sat_evts_region
		self.not_sat_evts_region = self.trans_grid.cuts_man.not_sat_evts_region

		self.neg_snr_phN_ch = self.trans_grid.cuts_man.neg_snr_phN_ch
		self.neg_snr_phN_h = self.trans_grid.cuts_man.neg_snr_phN_h
		self.not_neg_snr_phN_ch = self.trans_grid.cuts_man.not_neg_snr_phN_ch
		self.not_neg_snr_phN_h = self.trans_grid.cuts_man.not_neg_snr_phN_h

		self.neg_adc_phN_ch = self.trans_grid.cuts_man.neg_adc_phN_ch
		self.neg_adc_phN_h = self.trans_grid.cuts_man.neg_adc_phN_h
		self.not_neg_adc_phN_ch = self.trans_grid.cuts_man.not_neg_adc_phN_ch
		self.not_neg_adc_phN_h = self.trans_grid.cuts_man.not_neg_adc_phN_h

	def GetVarzFromTranspGrid(self):
		self.noise_varz = self.trans_grid.noise_varz
		self.ph_adc_ch_varz = self.trans_grid.ph_adc_ch_varz
		self.ph_snr_ch_varz = self.trans_grid.ph_snr_ch_varz
		self.ph_adc_h_varz = self.trans_grid.ph_adc_h_varz
		self.ph_snr_h_varz = self.trans_grid.ph_snr_h_varz
		self.phN_adc_ch_varz = self.trans_grid.phN_adc_ch_varz
		self.phN_snr_ch_varz = self.trans_grid.phN_snr_ch_varz
		self.phN_adc_h_varz = self.trans_grid.phN_adc_h_varz
		self.phN_snr_h_varz = self.trans_grid.phN_snr_h_varz

	def DefineSatRegion(self, before=0, after=1):
		if self.trans_grid.trans_tree.GetFriend('satRegions'):
			self.trans_grid.UnfriendTree(self.trans_grid.trans_tree.GetFriend('satRegions'))
		self.trans_grid.AddFriendWithSaturationRegions(skipAfter=after, skipBefore=before)

if __name__ == '__main__':
	c = FinalAnalysis(None, 0, 0)
