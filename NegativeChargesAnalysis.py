#!/usr/bin/env python
import ROOT as ro
import numpy as np
import ipdb, os
from ConfigParser import ConfigParser
from TransparentGrid import TransparentGrid
from optparse import OptionParser
from Utils import *

color_index = 10000

class NegativeChargesAnalysis:
	def __init__(self, trans_grid, numstrips, clustersize, noise_ana):
		self.window_shift = 3
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

		self.suffix = {'all': 'all', 'good': 'selection', 'bad': 'not_selection'}

		self.neg_ch_cut = self.trans_grid.cuts_man.GetNegPHChCut
		self.negN_chs_cut = self.trans_grid.cuts_man.GetNegPHNChsCut
		self.not_neg_ch_cut = self.trans_grid.cuts_man.GetNotNegPHChCut
		self.not_negN_chs_cut = self.trans_grid.cuts_man.GetNotNegPHNChsCut

		self.ph_ch_var = self.trans_grid.GetPHChVar
		self.phN_chs_var = self.trans_grid.GetPHNChsVar

		self.strips_for_analysis = np.arange(self.cluster_size, dtype='i4')


	def PosCanvas(self, canvas_name):
		self.w = PositionCanvas(self.trans_grid, canvas_name, self.w, self.window_shift)

	def OverlayNoiseDistribution(self, histo, cells='all', isFriend=False):
		if self.noise_ana:
			self.noise_ana.OverlayNoiseDistribution(histo, cells, isFriend)

	# def Do1DChSignalHistos(self, cells='all', doLog=False, typ='adc', isFriend=False):
	# 	def DrawHisto(name, histo_limits, plot_lims, deltax, varz, varname, cuts):
	# 		self.trans_grid.DrawHisto1D(name, histo_limits['min'], histo_limits['max'], deltax, varz, varname, cuts)
	# 		self.trans_grid.histo[name].GetXaxis().SetRangeUser(plot_lims['min'], plot_lims['max'])
	# 		SetX1X2NDC(self.trans_grid.histo[name], 0.15, 0.45, 'stats')
	# 		self.OverlayNoiseDistribution(self.trans_grid.histo[name], cells, isFriend)
	# 		if doLog: self.trans_grid.canvas[name].SetLogy()
	# 		legend = self.trans_grid.canvas[name].BuildLegend()
	# 		ro.gPad.Update()
	# 		SetLegendX1X2Y1Y2(legend, 0.15, 0.45, 0.5, 0.6)
	# 		self.PosCanvas(name)
    #
	# 	suffix = self.suffix[cells] if cells in self.suffix.keys() else ''
	# 	suffix = suffix + '_logScale' if doLog else suffix
    #
	# 	tempcutsadc = self.neg_adc_phN_h['PH{c}_H'.format(c=self.cluster_size)]
	# 	tempcutssnr = self.neg_snr_phN_h['PH{c}_H'.format(c=self.cluster_size)]
	# 	for ch in xrange(self.cluster_size):
	# 		if typ == 'snr':
	# 			if 'PH_Ch' + str(ch) in self.ph_snr_ch_varz.keys():
	# 				tempcuts = tempcutssnr
	# 				minz, maxz = self.min_snr_neg, self.max_snr_neg
	# 				hist_limits = Get1DLimits(minz, maxz, self.delta_snr, 2)
	# 				plot_limits = Get1DLimits(minz, maxz, self.delta_snr, 1)
	# 				DrawHisto('PH_Ch{c}_snr_neg_evts_{s}'.format(c=ch, s=suffix), hist_limits, plot_limits, self.delta_snr, self.ph_snr_ch_varz['PH_Ch{c}'.format(c=ch)], 'PH cluster ch{c} neg events [SNR]'.format(c=ch), tempcuts)
    #
	# 			if 'PH_H' + str(ch+1) in self.ph_snr_h_varz.keys():
	# 				tempcuts = tempcutssnr
	# 				minz, maxz = self.min_snr_neg, self.max_snr_neg
	# 				hist_limits = Get1DLimits(minz, maxz, self.delta_snr, 2)
	# 				plot_limits = Get1DLimits(minz, maxz, self.delta_snr, 1)
	# 				DrawHisto('PH_H{c}_snr_neg_evts_{s}'.format(c=ch+1, s=suffix), hist_limits, plot_limits, self.delta_snr, self.ph_snr_h_varz['PH_H{c}'.format(c=ch+1)], 'PH highest {c}{sf} ch neg events [SNR]'.format(c=ch+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), tempcuts)
    #
	# 		else:
	# 			if 'PH_Ch' + str(ch) in self.ph_adc_ch_varz.keys():
	# 				tempcuts = tempcutsadc
	# 				minz, maxz = self.min_adc_neg, self.max_adc_neg
	# 				hist_limits = Get1DLimits(minz, maxz, self.delta_adc, 2)
	# 				plot_limits = Get1DLimits(minz, maxz, self.delta_adc, 1)
	# 				DrawHisto('PH_Ch{c}_adc_neg_evts_{s}'.format(c=ch, s=suffix), hist_limits, plot_limits, self.delta_adc, self.ph_adc_ch_varz['PH_Ch{c}'.format(c=ch)], 'PH cluster ch{c} neg events [ADC]'.format(c=ch), tempcuts)
    #
	# 			if 'PH_H' + str(ch+1) in self.ph_adc_h_varz.keys():
	# 				tempcuts = tempcutsadc
	# 				minz, maxz = self.min_adc_neg, self.max_adc_neg
	# 				hist_limits = Get1DLimits(minz, maxz, self.delta_adc, 2)
	# 				plot_limits = Get1DLimits(minz, maxz, self.delta_adc, 1)
	# 				DrawHisto('PH_H{c}_adc_neg_evts_{s}'.format(c=ch+1, s=suffix), hist_limits, plot_limits, self.delta_adc, self.ph_adc_h_varz['PH_H{c}'.format(c=ch+1)], 'PH highest {c}{sf} ch neg events [ADC]'.format(c=ch+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), tempcuts)

	def DoProfileMaps(self, typ='adc', isFriend=False):
		xmin, xmax, deltax, xname = self.trans_grid.ch_ini - 1.5, self.trans_grid.ch_end + 1.5, 1.0/self.trans_grid.bins_per_ch_x, 'pred dia hit ch',
		ymin = min(self.trans_grid.row_cell_info_diamond['0_even'], self.trans_grid.row_cell_info_diamond['0_odd'])
		yup_max = max(self.trans_grid.row_cell_info_diamond['up_even'], self.trans_grid.row_cell_info_diamond['up_odd'])
		ymin = ymin - RoundInt(float(ymin) / self.trans_grid.row_cell_info_diamond['height'], 'f8') * self.trans_grid.row_cell_info_diamond['height']
		ymax, deltay, yname = ymin + 256 * 50., float(self.trans_grid.row_cell_info_diamond['height'])/self.trans_grid.bins_per_ch_y, 'sil pred dia hit in Y [#mum]'

		def DrawProfile2D(name, varz, varzname, cut, xdelt=deltax, namex=xname, varx='diaChXPred', getOccupancy=False):
			self.trans_grid.DrawProfile2D(name, xmin, xmax, xdelt, namex, ymin, ymax, deltay, yname, varx, 'diaChYPred', varz, varzname, cut)
			self.trans_grid.profile[name].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile[name].GetYaxis().SetRangeUser(ymin - int(self.trans_grid.row_cell_info_diamond['height'] / self.trans_grid.bins_per_ch_y), yup_max + int(self.trans_grid.row_cell_info_diamond['height'] / self.trans_grid.bins_per_ch_y))
			self.trans_grid.DrawTCutGs(name, 'diamond')
			self.trans_grid.DrawGoodAreasDiamondCenters(name)
			self.PosCanvas(name)
			if getOccupancy:
				self.trans_grid.GetOccupancyFromProfile(name)
				self.trans_grid.histo['hit_map_' + name].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
				self.trans_grid.histo['hit_map_' + name].GetYaxis().SetRangeUser(ymin - int(self.trans_grid.row_cell_info_diamond['height'] / self.trans_grid.bins_per_ch_y), yup_max + int(self.trans_grid.row_cell_info_diamond['height'] / self.trans_grid.bins_per_ch_y))
				self.trans_grid.DrawTCutGs('hit_map_' + name, 'diamond')
				self.trans_grid.DrawGoodAreasDiamondCenters('hit_map_' + name)
				self.PosCanvas('hit_map_' + name)

		for ch in self.strips_for_analysis:
			tempc = self.neg_ch_cut(ch, 'Ch', typ == 'snr', isFriend)
			hname = 'PH_Ch{c}_neg_map_pred'.format(c=ch) if not isFriend else 'PH_Ch{c}_buffer_{v}_neg_map_pred'.format(c=ch, v=self.trans_grid.noise_friend_buffer)
			DrawProfile2D('{n}_hit_snr'.format(n=hname), self.ph_ch_var(ch, 'Ch', typ == 'snr', isFriend), 'PH cluster ch{c} [{t}]'.format(c=ch, t=typ.upper()), tempc, getOccupancy=True)
			DrawProfile2D('{n}_ch_snr'.format(n=hname), self.ph_ch_var(ch, 'Ch', typ == 'snr', isFriend), 'PH cluster ch{c} [{t}]'.format(c=ch, t=typ.upper()), tempc, 1, 'cluster ch{c}'.format(c=ch), 'clusterChannel{c}'.format(c=ch))
			# tempc = self.neg_ch_cut(ch + 1, 'H', typ == 'snr', isFriend)
			# hname = 'PH_H{c}_neg_map_pred'.format(c=ch + 1) if not isFriend else 'PH_H{c}_buffer_{v}_neg_map_pred'.format(c=ch + 1, v=self.trans_grid.noise_friend_buffer)
			# DrawProfile2D('{n}_hit_snr'.format(n=hname), self.ph_ch_var(ch + 1, 'H', typ == 'snr', isFriend), 'PH highest {c}{sf} ch [{t}]'.format(c=ch+1, t=typ.upper(), sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), tempc, getOccupancy=True)
			# DrawProfile2D('{n}_hit_snr'.format(n=hname), self.ph_ch_var(ch + 1, 'H', typ == 'snr', isFriend), 'PH highest {c}{sf} ch [{t}]'.format(c=ch+1, t=typ.upper(), sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), tempc, 1, 'highest {c}{sf} ch'.format(c=ch+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), 'clusterChannelHighest{c}'.format(c=ch+1))

	def Do1DPHHistos(self, cells='all', typ='adc', isFriend=False):
		suffix = self.suffix[cells]
		def DrawPHHisto(name, varz, varzname, cuts):
			self.trans_grid.DrawPHInArea(name, varz, cells, cuts, varname=varzname, typ=typ)
			self.PosCanvas(name)

		tempc = self.negN_chs_cut(self.cluster_size, 'Ch', typ == 'snr', isFriend)
		for ch in (self.strips_for_analysis + 1):
			hname = 'PH{c}_Ch_neg_events_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower()) if not isFriend else 'PH{c}_Ch_buffer_{v}_neg_events_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
			DrawPHHisto(hname, self.phN_chs_var(ch, 'Ch', typ == 'snr', isFriend), 'PH{c} cluster chs [{t}]'.format(c=ch, t=typ.upper()), tempc)
			if ch != self.cluster_size:
				hname = 'PH{c}_H_neg_events_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower()) if not isFriend else 'PH{c}_H_buffer_{v}_neg_events_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
				DrawPHHisto(hname, self.phN_chs_var(ch, 'H', typ == 'snr', isFriend), 'PH{c} highest chs [{t}]'.format(c=ch, t=typ.upper()), tempc)

	def DoCellMaps(self, cells='all', typ='adc', isFriend=False):
		def PlotCellsProfiles(name, varz, varname, cut):
			self.trans_grid.DrawProfile2DDiamondCellOverlay(name, varz, cells, cut, varname=varname, typ=typ)
			self.PosCanvas(name)

		suffix = self.suffix[cells]
		tempc = self.negN_chs_cut(self.cluster_size, 'Ch', typ == 'snr', isFriend)
		for ch in (self.strips_for_analysis + 1):
			hname = 'PH{c}_H_cell_map_neg_events_snr_{s}'.format(c=ch, s=suffix) if not isFriend else 'PH{c}_H_buffer_{v}_cell_map_neg_events_snr_{s}'.format(c=ch, s=suffix, v=self.trans_grid.noise_friend_buffer)
			PlotCellsProfiles(hname, self.phN_chs_var(ch, 'H', typ == 'snr', isFriend), 'PH{c} highest chs [{t}]'.format(c=ch, t=typ.upper()), tempc)
			# hname = 'PH{c}_Ch_cell_map_neg_events_snr_{s}'.format(c=ch, s=suffix) if not isFriend else 'PH{c}_Ch_buffer_{v}_cell_map_neg_events_snr_{s}'.format(c=ch, s=suffix, v=self.trans_grid.noise_friend_buffer)
			# PlotCellsProfiles(hname, self.phN_chs_var(ch, 'Ch', typ == 'snr', isFriend), 'PH{c} cluster chs [{t}]'.format(c=ch, t=typ.upper()), tempc)

	def PlotStripHistograms(self, cells='all', typ='adc', isFriend=False):
		minx, maxx, deltax, xname, xvar = -0.5, 0.5, self.trans_grid.cell_resolution / float(self.trans_grid.row_cell_info_diamond['height']), 'dia pred. strip hit pos', 'diaChXPred-TMath::Floor(diaChXPred+0.5)'
		def Draw2DHistogram(name, zmin, zmax, yname, yvar, cuts, typ='adc'):
			histo_limits = Get1DLimits(zmin, zmax, 4 * self.delta_adc) if typ == 'adc' else Get1DLimits(zmin, zmax, 4 * self.delta_snr)
			deltay = 4 * self.delta_adc if typ == 'adc' else 4 * self.delta_snr
			self.trans_grid.DrawHisto2D(name, minx, maxx, deltax, xname, min(0, histo_limits['min']), histo_limits['max'], deltay, yname, xvar, yvar, cuts)
			self.PosCanvas(name)

		def Draw1DHistogram(name, cuts):
			self.trans_grid.DrawHisto1D(name, minx, maxx, deltax, xvar, xname, cuts)
			maxbin = self.trans_grid.histo[name].GetMaximumBin()
			self.trans_grid.histo[name].GetYaxis().SetRangeUser(0, self.trans_grid.histo[name].GetBinContent(maxbin) + self.trans_grid.histo[name].GetBinError(maxbin))
			self.PosCanvas(name)

		suffix = self.suffix[cells] if cells in self.suffix.keys() else ''
		tempc = self.negN_chs_cut(self.cluster_size, 'Ch', typ == 'snr', isFriend)
		for ch in (self.strips_for_analysis + 1):
			minz, maxz = self.trans_grid.minz['PH_Ch{c}_{t}'.format(c=ch-1, t=typ.lower())], self.trans_grid.maxz['PH_Ch{c}_{t}'.format(c=ch-1, t=typ.lower())]
			hname = 'PH_Ch{c}_Vs_strip_location_neg_events_{t}_{s}'.format(c=ch-1, s=suffix, t=typ.lower()) if not isFriend else 'PH_Ch{c}_buffer_{v}_Vs_strip_location_neg_events_{t}_{s}'.format(c=ch-1, s=suffix, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
			Draw2DHistogram(hname, minz, maxz, 'PH cluster ch{c} [{t}]'.format(c=ch-1, t=typ.upper()), self.ph_ch_var(ch - 1, 'Ch', typ == 'snr', isFriend), tempc, typ)
			minz, maxz = self.trans_grid.minz['PH_H{c}_{t}'.format(c=ch, t=typ.lower())], self.trans_grid.maxz['PH_H{c}_{t}'.format(c=ch, t=typ.lower())]
			hname = 'PH_H{c}_Vs_strip_location_neg_events_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower()) if not isFriend else 'PH_H{c}_Vs_strip_location_neg_events_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
			Draw2DHistogram(hname, minz, maxz, 'PH highest {c}{sf} ch [{t}]'.format(c=ch, t=typ.upper(), sf='st' if ch - 1 == 0 else 'nd' if ch - 1 == 1 else 'rd' if ch - 1 == 2 else 'th'), self.ph_ch_var(ch, 'H', typ == 'snr', isFriend), tempc, typ)
			minz, maxz = self.trans_grid.minz['PH{c}_H_{t}'.format(c=ch, t=typ.lower())], self.trans_grid.maxz['PH{c}_H_{t}'.format(c=ch, t=typ.lower())]
			hname = 'PH{c}_H_Vs_strip_location_neg_events_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower()) if not isFriend else 'PH{c}_H_buffer_{v}_Vs_strip_location_neg_events_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower(), v=self.trans_grid.noise_friend_buffer)
			Draw2DHistogram(hname, minz, maxz, 'PH{c} highest chs [{t}]'.format(c=ch, t=typ.upper()), self.phN_chs_var(ch, 'H', typ == 'snr', isFriend), tempc, typ)
		Draw1DHistogram('strip_location_neg_events_{t}_{s}'.format(s=suffix, t=typ.lower()), tempc)

	def DoNegativeAnalysis(self, cells='all', typ='adc', isFriend=False):
		self.SetStripsForAnalysis()
		# self.Do1DChSignalHistos(cells, False, typ, isFriend)
		self.DoProfileMaps(typ=typ, isFriend=isFriend)
		self.DoCellMaps(cells, typ, isFriend)
		self.PlotStripHistograms(cells, typ, isFriend)
		self.Do1DPHHistos(cells, typ, isFriend)

	def SetStripsForAnalysis(self, arr=None):
		if arr:
			self.strips_for_analysis = arr
		else:
			self.strips_for_analysis = np.arange(self.num_strips, dtype='i4')
			self.strips_for_analysis = self.strips_for_analysis if self.cluster_size > 3 else np.unique(np.append(self.strips_for_analysis, self.cluster_size - 1))


if __name__ == '__main__':
	c = NegativeChargesAnalysis(None, 0, 0)
