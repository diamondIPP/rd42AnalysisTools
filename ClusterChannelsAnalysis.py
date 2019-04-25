#!/usr/bin/env python
import ROOT as ro
import numpy as np
import ipdb, os
from ConfigParser import ConfigParser
from TransparentGrid import TransparentGrid
from optparse import OptionParser
from Utils import *

color_index = 10000

class ClusterChannelsAnalysis:
	def __init__(self, trans_grid, numstrips, clustersize, noise_ana):
		self.window_shift = 3
		# self.min_snr_neg, self.max_snr_neg, self.delta_snr = -64.25, 0.25, 0.125
		self.min_snr, self.max_snr = -650, 650
		self.min_adc, self.max_adc = -6500, 6500
		self.delta_adc, self.delta_snr = 20, 2
		self.min_adc_noise, self.max_adc_noise, self.delta_adc_noise = -322.5, 322.5, 0.5
		# self.min_adc_noise, self.max_adc_noise, self.delta_adc_noise = -32.25, 32.25, 0.5
		# self.min_snr_noise, self.max_snr_noise, self.delta_snr_noise = -3.225, 3.225, 0.05
		self.min_snr_noise, self.max_snr_noise, self.delta_snr_noise = -32.25, 32.25, 0.05
		self.neg_cut_lines = {}
		self.trash = []
		self.w = 0
		self.trans_grid = trans_grid
		self.minz = self.trans_grid.minz
		self.maxz = self.trans_grid.maxz
		self.num_strips = numstrips
		self.cluster_size = clustersize
		self.noise_ana = noise_ana

		self.suffix = {'all': 'all', 'good': 'selection', 'bad': 'not_selection'}

		self.noise_cuts = self.trans_grid.cuts_man.noise_cuts

		self.ph_cuts = self.trans_grid.cuts_man.GetPHCuts

		self.ph_ch_var = self.trans_grid.GetPHChVar

		self.phN_chs_var = self.trans_grid.GetPHNChsVar

		self.noise_varz = self.trans_grid.noise_varz

		self.strips_for_analysis = np.arange(self.cluster_size, 'i4')

	def PosCanvas(self, canvas_name):
		self.w = PositionCanvas(self.trans_grid, canvas_name, self.w, self.window_shift)

	def OverlayNoiseDistribution(self, histo, cells='all', isFriend=False):
		if self.noise_ana:
			self.noise_ana.OverlayNoiseDistribution(histo, cells, isFriend)

	def Do1DHistograms(self, cells='all', doLog=False, typ='adc', isFriend=False):
		def DrawHisto(name, histo_limits, plot_lims, deltax, varz, varname, cuts):
			tempc = self.trans_grid.cuts_man.ConcatenateCutWithCells(cut=cuts, cells=cells)
			self.trans_grid.DrawHisto1D(name, histo_limits['min'], histo_limits['max'], deltax, varz, varname, tempc)
			self.trans_grid.histo[name].GetXaxis().SetRangeUser(plot_lims['min'], plot_lims['max'])
			SetX1X2NDC(self.trans_grid.histo[name], 0.15, 0.45, 'stats')
			self.OverlayNoiseDistribution(self.trans_grid.histo[name], cells, isFriend)
			if doLog: self.trans_grid.canvas[name].SetLogy()
			legend = self.trans_grid.canvas[name].BuildLegend()
			ro.gPad.Update()
			SetLegendX1X2Y1Y2(legend, 0.15, 0.45, 0.5, 0.6)
			self.PosCanvas(name)

		suffix = self.suffix[cells] if cells in self.suffix.keys() else ''
		suffix = suffix + '_logScale' if doLog else suffix

		for ch in self.strips_for_analysis:
			tempcuts = self.ph_cuts('PH_Ch{i}'.format(i=ch), isFriend)
			minz, maxz = self.minz['PH_Ch{i}_{t}'.format(i=ch, t=typ.lower())], self.maxz['PH_Ch{i}_{t}'.format(i=ch, t=typ.lower())]
			delta = self.delta_adc if typ == 'adc' else self.delta_snr
			hist_limits = GetSymmetric1DLimits(minz, maxz, delta, 2)
			plot_limits = GetSymmetric1DLimits(minz, maxz, delta, 1, False)
			hname = 'PH_Ch{i}_{t}_{s}'.format(i=ch, s=suffix, t=typ.lower()) if not isFriend else 'PH_Ch{i}_buffer_{v}_{t}_{s}'.format(v=self.trans_grid.noise_friend_buffer, i=ch, s=suffix, t=typ.lower())
			DrawHisto(hname, hist_limits, plot_limits, delta, self.ph_ch_var(ch, 'Ch', typ=='snr', isFriend), 'PH cluster ch{i} [{t}]'.format(i=ch, t=typ.upper()), tempcuts)

			tempcuts = self.ph_cuts('PH_H{i}'.format(i=ch+1), isFriend)
			minz, maxz = self.minz['PH_H{i}_{t}'.format(i=ch+1, t=typ.lower())], self.maxz['PH_H{i}_{t}'.format(i=ch+1, t=typ.lower())]
			hist_limits = GetSymmetric1DLimits(minz, maxz, delta, 2)
			plot_limits = GetSymmetric1DLimits(minz, maxz, delta, 1, False)
			hname = 'PH_H{i}_{t}_{s}'.format(i=ch + 1, s=suffix, t=typ.lower()) if not isFriend else 'PH_H{i}_buffer_{v}_{t}_{s}'.format(v=self.trans_grid.noise_friend_buffer, i=ch+1, s=suffix, t=typ.lower())
			DrawHisto(hname, hist_limits, plot_limits, delta, self.ph_ch_var(ch+1, 'H', typ=='snr', isFriend), 'PH highest {i}{sf} ch [{t}]'.format(i=ch + 1, t=typ.upper(), sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), tempcuts)

	def Do2DProfileMaps(self, cells='all', isFriend=False):
		def DrawProfile(name, varx, xname, varz, varname, cuts, zmax, zmin):
			self.trans_grid.DrawProfile2DDiamondChannel(name, varx, xname, varz, varname, cells, cuts)
			self.trans_grid.profile[name].SetMaximum(zmax)
			self.trans_grid.profile[name].SetMinimum(zmin)
			self.trans_grid.DrawTCutGs(name, 'diamond')
			self.PosCanvas(name)

		suffix = self.suffix[cells] if cells in self.suffix.keys() else ''
		for ch in self.strips_for_analysis:
			tempcuts = self.ph_cuts('PH_Ch' + str(ch), isFriend)
			minz, maxz = self.minz['PH_Ch{i}_adc'.format(i=ch)], self.maxz['PH_Ch{i}_adc'.format(i=ch)]
			hname = 'PH_Ch{c}_map_{s}'.format(c=ch, s=suffix) if not isFriend else 'PH_Ch{c}_buffer_{v}_map_{s}'.format(c=ch, s=suffix, v=self.trans_grid.noise_friend_buffer)
			DrawProfile(hname, 'clusterChannel{c}'.format(c=ch), 'cluster ch{c}'.format(c=ch), self.ph_ch_var(ch, 'Ch', False, isFriend), 'PH [ADC]', tempcuts, maxz, minz)

			tempcuts = self.ph_cuts('PH_H' + str(ch + 1), isFriend)
			minz, maxz = self.minz['PH_H{i}_adc'.format(i=ch+1)], self.maxz['PH_H{i}_adc'.format(i=ch+1)]
			hname = 'PH_H{c}_map_{s}'.format(c=ch+1, s=suffix) if not isFriend else 'PH_H{c}_buffer_{v}_map_{s}'.format(c=ch+1, s=suffix, v=self.trans_grid.noise_friend_buffer)
			DrawProfile(hname, 'clusterChannelHighest{c}'.format(c=ch+1), 'highest {c}{sf} ch'.format(c=ch+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), self.ph_ch_var(ch + 1, 'H', False, isFriend), 'PH [ADC]', tempcuts, maxz, minz)

	def DoStrips2DHistograms(self, cells='all', typ='adc', isFriend=False):
		minx, maxx, deltax, xname, xvar = -0.5, 0.5, self.trans_grid.cell_resolution / float(self.trans_grid.row_info_diamond['pitch']), 'dia pred. strip hit pos', 'diaChXPred-TMath::Floor(diaChXPred+0.5)'

		def DrawHistogram(name, zmin, zmax, yname, yvar, cuts, typ='adc'):
			deltay = 4 * self.delta_adc if typ == 'adc' else 4 * self.delta_snr
			histo_limits = Get1DLimits(zmin, zmax, deltay)
			tempc = self.trans_grid.cuts_man.ConcatenateCutWithCells(cut=cuts, cells=cells)
			self.trans_grid.DrawHisto2D(name, minx, maxx, deltax, xname, min(0, histo_limits['min']), histo_limits['max'], deltay, yname, xvar, yvar, tempc)
			self.PosCanvas(name)

		suffix = self.suffix[cells] if cells in self.suffix.keys() else ''
		for ch in self.strips_for_analysis:
			tempcuts = self.ph_cuts('PH_Ch{i}'.format(i=ch), isFriend)
			minz, maxz = self.minz['PH_Ch{c}_{t}'.format(c=ch, t=typ.lower())], self.maxz['PH_Ch{c}_{t}'.format(c=ch, t=typ.lower())]
			hname = 'PH_Ch{c}_Vs_strip_location_{t}_{s}'.format(c=ch, s=suffix, t=typ.lower()) if not isFriend else 'PH_Ch{c}_buffer_{v}_Vs_strip_location_{t}_{s}'.format(v=self.trans_grid.noise_friend_buffer, c=ch, s=suffix, t=typ.lower())
			DrawHistogram(hname, minz, maxz, 'PH cluster ch{c} [{t}]'.format(c=ch, t=typ.upper()), self.ph_ch_var(ch, 'Ch', typ=='snr', isFriend), tempcuts)

			tempcuts = self.ph_cuts('PH_H{i}'.format(i=ch+1), isFriend)
			minz, maxz = self.minz['PH_H{c}_{t}'.format(c=ch+1, t=typ.lower())], self.maxz['PH_H{c}_{t}'.format(c=ch+1, t=typ.lower())]
			hname = 'PH_H{c}_Vs_strip_location_{t}_{s}'.format(c=ch+1, s=suffix, t=typ.lower()) if not isFriend else 'PH_H{c}_buffer_{v}_Vs_strip_location_{t}_{s}'.format(v=self.trans_grid.noise_friend_buffer, c=ch+1, s=suffix, t=typ.lower())
			DrawHistogram(hname, minz, maxz, 'PH highest {c}{sf} ch [{t}]'.format(t=typ.upper(), c=ch+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), self.ph_ch_var(ch+1, 'H', typ=='snr', isFriend), tempcuts)

			tempcuts = self.ph_cuts('PH{i}_Ch'.format(i=ch+1), isFriend)
			minz, maxz = self.minz['PH{c}_Ch_{t}'.format(c=ch+1, t=typ.lower())], self.maxz['PH{c}_Ch_{t}'.format(c=ch+1, t=typ.lower())]
			hname = 'PH{c}_Ch_Vs_strip_location_{t}_{s}'.format(c=ch+1, s=suffix, t=typ.lower()) if not isFriend else 'PH{c}_Ch_buffer_{v}_Vs_strip_location_{t}_{s}'.format(v=self.trans_grid.noise_friend_buffer, c=ch+1, s=suffix, t=typ.lower())
			DrawHistogram(hname, minz, maxz, 'PH{c} cluster chs [{t}]'.format(c=ch+1, t=typ.upper()), self.phN_chs_var(ch + 1, 'Ch', typ=='snr', isFriend), tempcuts)

			if ch != self.cluster_size - 1:
				# these will be the same as PH{c1}_Ch when c1 == self.cluster_size - 1, because it is all the channels in the cluster
				tempcuts = self.ph_cuts('PH{i}_H'.format(i=ch+1), isFriend)
				minz, maxz = self.minz['PH{c}_H_{t}'.format(c=ch+1, t=typ.lower())], self.maxz['PH{c}_H_{t}'.format(c=ch+1, t=typ.lower())]
				hname = 'PH{c}_H_Vs_strip_location_{t}_{s}'.format(c=ch+1, s=suffix, t=typ.lower()) if not isFriend else 'PH{c}_H_buffer_{v}_Vs_strip_location_{t}_{s}'.format(v=self.trans_grid.noise_friend_buffer, c=ch+1, s=suffix, t=typ.lower())
				DrawHistogram(hname, minz, maxz, 'PH{c} highest chs [{t}]'.format(c=ch + 1, t=typ.upper()), self.phN_chs_var(ch + 1, 'H', typ=='snr', isFriend), tempcuts)

	def DoPHStripCorrelations(self, cells='all', typ='adc', isFriend=False):
		xmin_adc, xmax_adc = -550, 2650
		xmin_snr, xmax_snr = -55, 265

		def DrawHistogram(name, xname, ymin, ymax, yname, varx, vary, cuts, typ='adc'):
			if typ == 'adc':
				self.trans_grid.DrawHisto2D(name, xmin_adc, xmax_adc, self.delta_adc * 4, xname, ymin, ymax, self.delta_adc * 4, yname, varx, vary, cuts)
			else:
				self.trans_grid.DrawHisto2D(name, xmin_snr, xmax_snr, self.delta_snr * 4, xname, ymin, ymax, self.delta_snr * 4, yname, varx, vary, cuts)
			self.PosCanvas(name)

		for ch in self.strips_for_analysis:
			for ch2 in self.strips_for_analysis:
				if ch != self.cluster_size - 1:
					# these will be the same as PH{c1}_Ch when c1 == self.cluster_size - 1, because it is all the channels in the cluster
					ymin, ymax = (0, 4200) if typ == 'adc' else (0, 420)
					tempcuts = self.trans_grid.cuts_man.ConcatenateCutWithCells(cut=self.trans_grid.cuts_man.AndCuts([self.ph_cuts('PH{c1}_H'.format(c1=ch+1), isFriend), self.ph_cuts('PH_H{c2}'.format(c2=ch2+1),isFriend)]), cells=cells)
					DrawHistogram('PH{c1}_H_Vs_PH_H{c2}_{t}'.format(c1=ch+1, c2=ch2+1, t=typ.lower()), 'PH highest {c}{sf} ch [{t}]'.format(c=ch2+1, t=typ.upper(), sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), ymin, ymax, 'PH{c1} highest chs [{t}]'.format(c1=ch+1, t=typ.upper()), self.ph_ch_var(ch2+1, 'H', typ=='snr', isFriend), self.phN_chs_var(ch+1, 'H', typ=='snr', isFriend), tempcuts, typ)

					tempcuts = self.trans_grid.cuts_man.ConcatenateCutWithCells(cut=self.trans_grid.cuts_man.AndCuts([self.ph_cuts('PH{c1}_H'.format(c1=ch+1), isFriend), self.ph_cuts('PH_Ch{c2}'.format(c2=ch2), isFriend)]), cells=cells)
					DrawHistogram('PH{c1}_H_Vs_PH_Ch{c2}_{t}'.format(c1=ch+1, c2=ch2, t=typ.lower()), 'PH cluster ch{c2} [{t}]'.format(c2=ch2, t=typ.upper()), ymin, ymax, 'PH{c1} highest chs [{t}]'.format(c1=ch+1, t=typ.upper()), self.ph_ch_var(ch2, 'Ch', typ=='snr', isFriend), self.phN_chs_var(ch+1, 'H', typ=='snr', isFriend), tempcuts, typ)

					ymin2, ymax2 = (-550, 2650) if typ == 'adc' else (-55, 265)
					tempcuts = self.trans_grid.cuts_man.ConcatenateCutWithCells(cut=self.trans_grid.cuts_man.AndCuts([self.ph_cuts('PH{c1}_Ch'.format(c1=ch+1), isFriend), self.ph_cuts('PH_Ch{c2}'.format(c2=ch2), isFriend)]), cells=cells)
					DrawHistogram('PH{c1}_Ch_Vs_PH_Ch{c2}_{t}'.format(c1=ch+1, c2=ch2, t=typ.lower()), 'PH cluster ch{c2} [{t}]'.format(c2=ch2, t=typ.upper()), ymin, ymax, 'PH{c1} cluster chs [{t}]'.format(c1=ch+1, t=typ.upper()), self.ph_ch_var(ch2, 'Ch', typ=='snr', isFriend), self.phN_chs_var(ch+1, 'Ch', typ=='snr', isFriend), tempcuts, typ)

					tempcuts = self.trans_grid.cuts_man.ConcatenateCutWithCells(cut=self.trans_grid.cuts_man.AndCuts([self.ph_cuts('PH{c1}_Ch'.format(c1=ch+1), isFriend), self.ph_cuts('PH_H{c2}'.format(c2=ch2+1), isFriend)]), cells=cells)
					DrawHistogram('PH{c1}_Ch_Vs_PH_H{c2}_{t}'.format(c1=ch+1, c2=ch2+1, t=typ.lower()), 'PH highest {c2}{sf} ch [{t}]'.format(c2=ch2+1, t=typ.upper(), sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), ymin, ymax, 'PH{c1} cluster chs [{t}]'.format(c1=ch+1, t=typ.upper()), self.ph_ch_var(ch2+1, 'H', typ=='snr', isFriend), self.phN_chs_var(ch+1, 'Ch', typ=='snr', isFriend), tempcuts, typ)

					tempcuts = self.trans_grid.cuts_man.ConcatenateCutWithCells(cut=self.trans_grid.cuts_man.AndCuts([self.ph_cuts('PH_H{c1}'.format(c1=ch+1), isFriend), self.ph_cuts('PH_Ch{c2}'.format(c2=ch2), isFriend)]), cells=cells)
					DrawHistogram('PH_H{c1}_Vs_PH_Ch{c2}_{t}'.format(c1=ch+1, c2=ch2, t=typ.lower()), 'PH cluster ch{c2} [{t}]'.format(c2=ch2, t=typ.upper()), ymin2, ymax2, 'PH highest {c1}{sf} ch [{t}]'.format(c1=ch+1, t=typ.upper(), sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), self.ph_ch_var(ch2, 'Ch', typ=='snr', isFriend), self.ph_ch_var(ch+1, 'H', typ=='snr', isFriend), tempcuts, typ)

	def DoClusterStudies(self, cells='all', typ='adc', isFriend=False):
		self.SetStripsForAnalysis()
		self.Do2DProfileMaps(cells)
		self.DoStrips2DHistograms(cells, typ=typ, isFriend=isFriend)
		self.DoPHStripCorrelations(cells, typ=typ, isFriend=isFriend)
		self.Do1DHistograms(cells, False, typ=typ, isFriend=isFriend)
		# self.Do1DHistograms(cells, True, typ='adc')

	def SetStripsForAnalysis(self, arr=None):
		if arr:
			self.strips_for_analysis = arr
		else:
			self.strips_for_analysis = np.arange(self.num_strips, 'i4')
			self.strips_for_analysis = self.strips_for_analysis if self.cluster_size > 3 else np.unique(np.append(self.strips_for_analysis, self.cluster_size - 1))

if __name__ == '__main__':
	c = ClusterChannelsAnalysis(None, 0, 0)
