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
		self.trans_grid = trans_grid
		self.delta_adc, self.delta_snr = self.trans_grid.delta_adc_cluster_ch, self.trans_grid.delta_adc_cluster_ch / 10.
		self.min_adc_noise, self.max_adc_noise, self.delta_adc_noise = -322.5, 322.5, 0.5
		# self.min_adc_noise, self.max_adc_noise, self.delta_adc_noise = -32.25, 32.25, 0.5
		# self.min_snr_noise, self.max_snr_noise, self.delta_snr_noise = -3.225, 3.225, 0.05
		self.min_snr_noise, self.max_snr_noise, self.delta_snr_noise = -32.25, 32.25, 0.05
		self.neg_cut_lines = {}
		self.trash = []
		self.w = 0
		self.num_strips = numstrips
		self.cluster_size = clustersize
		self.noise_ana = noise_ana

		self.suffix = GetSuffixDictionary(self.trans_grid)

		self.noise_cuts = self.trans_grid.cuts_man.noise_cuts

		self.ph_cuts = self.trans_grid.cuts_man.GetPHCuts

		self.ph_ch_var = self.trans_grid.GetPHChVar

		self.phN_chs_var = self.trans_grid.GetPHNChsVar

		self.noise_varz = self.trans_grid.noise_varz

		self.strips_for_analysis = np.arange(self.cluster_size, dtype='i4')

	def PosCanvas(self, canvas_name):
		self.w = PositionCanvas(self.trans_grid, canvas_name, self.w, self.window_shift)

	def OverlayNoiseDistribution(self, histo, cells='all', isFriend=False):
		if self.noise_ana:
			self.noise_ana.OverlayNoiseDistribution(histo, cells, isFriend)

	def Do1DHistograms(self, cells='all', doLog=False, typ='adc', isFriend=False):
		"""
		Creates 1D histograms with the distribution of the PH of different channels ordered by proximity or PH magnitude. This is usefull to identify excess in negative charges, mostly in PH_CH1 and make a cut if necessary
		:param cells: Which cells to take into account for the profile
		:param doLog: Show plots in Logarithmic plot
		:param typ: indicates either to show PH in adc ('adc') or in sigmas ('snr')
		:param isFriend: if true, it will use the data from a pedTree friend
		:return:
		"""
		def DrawHisto(name, histo_limits, plot_lims, deltax, varz, varname, cuts):
			tempc = self.trans_grid.cuts_man.ConcatenateCutWithCells(cut=cuts, cells=cells)
			self.trans_grid.DrawHisto1D(name, histo_limits['min'], histo_limits['max'], deltax, varz, varname, tempc)

			neg_cut = self.trans_grid.neg_cut_adc if typ == 'adc' else self.trans_grid.neg_cut_snr
			self.trans_grid.line[name + '_neg_cut'] = ro.TLine(-neg_cut, 0, -neg_cut, self.trans_grid.histo[name].GetMaximum())
			self.trans_grid.line[name + '_neg_cut'].SetLineColor(ro.kBlue)
			self.trans_grid.line[name + '_neg_cut'].SetLineWidth(2)
			self.trans_grid.line[name + '_neg_cut'].SetLineStyle(2)
			self.trans_grid.line[name + '_neg_cut'].Draw('same')

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
			hist_limits = GetSymmetric1DLimits(minz, maxz, delta, 2) if ch == 0 else GetSymmetric1DLimits(-min(abs(minz), abs(maxz)), min(abs(minz), abs(maxz)), delta, 2)
			plot_limits = GetSymmetric1DLimits(minz, maxz, delta, 1, False) if ch == 0 else GetSymmetric1DLimits(-min(abs(minz), abs(maxz)), min(abs(minz), abs(maxz)), delta, 1)
			hname = 'PH_Ch{i}_{t}_{s}'.format(i=ch, s=suffix, t=typ.lower()) if not isFriend else 'PH_Ch{i}_buffer_{v}_{t}_{s}'.format(v=self.trans_grid.noise_friend_buffer, i=ch, s=suffix, t=typ.lower())
			DrawHisto(hname, hist_limits, plot_limits, delta, self.ph_ch_var(ch, 'Ch', typ=='snr', isFriend), 'PH cluster ch{i} [{t}]'.format(i=ch, t=typ.upper()), tempcuts)

			tempcuts = self.ph_cuts('PH_H{i}'.format(i=ch+1), isFriend)
			minz, maxz = self.minz['PH_H{i}_{t}'.format(i=ch+1, t=typ.lower())], self.maxz['PH_H{i}_{t}'.format(i=ch+1, t=typ.lower())]
			hist_limits = GetSymmetric1DLimits(minz, maxz, delta, 2) if ch == 0 else GetSymmetric1DLimits(-min(abs(minz), abs(maxz)), min(abs(minz), abs(maxz)), delta, 2)
			plot_limits = GetSymmetric1DLimits(minz, maxz, delta, 1, False) if ch == 0 else GetSymmetric1DLimits(-min(abs(minz), abs(maxz)), min(abs(minz), abs(maxz)), delta, 1)
			hname = 'PH_H{i}_{t}_{s}'.format(i=ch + 1, s=suffix, t=typ.lower()) if not isFriend else 'PH_H{i}_buffer_{v}_{t}_{s}'.format(v=self.trans_grid.noise_friend_buffer, i=ch+1, s=suffix, t=typ.lower())
			DrawHisto(hname, hist_limits, plot_limits, delta, self.ph_ch_var(ch+1, 'H', typ=='snr', isFriend), 'PH highest {i}{sf} ch [{t}]'.format(i=ch + 1, t=typ.upper(), sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), tempcuts)

	def Do2DProfileMaps(self, cells='all', isFriend=False):
		"""
		Makes profile map of the PH of different channels ordered by proximity or PH magnitude. The Y axis is the predicted hit position, while the X axis is either the clusterChannel{i} or the clusterChannelHighest{i}
		:param cells: Which cells to take into account for the profile
		:param isFriend: if true, it will use the data from a pedTree friend
		:return:
		"""
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
		"""
		Make a 2D histogram with the ph of a certain channel or group of channels vs the predicted hit position in the strip
		:param cells: Which cells to take into account for the profile
		:param typ: indicates either to show PH in adc ('adc') or in sigmas ('snr')
		:param isFriend: if true, it will use the data from a pedTree friend
		:return:
		"""
		minx, maxx, deltax, xname, xvar = -self.trans_grid.row_cell_info_diamond['width'] / (self.trans_grid.col_pitch * 2.0), self.trans_grid.row_cell_info_diamond['width'] / (self.trans_grid.col_pitch * 2.0), self.trans_grid.cell_resolution / float(self.trans_grid.row_cell_info_diamond['width']), 'dia pred. strip hit pos in strip', 'x0'

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
		"""
		Makes correlation plots between the different PH of different channels or group of channels ordered by proximity to hit position or by PH height
		:param cells: Which cells to take into account for the profile
		:param typ: indicates either to show PH in adc ('adc') or in sigmas ('snr')
		:param isFriend: if true, it will use the data from a pedTree friend
		:return:
		"""
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

	def DoStudyBla(self):
		deltax = 5
		xmin = -500
		xmax = 500
		ymin = -500
		ymax = 3000
		deltay = 50
		hlims1D = GetSymmetric1DLimits(xmin, xmax, deltax, 1)
		hlimsy2D = Get1DLimits(ymin, ymax, deltay)
		colorsopt = {0: ro.kBlack, 1: ro.kBlue, 2: ro.kRed, 3: ro.kOrange, 4: ro.kMagenta, 5: ro.kViolet, 6: ro.kGreen, 7: ro.kBlue+2, 8: ro.kRed+2, 9: ro.kOrange+7, 10: ro.kMagenta+2}

		tempch = {}
		hname = {}
		hnamesc = {}
		hnameshift = {}
		hbothname = {}
		hbothnameshift = {}
		for cells in ['bad', '', 'all', 'good']:
			tempch[cells] = {}
			hname[cells] = {}
			hnamesc[cells] = {}
			hnameshift[cells] = {}
			for ch in list(set(range(1, self.num_strips) + [self.cluster_size - 1])):
				tempc = self.trans_grid.cuts_man.AndCuts([self.ph_cuts('PH_Ch{c}'.format(c=ch), False), self.trans_grid.cuts_man.cluster_ch_lowest[ch]])
				tempch[cells][ch] = self.trans_grid.cuts_man.ConcatenateCutWithCells(cut=tempc, cells=cells)
				hname[cells][ch] = 'PH_Ch{c}_with_Ch{c}_lowest_{s}'.format(c=ch, s=self.suffix[cells])
				hnamesc[cells][ch] = 'PH_Ch{c}_with_Ch{c}_lowest_scaled_{s}'.format(c=ch, s=self.suffix[cells])
				hnameshift[cells][ch] = 'PH_Ch{c}_with_Ch{c}_lowest_shifted_{s}'.format(c=ch, s=self.suffix[cells])
				self.trans_grid.DrawHisto1D(hname[cells][ch], hlims1D['min'], hlims1D['max'], deltax, self.ph_ch_var(ch, 'Ch', False, False), 'PH cluster ch{c} [ADC]'.format(c=ch), tempch[cells][ch])
				self.PosCanvas(hname[cells][ch])
				self.trans_grid.histo[hnamesc[cells][ch]] = self.trans_grid.histo[hname[cells][ch]].Clone('h_' + hnamesc[cells][ch])
				self.trans_grid.histo[hnamesc[cells][ch]].SetTitle('h_' + hnamesc[cells][ch])
				self.trans_grid.histo[hnamesc[cells][ch]].SetLineColor(colorsopt[ch])
				self.trans_grid.histo[hnameshift[cells][ch]] = self.trans_grid.histo[hnamesc[cells][ch]].Clone('h_' + hnameshift[cells][ch])
				self.trans_grid.histo[hnameshift[cells][ch]].SetStats(0)
				self.trans_grid.histo[hnamesc[cells][ch]].Scale(1.0 / self.trans_grid.histo[hnamesc[cells][ch]].GetMaximum())

                        if self.cluster_size > 2:
			        hbothname[cells] = 'PH_Ch1_and_PH_Ch2_lowest_overlaid_and_scaled_{s}'.format(s=self.suffix[cells])
			        self.trans_grid.canvas[hbothname[cells]] = ro.TCanvas('c_' + hbothname[cells], 'c_' + hbothname[cells], 1)
			        self.trans_grid.histo[hnamesc[cells][2]].SetTitle(hbothname[cells])
        			self.trans_grid.histo[hnamesc[cells][2]].Draw()
	        		# self.trans_grid.histo[hnamesc[1]].Scale(self.trans_grid.histo[hnamesc[2]].GetBinContent(self.trans_grid.histo[hnamesc[1]].GetMaximumBin()) / self.trans_grid.histo[hnamesc[1]].GetMaximum())
		        	self.trans_grid.histo[hnamesc[cells][1]].SetTitle(hbothname[cells])
			        self.trans_grid.histo[hnamesc[cells][1]].Draw('same')
			        SetDefault1DCanvasSettings(self.trans_grid.canvas[hbothname[cells]])
			        self.trans_grid.canvas[hbothname[cells]].SetLogy()
			        self.PosCanvas(hbothname[cells])

			        hbothnameshift[cells] = 'PH_Ch1_shifted_and_PH_Ch2_scaled_lowest_overalid_{s}'.format(s=self.suffix[cells])
			        maxbinch1 = self.trans_grid.histo[hnameshift[cells][1]].GetMaximumBin()
			        maxbinch2 = self.trans_grid.histo[hnameshift[cells][2]].GetMaximumBin()
			        if maxbinch2 > maxbinch1:
				        for i in np.arange(int(self.trans_grid.histo[hnameshift[cells][1]].GetNbinsX()), int(maxbinch2 - maxbinch1), -1):
					        self.trans_grid.histo[hnameshift[cells][1]].SetBinContent(i, self.trans_grid.histo[hnameshift[cells][1]].GetBinContent(i - int(maxbinch2 - maxbinch1)))
				        for i in np.arange(1, int(maxbinch2 - maxbinch1) + 1):
					        self.trans_grid.histo[hnameshift[cells][1]].SetBinContent(i, 0)
			        elif maxbinch1 > maxbinch2:
				        for i in np.arange(1, int(self.trans_grid.histo[hnameshift[cells][1]].GetNbinsX() - int(maxbinch1 - maxbinch2)) + 1, 1):
					        self.trans_grid.histo[hnameshift[cells][1]].SetBinContent(i, self.trans_grid.histo[hnameshift[cells][1]].GetBinContent(i + int(maxbinch1 - maxbinch2)))
				        for i in np.arange(int(self.trans_grid.histo[hnameshift[cells][1]].GetNbinsX()) - int(maxbinch1 - maxbinch2) + 1, int(self.trans_grid.histo[hnameshift[cells][1]].GetNbinsX()) + 1, 1):
					        self.trans_grid.histo[hnameshift[cells][1]].SetBinContent(i, 0)
			        self.trans_grid.histo[hnameshift[cells][2]].Scale(self.trans_grid.histo[hnameshift[cells][1]].GetMaximum() / self.trans_grid.histo[hnameshift[cells][2]].GetMaximum())
			        self.trans_grid.canvas[hbothnameshift[cells]] = ro.TCanvas('c_' + hbothnameshift[cells], 'c_' + hbothnameshift[cells], 1)
			        self.trans_grid.histo[hnameshift[cells][1]].SetTitle('PH_Ch{c}_with_Ch{c}_lowest_shifted_{s}'.format(c=1, s=self.suffix[cells]))
			        self.trans_grid.histo[hnameshift[cells][2]].SetTitle('PH_Ch{c}_with_Ch{c}_lowest_scaled_{s}'.format(c=2, s=self.suffix[cells]))
			        self.trans_grid.histo[hnameshift[cells][1]].Draw()
			        self.trans_grid.histo[hnameshift[cells][2]].Draw('same')
			        legend = self.trans_grid.canvas[hbothnameshift[cells]].BuildLegend()
			        self.trans_grid.canvas[hbothnameshift[cells]].SetLogy()
			        SetDefault1DCanvasSettings(self.trans_grid.canvas[hbothnameshift[cells]])
			        self.PosCanvas(hbothnameshift[cells])

			satOpts = [True, False]
			for satOpt in satOpts:
				for ch in list(set(range(1, self.num_strips) + [self.cluster_size - 1])):
					tempc = tempch[cells][ch] if satOpt else self.trans_grid.cuts_man.AndCuts([tempch[cells][ch], '(diaChADC[clusterChannelHighest1]!=4095)'])
					hname1 = 'PH_H1_vs_PH_Ch{c}_lowest_{s}'.format(s=self.suffix[cells], c=ch) if satOpt else 'PH_H1_vs_PH_Ch{c}_lowest_no_sat_{s}'.format(s=self.suffix[cells], c=ch)
					self.trans_grid.DrawHisto2D(hname1, hlims1D['min'], hlims1D['max'], deltax, 'PH cluster ch{c} when lowest [ADC]'.format(c=ch), hlimsy2D['min'], hlimsy2D['max'], deltay, 'PH highest ch [ADC]', self.ph_ch_var(ch, 'Ch', False, False), self.ph_ch_var(1, 'H', False, False), tempc)
					self.PosCanvas(hname1)
                                        if self.num_strips > 1:
					        hname1 = 'PH2_H_vs_PH_Ch{c}_lowest_{s}'.format(s=self.suffix[cells], c=ch) if satOpt else 'PH2_H_vs_PH_Ch{c}_lowest_no_sat_{s}'.format(s=self.suffix[cells], c=ch)
					        self.trans_grid.DrawHisto2D(hname1, hlims1D['min'], hlims1D['max'], deltax, 'PH cluster ch{c} when lowest [ADC]'.format(c=ch), hlimsy2D['min'], hlimsy2D['max'], deltay, 'PH2 highest chs [ADC]', self.ph_ch_var(ch, 'Ch', False, False), self.phN_chs_var(2, 'H', False, False), tempc)
					        self.PosCanvas(hname1)
					deltax2 = 0.01
					deltay2 = 5
					hlims1V2 = GetSymmetric1DLimits(-1, 1, deltax2)
					hlimsyV2 = GetSymmetric1DLimits(-500, 500, deltay2)
					tempc = tempch[cells][ch] if satOpt else self.trans_grid.cuts_man.AndCuts([tempch[cells][ch], '(diaChADC[clusterChannelHighest1]!=4095)'])
					tempc = self.trans_grid.cuts_man.ConcatenateCutWithCells(cut=tempc, cells=cells)
					hname1 = 'PH_Ch{c}_Vs_Ratio_PH_Ch{c}_PH_Ch0_with_PH_Ch{c}_lowest_with_sat_{s}'.format(c=ch, s=self.suffix[cells]) if satOpt else 'PH_Ch{c}_Vs_Ratio_PH_Ch{c}_PH_Ch0_with_PH_Ch{c}_lowest_no_sat_{s}'.format(c=ch, s=self.suffix[cells])
					self.trans_grid.DrawHisto2D(hname1, hlims1V2['min'], hlims1V2['max'], deltax2, 'PH_Ch{c}/PH_Ch0'.format(c=ch), hlimsyV2['min'], hlimsyV2['max'], deltay2, 'PH_Ch'+str(ch), 'diaChSignal[clusterChannel{c}]/diaChSignal[clusterChannel0]'.format(c=ch), 'diaChSignal[clusterChannel{c}]'.format(c=ch), tempc)
					SetDefault1DCanvasSettings(self.trans_grid.canvas[hname1])
					self.trans_grid.canvas[hname1].SetLogz()
					# self.trans_grid.FitPol(hname1, 1, 0, True)
					self.PosCanvas(hname1)



	def DoClusterStudies(self, cells='all', typ='adc', isFriend=False):
		self.SetStripsForAnalysis()
		self.Do2DProfileMaps(cells)
		self.DoStrips2DHistograms(cells, typ=typ, isFriend=isFriend)
		self.DoPHStripCorrelations(cells, typ=typ, isFriend=isFriend)
		self.Do1DHistograms(cells, False, typ=typ, isFriend=isFriend)
		self.DoStudyBla()
		# self.Do1DHistograms(cells, True, typ='adc')

	def SetStripsForAnalysis(self, arr=None):
		if arr:
			self.strips_for_analysis = arr
		else:
			self.strips_for_analysis = np.arange(self.num_strips, dtype='i4')
			self.strips_for_analysis = self.strips_for_analysis if self.cluster_size > 3 else np.unique(np.append(self.strips_for_analysis, self.cluster_size - 1))
		self.minz = self.trans_grid.minz
		self.maxz = self.trans_grid.maxz

if __name__ == '__main__':
	c = ClusterChannelsAnalysis(None, 0, 0)
