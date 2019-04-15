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
		self.minz = {t: {} for t in ['all', 'good', 'bad']}
		self.maxz = {t: {} for t in ['all', 'good', 'bad']}
		self.neg_cut_lines = {}
		self.trash = []
		self.w = 0
		self.trans_grid = trans_grid
		self.num_strips = numstrips
		self.cluster_size = clustersize
		self.noise_ana = noise_ana

		self.suffix = {'all': 'all', 'good': 'selection', 'bad': 'not_selection'}

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

	def PosCanvas(self, canvas_name):
		self.w = PositionCanvas(self.trans_grid, canvas_name, self.w, self.window_shift)

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
		self.PosCanvas('signal_noise_{c}_snr'.format(c=suffix))
		self.trans_grid.DrawHisto1D('signal_noise_{c}_adc'.format(c=suffix), self.min_adc_noise, self.max_adc_noise, self.delta_adc_noise, self.noise_varz['adc'], varname='Signal not in cluster (SNR)', cuts=temp_cut_noise, option='e hist')
		self.trans_grid.FitGaus('signal_noise_{c}_adc'.format(c=suffix))
		self.trans_grid.histo['signal_noise_{c}_adc'.format(c=suffix)].GetXaxis().SetRangeUser(-32, 32)
		self.PosCanvas('signal_noise_{c}_adc'.format(c=suffix))

	def OverlayNoiseDistribution(self, histo, cells='all'):
		if self.noise_ana:
			self.noise_ana.OverlayNoiseDistribution(histo, cells)

	def Do1DHistograms(self, cells='all', doLog=False, typ='adc'):
		def DrawHisto(name, histo_limits, plot_lims, deltax, varz, varname, cuts):
			self.trans_grid.DrawHisto1D(name, histo_limits['min'], histo_limits['max'], deltax, varz, varname, cuts)
			self.trans_grid.histo[name].GetXaxis().SetRangeUser(plot_lims['min'], plot_lims['max'])
			SetX1X2NDC(self.trans_grid.histo[name], 0.15, 0.45, 'stats')
			self.OverlayNoiseDistribution(self.trans_grid.histo[name], cells)
			if doLog: self.trans_grid.canvas[name].SetLogy()
			legend = self.trans_grid.canvas[name].BuildLegend()
			ro.gPad.Update()
			SetLegendX1X2Y1Y2(legend, 0.15, 0.45, 0.5, 0.6)
			self.PosCanvas(name)

		suffix = self.suffix[cells] if cells in self.suffix.keys() else ''
		suffix = suffix + '_logScale' if doLog else suffix

		for ch in xrange(self.cluster_size):
			if typ == 'adc':
				if 'PH_Ch' + str(ch) in self.ph_adc_ch_varz.keys():
					tempcuts = self.ph_adc_ch_cuts[cells]['PH_Ch{i}'.format(i=ch)]
					minz, maxz = self.minz[cells]['PH_Ch{i}_adc'.format(i=ch)], self.maxz[cells]['PH_Ch{i}_adc'.format(i=ch)]
					hist_limits = GetSymmetric1DLimits(minz, maxz, self.delta_adc, 2)
					plot_limits = GetSymmetric1DLimits(minz, maxz, self.delta_adc, 1, False)
					DrawHisto('PH_Ch{i}_adc_{s}'.format(i=ch, s=suffix), hist_limits, plot_limits, self.delta_adc, self.ph_adc_ch_varz['PH_Ch{i}'.format(i=ch)], 'PH cluster ch{i} [SNR]'.format(i=ch), tempcuts)

				if 'PH_H' + str(ch + 1) in self.ph_adc_h_varz.keys():
					tempcuts = self.ph_adc_h_cuts[cells]['PH_H{i}'.format(i=ch+1)]
					minz, maxz = self.minz[cells]['PH_H{i}_adc'.format(i=ch+1)], self.maxz[cells]['PH_H{i}_adc'.format(i=ch+1)]
					hist_limits = GetSymmetric1DLimits(minz, maxz, self.delta_adc, 2)
					plot_limits = GetSymmetric1DLimits(minz, maxz, self.delta_adc, 1, False)
					DrawHisto('PH_H{i}_adc_{s}'.format(i=ch + 1, s=suffix), hist_limits, plot_limits, self.delta_adc, self.ph_adc_h_varz['PH_H{i}'.format(i=ch + 1)], 'PH highest {i}{sf} ch [SNR]'.format(i=ch + 1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), tempcuts)
			else:
				if 'PH_Ch' + str(ch) in self.ph_snr_ch_varz.keys():
					tempcuts = self.ph_snr_ch_cuts[cells]['PH_Ch{i}'.format(i=ch)]
					minz, maxz = self.minz[cells]['PH_Ch{i}_snr'.format(i=ch)], self.maxz[cells]['PH_Ch{i}_snr'.format(i=ch)]
					hist_limits = GetSymmetric1DLimits(minz, maxz, self.delta_snr, 2)
					plot_limits = GetSymmetric1DLimits(minz, maxz, self.delta_snr, 1, False)
					DrawHisto('PH_Ch{i}_snr_{s}'.format(i=ch, s=suffix), hist_limits, plot_limits, self.delta_snr, self.ph_snr_ch_varz['PH_Ch{i}'.format(i=ch)], 'PH cluster ch{i} [SNR]'.format(i=ch), tempcuts)

				if 'PH_H' + str(ch + 1) in self.ph_snr_h_varz.keys():
					tempcuts = self.ph_snr_h_cuts[cells]['PH_H{i}'.format(i=ch+1)]
					minz, maxz = self.minz[cells]['PH_H{i}_snr'.format(i=ch+1)], self.maxz[cells]['PH_H{i}_snr'.format(i=ch+1)]
					hist_limits = GetSymmetric1DLimits(minz, maxz, self.delta_snr, 2)
					plot_limits = GetSymmetric1DLimits(minz, maxz, self.delta_snr, 1, False)
					DrawHisto('PH_H{i}_snr_{s}'.format(i=ch + 1, s=suffix), hist_limits, plot_limits, self.delta_snr, self.ph_snr_h_varz['PH_H{i}'.format(i=ch + 1)], 'PH highest {i}{sf} ch [SNR]'.format(i=ch + 1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), tempcuts)

	def Do2DProfileMaps(self, cells='all'):
		def DrawProfile(name, varx, xname, varz, cuts, zmax, zmin):
			self.trans_grid.DrawProfile2DDiamondChannel(name, varx, xname, varz, cuts)
			self.trans_grid.profile[name].SetMaximum(zmax)
			self.trans_grid.profile[name].SetMinimum(zmin)
			self.trans_grid.DrawTCutGs(name, 'diamond')
			self.PosCanvas(name)

		suffix = self.suffix[cells] if cells in self.suffix.keys() else ''
		for ch in xrange(self.cluster_size):
			if 'PH_Ch' + str(ch) in self.ph_adc_ch_varz.keys():
				tempcuts = self.ph_adc_ch_cuts[cells]['PH_Ch{i}'.format(i=ch)]
				minz, maxz = self.minz[cells]['PH_Ch{i}_adc'.format(i=ch)], self.maxz[cells]['PH_Ch{i}_adc'.format(i=ch)]
				DrawProfile('PH_Ch{c}_map_{s}'.format(c=ch, s=suffix), 'clusterChannel{c}'.format(c=ch), 'cluster ch{c}'.format(c=ch), self.ph_adc_ch_varz['PH_Ch{i}'.format(i=ch)], tempcuts, maxz, minz)

			if 'PH_H' + str(ch + 1) in self.ph_adc_h_varz.keys():
				tempcuts = self.ph_snr_h_cuts[cells]['PH_H{i}'.format(i=ch+1)]
				minz, maxz = self.minz[cells]['PH_H{i}_adc'.format(i=ch+1)], self.maxz[cells]['PH_H{i}_adc'.format(i=ch+1)]
				DrawProfile('PH_H{c}_map_{s}'.format(c=ch+1, s=suffix), 'clusterChannelHighest{c}'.format(c=ch+1), 'highest {c}{sf} ch'.format(c=ch+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), self.ph_adc_h_varz['PH_H{i}'.format(i=ch+1)], tempcuts, maxz, minz)

	def DoStrips2DHistograms(self, cells='all', typ='adc'):
		minx, maxx, deltax, xname, xvar = -0.5, 0.5, self.trans_grid.cell_resolution / float(self.trans_grid.row_info_diamond['pitch']), 'dia pred. strip hit pos', 'diaChXPred-TMath::Floor(diaChXPred+0.5)'

		def DrawHistogram(name, zmin, zmax, yname, yvar, cuts, typ='adc'):
			histo_limits = Get1DLimits(zmin, zmax, 4 * self.delta_adc) if typ == 'adc' else Get1DLimits(zmin, zmax, 4 * self.delta_snr)
			deltay = 4 * self.delta_adc if typ == 'adc' else 4 * self.delta_snr
			self.trans_grid.DrawHisto2D(name, minx, maxx, deltax, xname, min(0, histo_limits['min']), histo_limits['max'], deltay, yname, xvar, yvar, cuts)
			self.PosCanvas(name)

		suffix = self.suffix[cells] if cells in self.suffix.keys() else ''
		for ch in xrange(self.cluster_size):
			if typ == 'adc':
				if 'PH_Ch' + str(ch) in self.ph_adc_ch_varz.keys():
					tempcuts = self.ph_adc_ch_cuts[cells]['PH_Ch{i}'.format(i=ch)]
					minz, maxz = self.minz[cells]['PH_Ch{c}_adc'.format(c=ch)], self.maxz[cells]['PH_Ch{c}_adc'.format(c=ch)]
					DrawHistogram('PH_Ch{c}_Vs_strip_location_adc_{s}'.format(c=ch, s=suffix), minz, maxz, 'PH cluster ch{c} [SNR]'.format(c=ch), self.ph_adc_ch_varz['PH_Ch{i}'.format(i=ch)], tempcuts)

				if 'PH_H' + str(ch+1) in self.ph_adc_h_varz.keys():
					tempcuts = self.ph_adc_h_cuts[cells]['PH_H{i}'.format(i=ch+1)]
					minz, maxz = self.minz[cells]['PH_H{c}_adc'.format(c=ch+1)], self.maxz[cells]['PH_H{c}_adc'.format(c=ch+1)]
					DrawHistogram('PH_H{c}_Vs_strip_location_adc_{s}'.format(c=ch+1, s=suffix), minz, maxz, 'PH highest {c}{sf} ch [SNR]'.format(c=ch+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), self.ph_adc_h_varz['PH_H{i}'.format(i=ch+1)], tempcuts)

				if 'PH{ch}_Ch'.format(ch=ch+1) in self.phN_adc_ch_varz.keys():
					tempcuts = self.phN_adc_ch_cuts[cells]['PH{i}_Ch'.format(i=ch+1)]
					minz, maxz = self.minz[cells]['PH{c}_Ch_adc'.format(c=ch+1)], self.maxz[cells]['PH{c}_Ch_adc'.format(c=ch+1)]
					DrawHistogram('PH{c}_Ch_Vs_strip_location_adc_{s}'.format(c=ch+1, s=suffix), minz, maxz, 'PH{c} cluster chs [SNR]'.format(c=ch+1), self.phN_adc_ch_varz['PH{i}_Ch'.format(i=ch+1)], tempcuts)

				if ch != self.cluster_size - 1:
					# these will be the same as PH{c1}_Ch when c1 == self.cluster_size - 1, because it is all the channels in the cluster
					if 'PH{ch}_H'.format(ch=ch+1) in self.phN_adc_h_varz.keys():
						tempcuts = self.phN_adc_h_cuts[cells]['PH{i}_H'.format(i=ch+1)]
						minz, maxz = self.minz[cells]['PH{c}_H_adc'.format(c=ch+1)], self.maxz[cells]['PH{c}_H_adc'.format(c=ch+1)]
						DrawHistogram('PH{c}_H_Vs_strip_location_adc_{s}'.format(c=ch+1, s=suffix), minz, maxz, 'PH{c} highest chs [SNR]'.format(c=ch + 1), self.phN_adc_h_varz['PH{i}_H'.format(i=ch + 1)], tempcuts)
			else:
				if 'PH_Ch' + str(ch) in self.ph_snr_ch_varz.keys():
					tempcuts = self.ph_snr_ch_cuts[cells]['PH_Ch{i}'.format(i=ch)]
					minz, maxz = self.minz[cells]['PH_Ch{c}_snr'.format(c=ch)], self.maxz[cells]['PH_Ch{c}_snr'.format(c=ch)]
					DrawHistogram('PH_Ch{c}_Vs_strip_location_snr_{s}'.format(c=ch, s=suffix), minz, maxz, 'PH cluster ch{c} [SNR]'.format(c=ch), self.ph_snr_ch_varz['PH_Ch{i}'.format(i=ch)], tempcuts, 'snr')

				if 'PH_H' + str(ch+1) in self.ph_snr_h_varz.keys():
					tempcuts = self.ph_snr_h_cuts[cells]['PH_H{i}'.format(i=ch+1)]
					minz, maxz = self.minz[cells]['PH_H{c}_snr'.format(c=ch+1)], self.maxz[cells]['PH_H{c}_snr'.format(c=ch+1)]
					DrawHistogram('PH_H{c}_Vs_strip_location_snr_{s}'.format(c=ch+1, s=suffix), minz, maxz, 'PH highest {c}{sf} ch [SNR]'.format(c=ch+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), self.ph_snr_h_varz['PH_H{i}'.format(i=ch+1)], tempcuts, 'snr')

				if 'PH{ch}_Ch'.format(ch=ch+1) in self.phN_snr_ch_varz.keys():
					tempcuts = self.phN_snr_ch_cuts[cells]['PH{i}_Ch'.format(i=ch+1)]
					minz, maxz = self.minz[cells]['PH{c}_Ch_snr'.format(c=ch+1)], self.maxz[cells]['PH{c}_Ch_snr'.format(c=ch+1)]
					DrawHistogram('PH{c}_Ch_Vs_strip_location_snr_{s}'.format(c=ch+1, s=suffix), minz, maxz, 'PH{c} cluster chs [SNR]'.format(c=ch+1), self.phN_snr_ch_varz['PH{i}_Ch'.format(i=ch+1)], tempcuts, 'snr')

				if ch != self.cluster_size - 1:
					# these will be the same as PH{c1}_Ch when c1 == self.cluster_size - 1, because it is all the channels in the cluster
					if 'PH{ch}_H'.format(ch=ch+1) in self.phN_snr_h_varz.keys():
						tempcuts = self.phN_snr_h_cuts[cells]['PH{i}_H'.format(i=ch+1)]
						minz, maxz = self.minz[cells]['PH{c}_H_snr'.format(c=ch+1)], self.maxz[cells]['PH{c}_H_snr'.format(c=ch+1)]
						DrawHistogram('PH{c}_H_Vs_strip_location_snr_{s}'.format(c=ch+1, s=suffix), minz, maxz, 'PH{c} highest chs [SNR]'.format(c=ch + 1), self.phN_snr_h_varz['PH{i}_H'.format(i=ch + 1)], tempcuts, 'snr')

	def DoPHStripCorrelations(self, cells='all', typ='adc'):
		xmin_adc, xmax_adc = -550, 2650
		xmin_snr, xmax_snr = -55, 265

		def DrawHistogram(name, xname, ymin, ymax, yname, varx, vary, cuts, typ='adc'):
			if typ == 'adc':
				self.trans_grid.DrawHisto2D(name, xmin_adc, xmax_adc, self.delta_adc * 4, xname, ymin, ymax, self.delta_adc * 4, yname, varx, vary, cuts)
			else:
				self.trans_grid.DrawHisto2D(name, xmin_snr, xmax_snr, self.delta_snr * 4, xname, ymin, ymax, self.delta_snr * 4, yname, varx, vary, cuts)
			self.PosCanvas(name)

		for ch in xrange(self.cluster_size):
			for ch2 in xrange(self.cluster_size):
				if ch != self.cluster_size - 1:
					# these will be the same as PH{c1}_Ch when c1 == self.cluster_size - 1, because it is all the channels in the cluster
					if typ == 'adc':
						if 'PH{c1}_H'.format(c1=ch+1) in self.phN_adc_h_varz.keys() and 'PH_H{c2}'.format(c2=ch2+1) in self.ph_adc_h_varz.keys():
							tempcuts = self.trans_grid.cuts_man.AndCuts([self.phN_adc_h_cuts[cells]['PH{c1}_H'.format(c1=ch+1)], self.ph_adc_h_cuts[cells]['PH_H{c2}'.format(c2=ch2+1)]])
							# tempcuts = self.trans_grid.cuts_man.ConcatenateCuts(self.phN_adc_h_cuts[cells]['PH{c1}_H'.format(c1=ch+1)], self.ph_adc_h_cuts[cells]['PH_H{c2}'.format(c2=ch2+1)])
							DrawHistogram('PH{c1}_H_Vs_PH_H{c2}_adc'.format(c1=ch+1, c2=ch2+1), 'PH highest {c}{sf} ch [SNR]'.format(c=ch2+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), 0, 4200, 'PH{c1} highest chs [SNR]'.format(c1=ch+1), self.ph_adc_h_varz['PH_H{i}'.format(i=ch2+1)], self.phN_adc_h_varz['PH{i}_H'.format(i=ch+1)], tempcuts)

						if 'PH{c1}_H'.format(c1=ch+1) in self.phN_adc_h_varz.keys() and 'PH_Ch{c2}'.format(c2=ch2) in self.ph_adc_ch_varz.keys():
							tempcuts = self.trans_grid.cuts_man.AndCuts([self.phN_adc_h_cuts[cells]['PH{c1}_H'.format(c1=ch+1)], self.ph_adc_ch_cuts[cells]['PH_Ch{c2}'.format(c2=ch2)]])
							# tempcuts = self.trans_grid.cuts_man.ConcatenateCuts(self.phN_adc_h_cuts[cells]['PH{c1}_H'.format(c1=ch+1)], self.ph_adc_ch_cuts[cells]['PH_Ch{c2}'.format(c2=ch2)])
							DrawHistogram('PH{c1}_H_Vs_PH_Ch{c2}_adc'.format(c1=ch+1, c2=ch2), 'PH cluster ch{c2} [SNR]'.format(c2=ch2), 0, 4200, 'PH{c1} highest chs [SNR]'.format(c1=ch+1), self.ph_adc_ch_varz['PH_Ch{i}'.format(i=ch2)], self.phN_adc_h_varz['PH{i}_H'.format(i=ch+1)], tempcuts)
					else:
						if 'PH{c1}_H'.format(c1=ch+1) in self.phN_snr_h_varz.keys() and 'PH_H{c2}'.format(c2=ch2+1) in self.ph_snr_h_varz.keys():
							tempcuts = self.trans_grid.cuts_man.AndCuts([self.phN_snr_h_cuts[cells]['PH{c1}_H'.format(c1=ch+1)], self.ph_snr_h_cuts[cells]['PH_H{c2}'.format(c2=ch2+1)]])
							# tempcuts = self.trans_grid.cuts_man.ConcatenateCuts(self.phN_snr_h_cuts[cells]['PH{c1}_H'.format(c1=ch+1)], self.ph_snr_h_cuts[cells]['PH_H{c2}'.format(c2=ch2+1)])
							DrawHistogram('PH{c1}_H_Vs_PH_H{c2}_snr'.format(c1=ch+1, c2=ch2+1), 'PH highest {c}{sf} ch [SNR]'.format(c=ch2+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), 0, 420, 'PH{c1} highest chs [SNR]'.format(c1=ch+1), self.ph_snr_h_varz['PH_H{i}'.format(i=ch2+1)], self.phN_snr_h_varz['PH{i}_H'.format(i=ch+1)], tempcuts, 'snr')

						if 'PH{c1}_H'.format(c1=ch+1) in self.phN_snr_h_varz.keys() and 'PH_Ch{c2}'.format(c2=ch2) in self.ph_snr_ch_varz.keys():
							tempcuts = self.trans_grid.cuts_man.AndCuts([self.phN_snr_h_cuts[cells]['PH{c1}_H'.format(c1=ch+1)], self.ph_snr_ch_cuts[cells]['PH_Ch{c2}'.format(c2=ch2)]])
							# tempcuts = self.trans_grid.cuts_man.ConcatenateCuts(self.phN_snr_h_cuts[cells]['PH{c1}_H'.format(c1=ch+1)], self.ph_snr_ch_cuts[cells]['PH_Ch{c2}'.format(c2=ch2)])
							DrawHistogram('PH{c1}_H_Vs_PH_Ch{c2}_snr'.format(c1=ch+1, c2=ch2), 'PH cluster ch{c2} [SNR]'.format(c2=ch2), 0, 420, 'PH{c1} highest chs [SNR]'.format(c1=ch+1), self.ph_snr_ch_varz['PH_Ch{i}'.format(i=ch2)], self.phN_snr_h_varz['PH{i}_H'.format(i=ch+1)], tempcuts, 'snr')
				if typ == 'adc':
					if 'PH{c1}_Ch'.format(c1=ch+1) in self.phN_adc_ch_varz.keys() and 'PH_Ch{c2}'.format(c2=ch2) in self.ph_adc_ch_varz.keys():
						tempcuts = self.trans_grid.cuts_man.AndCuts([self.phN_adc_ch_cuts[cells]['PH{c1}_Ch'.format(c1=ch+1)], self.ph_adc_ch_cuts[cells]['PH_Ch{c2}'.format(c2=ch2)]])
						# tempcuts = self.trans_grid.cuts_man.ConcatenateCuts(self.phN_adc_ch_cuts[cells]['PH{c1}_Ch'.format(c1=ch+1)], self.ph_adc_ch_cuts[cells]['PH_Ch{c2}'.format(c2=ch2)])
						DrawHistogram('PH{c1}_Ch_Vs_PH_Ch{c2}_adc'.format(c1=ch+1, c2=ch2), 'PH cluster ch{c2} [SNR]'.format(c2=ch2), 0, 4200, 'PH{c1} cluster chs [SNR]'.format(c1=ch+1), self.ph_adc_ch_varz['PH_Ch{i}'.format(i=ch2)], self.phN_adc_ch_varz['PH{i}_Ch'.format(i=ch+1)], tempcuts)

					if 'PH{c1}_Ch'.format(c1=ch+1) in self.phN_adc_ch_varz.keys() and 'PH_H{c2}'.format(c2=ch2+1) in self.ph_adc_h_varz.keys():
						tempcuts = self.trans_grid.cuts_man.AndCuts([self.phN_adc_ch_cuts[cells]['PH{c1}_Ch'.format(c1=ch+1)], self.ph_adc_h_cuts[cells]['PH_H{c2}'.format(c2=ch2+1)]])
						# tempcuts = self.trans_grid.cuts_man.ConcatenateCuts(self.phN_adc_ch_cuts[cells]['PH{c1}_Ch'.format(c1=ch+1)], self.ph_adc_h_cuts[cells]['PH_H{c2}'.format(c2=ch2+1)])
						DrawHistogram('PH{c1}_Ch_Vs_PH_H{c2}_adc'.format(c1=ch+1, c2=ch2+1), 'PH highest {c2}{sf} ch [SNR]'.format(c2=ch2+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), 0, 4200, 'PH{c1} cluster chs [SNR]'.format(c1=ch+1), self.ph_adc_h_varz['PH_H{i}'.format(i=ch2+1)], self.phN_adc_ch_varz['PH{i}_Ch'.format(i=ch+1)], tempcuts)

					if 'PH_H{c1}'.format(c1=ch+1) in self.ph_adc_h_varz.keys() and 'PH_Ch{c2}'.format(c2=ch2) in self.ph_adc_ch_varz.keys():
						tempcuts = self.trans_grid.cuts_man.AndCuts([self.ph_adc_h_cuts[cells]['PH_H{c1}'.format(c1=ch+1)], self.ph_adc_ch_cuts[cells]['PH_Ch{c2}'.format(c2=ch2)]])
						# tempcuts = self.trans_grid.cuts_man.ConcatenateCuts(self.ph_adc_h_cuts[cells]['PH_H{c1}'.format(c1=ch+1)], self.ph_adc_ch_cuts[cells]['PH_Ch{c2}'.format(c2=ch2)])
						DrawHistogram('PH_H{c1}_Vs_PH_Ch{c2}_adc'.format(c1=ch+1, c2=ch2), 'PH cluster ch{c2} [SNR]'.format(c2=ch2), -550, 2650, 'PH highest {c1}{sf} ch [SNR]'.format(c1=ch+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), self.ph_adc_ch_varz['PH_Ch{i}'.format(i=ch2)], self.ph_adc_h_varz['PH_H{i}'.format(i=ch+1)], tempcuts)
				else:
					if 'PH{c1}_Ch'.format(c1=ch+1) in self.phN_snr_ch_varz.keys() and 'PH_Ch{c2}'.format(c2=ch2) in self.ph_snr_ch_varz.keys():
						tempcuts = self.trans_grid.cuts_man.AndCuts([self.phN_snr_ch_cuts[cells]['PH{c1}_Ch'.format(c1=ch+1)], self.ph_snr_ch_cuts[cells]['PH_Ch{c2}'.format(c2=ch2)]])
						# tempcuts = self.trans_grid.cuts_man.ConcatenateCuts(self.phN_snr_ch_cuts[cells]['PH{c1}_Ch'.format(c1=ch+1)], self.ph_snr_ch_cuts[cells]['PH_Ch{c2}'.format(c2=ch2)])
						DrawHistogram('PH{c1}_Ch_Vs_PH_Ch{c2}_snr'.format(c1=ch+1, c2=ch2), 'PH cluster ch{c2} [SNR]'.format(c2=ch2), 0, 420, 'PH{c1} cluster chs [SNR]'.format(c1=ch+1), self.ph_snr_ch_varz['PH_Ch{i}'.format(i=ch2)], self.phN_snr_ch_varz['PH{i}_Ch'.format(i=ch+1)], tempcuts, 'snr')

					if 'PH{c1}_Ch'.format(c1=ch+1) in self.phN_snr_ch_varz.keys() and 'PH_H{c2}'.format(c2=ch2+1) in self.ph_snr_h_varz.keys():
						tempcuts = self.trans_grid.cuts_man.AndCuts([self.phN_snr_ch_cuts[cells]['PH{c1}_Ch'.format(c1=ch+1)], self.ph_snr_h_cuts[cells]['PH_H{c2}'.format(c2=ch2+1)]])
						# tempcuts = self.trans_grid.cuts_man.ConcatenateCuts(self.phN_snr_ch_cuts[cells]['PH{c1}_Ch'.format(c1=ch+1)], self.ph_snr_h_cuts[cells]['PH_H{c2}'.format(c2=ch2+1)])
						DrawHistogram('PH{c1}_Ch_Vs_PH_H{c2}_snr'.format(c1=ch+1, c2=ch2+1), 'PH highest {c2}{sf} ch [SNR]'.format(c2=ch2+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), 0, 420, 'PH{c1} cluster chs [SNR]'.format(c1=ch+1), self.ph_snr_h_varz['PH_H{i}'.format(i=ch2+1)], self.phN_snr_ch_varz['PH{i}_Ch'.format(i=ch+1)], tempcuts, 'snr')

					if 'PH_H{c1}'.format(c1=ch+1) in self.ph_snr_h_varz.keys() and 'PH_Ch{c2}'.format(c2=ch2) in self.ph_snr_ch_varz.keys():
						tempcuts = self.trans_grid.cuts_man.AndCuts([self.ph_snr_h_cuts[cells]['PH_H{c1}'.format(c1=ch+1)], self.ph_snr_ch_cuts[cells]['PH_Ch{c2}'.format(c2=ch2)]])
						# tempcuts = self.trans_grid.cuts_man.ConcatenateCuts(self.ph_snr_h_cuts[cells]['PH_H{c1}'.format(c1=ch+1)], self.ph_snr_ch_cuts[cells]['PH_Ch{c2}'.format(c2=ch2)])
						DrawHistogram('PH_H{c1}_Vs_PH_Ch{c2}_snr'.format(c1=ch+1, c2=ch2), 'PH cluster ch{c2} [SNR]'.format(c2=ch2), -55, 265, 'PH highest {c1}{sf} ch [SNR]'.format(c1=ch+1, sf='st' if ch == 0 else 'nd' if ch == 1 else 'rd' if ch == 2 else 'th'), self.ph_snr_ch_varz['PH_Ch{i}'.format(i=ch2)], self.ph_snr_h_varz['PH_H{i}'.format(i=ch+1)], tempcuts, 'snr')

	def DoClusterStudies(self, cells='all', typ='adc'):
		self.SetClusterStudies(cells)
		self.Do2DProfileMaps(cells)
		self.DoStrips2DHistograms(cells, typ=typ)
		self.DoPHStripCorrelations(cells, typ=typ)
		self.Do1DHistograms(cells, False, typ=typ)
		# self.Do1DHistograms(cells, True, typ='adc')

	def SetClusterStudies(self, cells='all'):
		self.minz[cells] = self.trans_grid.minz[cells]
		self.maxz[cells] = self.trans_grid.maxz[cells]
		self.GetCutsFromCutManager(cells)
		self.GetVarzFromTranspGrid()

	def GetCutsFromCutManager(self, cells):
		self.noise_cuts[cells] = self.trans_grid.cuts_man.noise_cuts[cells]
		self.ph_adc_ch_cuts[cells] = self.trans_grid.cuts_man.ph_adc_ch_cuts[cells]
		self.ph_snr_ch_cuts[cells] = self.trans_grid.cuts_man.ph_snr_ch_cuts[cells]
		self.ph_adc_h_cuts[cells] = self.trans_grid.cuts_man.ph_adc_h_cuts[cells]
		self.ph_snr_h_cuts[cells] = self.trans_grid.cuts_man.ph_snr_h_cuts[cells]
		self.phN_adc_ch_cuts[cells] = self.trans_grid.cuts_man.phN_adc_ch_cuts[cells]
		self.phN_snr_ch_cuts[cells] = self.trans_grid.cuts_man.phN_snr_ch_cuts[cells]
		self.phN_adc_h_cuts[cells] = self.trans_grid.cuts_man.phN_adc_h_cuts[cells]
		self.phN_snr_h_cuts[cells] = self.trans_grid.cuts_man.phN_snr_h_cuts[cells]

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


if __name__ == '__main__':
	c = ClusterChannelsAnalysis(None, 0, 0)
