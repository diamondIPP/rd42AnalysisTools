#!/usr/bin/env python
import ROOT as ro
import numpy as np
import ipdb, os
from ConfigParser import ConfigParser
from TransparentGrid import TransparentGrid
from optparse import OptionParser
from Utils import *

color_index = 10000

class NoiseAnalysis:
	def __init__(self, trans_grid, numstrips, clustersize):
		self.window_shift = 3
		self.delta_ev = 100
		self.min_adc_noise, self.max_adc_noise, self.delta_adc_noise = -322.25, 322.25, 0.5
		self.min_snr_noise, self.max_snr_noise, self.delta_snr_noise = -32.225, 32.225, 0.05
		self.trash = []
		self.w = 0
		self.trans_grid = trans_grid
		self.num_strips = numstrips
		self.cluster_size = clustersize

		self.suffix = GetSuffixDictionary(self.trans_grid)

		self.noise_cuts = self.trans_grid.cuts_man.noise_cuts
		self.noise_friend_cuts = self.trans_grid.cuts_man.noise_friend_cuts
		self.noise_nc_cuts = self.trans_grid.cuts_man.noise_nc_cuts
		self.noise_nc_friend_cuts = self.trans_grid.cuts_man.noise_nc_friend_cuts

		self.in_transp_cluster = self.trans_grid.cuts_man.AndCuts([self.trans_grid.cuts_man.transp_ev, self.trans_grid.cuts_man.in_transp_cluster])

		self.noise_varz = self.trans_grid.noise_varz
		self.noise_friend_varz = self.trans_grid.noise_friend_varz

	def PosCanvas(self, canvas_name):
		self.w = PositionCanvas(self.trans_grid, canvas_name, self.w, self.window_shift)

	def OverlayNoiseDistribution(self, histo, cells='all', isFriend=False):
		suffix = self.suffix[cells]
		hname = histo.GetName().split('h_')[1]
		typ = 'adc' if 'adc' in hname.lower() else 'snr'
		noise_name0 = 'signal_noise_{s}_{t}'.format(s=suffix, t=typ) if not isFriend else 'signal_noise_buffer_{b}_{s}_{t}'.format(s=suffix, t=typ, b=self.trans_grid.noise_friend_buffer)
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

	def PlotNoise1D(self, varzdic, name, cut, typ='adc'):
		temp_cut_noise = cut
		temph = ro.TH1F('temph0', 'temph0', int(RoundInt((self.max_adc_noise - self.min_adc_noise) / float(self.delta_adc_noise))), self.min_adc_noise, self.max_adc_noise)
		self.trans_grid.trans_tree.Draw('{v}>>temph0'.format(v=varzdic['adc']), temp_cut_noise, 'goff')
		mean, sigma = temph.GetMean(), temph.GetRMS()
		temph.Delete()
		if sigma > 0:
			self.min_snr_noise, self.max_snr_noise, self.delta_snr_noise = (ni / float(sigma) for ni in [self.min_adc_noise, self.max_adc_noise, self.delta_adc_noise])
			if typ == 'snr':
				self.trans_grid.DrawHisto1D(name + '_snr', self.min_snr_noise, self.max_snr_noise, self.delta_snr_noise, varzdic['snr'], varname='Signal not in cluster [SNR]', cuts=temp_cut_noise, option='e hist')
				self.trans_grid.FitGaus(name + '_snr')
				self.trans_grid.histo[name + '_snr'].GetXaxis().SetRangeUser(-3.2, 3.2)
				self.PosCanvas(name + '_snr')
			else:
				self.trans_grid.DrawHisto1D(name + '_adc', self.min_adc_noise, self.max_adc_noise, self.delta_adc_noise, varzdic['adc'], varname='Signal not in cluster [ADC]', cuts=temp_cut_noise, option='e hist')
				self.trans_grid.FitGaus(name + '_adc')
				self.trans_grid.histo[name + '_adc'].GetXaxis().SetRangeUser(-32, 32)
				self.PosCanvas(name + '_adc')

	def PlotNoiseNotInCluster(self, cells='all', typ='adc', doNC=False, isFriend=False):
		suffix = self.suffix[cells]
		if doNC:
			nameh = 'signal_noise_NC_chs_{c}'.format(c=suffix) if not isFriend else 'signal_noise_NC_chs_buffer_{b}_{c}'.format(c=suffix, b=self.trans_grid.noise_friend_buffer)
		else:
			nameh = 'signal_noise_{c}'.format(c=suffix) if not isFriend else 'signal_noise_buffer_{b}_{c}'.format(c=suffix, b=self.trans_grid.noise_friend_buffer)
		temp_cut_noise = self.noise_cuts[cells] if not isFriend and not doNC else self.noise_friend_cuts[cells] if isFriend and not doNC else self.noise_nc_cuts[cells] if not isFriend and doNC else self.noise_nc_friend_cuts[cells]
		var = self.noise_varz if not isFriend else self.noise_friend_varz
		self.PlotNoise1D(var, nameh, temp_cut_noise, typ=typ)

	def DoCommonMode(self, isFriend=False):
		minev, maxev = self.trans_grid.trans_tree.GetMinimum('event'), self.trans_grid.trans_tree.GetMaximum('event')
		var = 'cmn' if not isFriend else 'pedTree.cmn'
		hname = 'cm_event_profile' if not isFriend else 'cm_event_profile_buffer_{v}'.format(v=RoundInt(self.trans_grid.noise_friend_buffer))
		self.trans_grid.DrawProfile1D(hname, minev + self.delta_ev / 2.0, maxev - self.delta_ev / 2.0, self.delta_ev, 'event', 'event', var, 'cm [ADC]', self.trans_grid.cuts_man.transp_ev)
		self.PosCanvas(hname)

		hname = 'cm_histo' if not isFriend else 'cm_histo_buffer_{v}'.format(v=RoundInt(self.trans_grid.noise_friend_buffer))
		self.trans_grid.DrawHisto1D(hname, -31, 31, 2, var, 'cm [ADC]', self.trans_grid.cuts_man.transp_ev)
		self.PosCanvas(hname)
	# self.trans_grid.DrawProfile2D('adc_channel_event_profile', minev + self.delta_ev / 2.0, maxev - self.delta_ev / 2.0, self.delta_ev, 'event', 0, 127, 1, 'VA channel', 'event', 'diaChannels', 'diaChADC', 'ADC', self.trans_grid.cuts_man.transp_ev)
	# self.PosCanvas('adc_channel_event_profile')
	# self.trans_grid.DrawProfile2D('signal_channel_event_profile', minev + self.delta_ev / 2.0, maxev - self.delta_ev / 2.0, self.delta_ev, 'event', 0, 127, 1, 'VA channel', 'event', 'diaChannels', 'diaChSignal', 'Signal [ADC]', self.trans_grid.cuts_man.transp_ev)
	# self.PosCanvas('signal_channel_event_profile')
	# self.trans_grid.DrawProfile2D('noise_cmc_channel_event_profile_all', minev + self.delta_ev / 2.0, maxev - self.delta_ev / 2.0, self.delta_ev, 'event', 0, 127, 1, 'VA channel', 'event', 'diaChannels', 'diaChPedSigmaCmc', 'Noise [ADC]', self.trans_grid.cuts_man.transp_ev)
	# self.PosCanvas('noise_cmc_channel_event_profile_all')

	def DoStrips2DHistograms(self, typ='adc', isFriend=False):
		minch, maxch, deltach, xname, xvar = -0.5, 127.5, 1, 'VA channel', 'diaChannels'
		minch_plot, maxch_plot = int(self.trans_grid.ch_ini - np.ceil((self.cluster_size - 1) / 2.0)), int(self.trans_grid.ch_end + np.ceil((self.cluster_size - 1) / 2.0))

		def DrawHistogram(name, zmin, zmax, yname, yvar, cuts, typ='adc'):
			histo_limits = Get1DLimits(zmin, zmax, self.delta_adc_noise) if typ == 'adc' else Get1DLimits(zmin, zmax, self.delta_snr_noise)
			deltay = self.delta_adc_noise if typ == 'adc' else self.delta_snr_noise
			miny_plot, maxy_plot = (-50, 50) if typ == 'adc' else (-5, 5)
			self.trans_grid.DrawHisto2D(name, minch, maxch, deltach, xname, histo_limits['min'], histo_limits['max'], deltay, yname, xvar, yvar, cuts)
			self.trans_grid.histo[name].GetXaxis().SetRangeUser(minch_plot, maxch_plot)
			self.trans_grid.histo[name].GetYaxis().SetRangeUser(miny_plot, maxy_plot)
			self.PosCanvas(name)

		tempcuts = self.trans_grid.cuts_man.AndCuts([self.trans_grid.cuts_man.not_in_transp_cluster, self.trans_grid.cuts_man.valid_ped_sigma])
		# tempcuts = self.trans_grid.cuts_man.ConcatenateCuts(self.trans_grid.cuts_man.not_in_cluster, self.trans_grid.cuts_man.valid_ped_sigma)
		minz, maxz = (self.min_adc_noise, self.max_adc_noise) if typ == 'adc' else (self.min_adc_noise / 10., self.max_adc_noise / 10.)
		hname = 'noise_Vs_channel_{t}'.format(t=typ) if not isFriend else 'noise_buffer_{v}_Vs_channel_{t}'.format(t=typ, v=self.trans_grid.noise_friend_buffer)
		var = self.noise_varz['adc'] if typ == 'adc' and not isFriend else self.noise_friend_varz['adc'] if typ == 'adc' and isFriend else self.noise_friend_varz['snr'] if typ == 'snr' and isFriend else self.noise_varz['snr']
		DrawHistogram(hname, minz, maxz, 'signal noise [{t}]'.format(t=typ.upper()), var, tempcuts, typ)

	def DoPedestalEventHistograms(self, isFriend=False):
		tempcuts = self.in_transp_cluster
		minev, maxev = self.trans_grid.trans_tree.GetMinimum('event'), self.trans_grid.trans_tree.GetMaximum('event')
		deltaev = 100.
		if not isFriend:
			varz = 'diaChPedMeanCmc'
			maxy = GetMaximumFromTree(self.trans_grid.trans_tree, 'diaChPedMeanCmc', tempcuts)
			miny = GetMinimumFromTree(self.trans_grid.trans_tree, 'diaChPedMeanCmc', tempcuts)
			nameh = 'pedestal_mean_vs_event'
		else:
			if not self.trans_grid.trans_tree.GetFriend('pedTree'):
				print 'The transparent tree has no pedTree friend. Cannot do these plots'
				return
			optending = 'buffer_{v}'.format(v=int(RoundInt(self.trans_grid.trans_tree.GetMaximum('pedTree.slidingLength'))))
			varz = 'pedTree.diaChPedMeanCmc'
			maxy = GetMaximumFromTree(self.trans_grid.trans_tree, 'pedTree.diaChPedMeanCmc', tempcuts)
			miny = GetMinimumFromTree(self.trans_grid.trans_tree, 'pedTree.diaChPedMeanCmc', tempcuts)
			nameh = 'pedestal_mean_{s}_vs_event'.format(s=optending)

		self.trans_grid.DrawHisto2D(nameh, minev, maxev, deltaev, 'event', miny, maxy, 1.0, 'pedestal mean cluster chs [ADC]', 'event', varz, tempcuts)
		self.PosCanvas(nameh)

	def DoNoiseVsEventStudies(self, cells='all', num_delta_ev=100, doNC=False, isFriend=False):
		suffix = self.suffix[cells]
		xmin, xmax, deltax = self.trans_grid.trans_tree.GetMinimum('event'), self.trans_grid.trans_tree.GetMaximum('event'), num_delta_ev * self.delta_ev
		hlimitsx = Get1DLimits(xmin, xmax, deltax, oddbins=False)
		deltay = self.delta_adc_noise * 4
		hlimitsy = Get1DLimits(RoundInt(self.min_adc_noise / 10.0), RoundInt(self.max_adc_noise / 10.0), deltay)
		if not isFriend:
			nameh = 'signal_noise_Vs_event_{c}'.format(c=suffix) if not doNC else 'signal_noise_NC_chs_Vs_event_{c}'.format(c=suffix)
			temp_cut_noise = self.noise_cuts[cells] if not doNC else self.noise_nc_cuts[cells]
			varz = 'diaChSignal'
		else:
			if not self.trans_grid.trans_tree.GetFriend('pedTree'):
				print 'The transparent tree has no pedTree friend. Cannot do these plots'
				return
			optending = 'buffer_{v}'.format(v=int(RoundInt(self.trans_grid.trans_tree.GetMaximum('pedTree.slidingLength'))))
			nameh = 'signal_noise_{o}_Vs_event_{c}'.format(o=optending, c=suffix) if not doNC else 'signal_noise_NC_chs_{o}_Vs_event_{c}'.format(o=optending, c=suffix)
			temp_cut_noise = self.noise_friend_cuts[cells] if not doNC else self.noise_nc_friend_cuts[cells]
			varz = 'pedTree.diaChSignal'

		self.trans_grid.DrawHisto2D(nameh, hlimitsx['min'], hlimitsx['max'], deltax, 'event', hlimitsy['min'], hlimitsy['max'], deltay, 'signal noise [ADC]', 'event', varz, temp_cut_noise)
		self.PosCanvas(nameh)

		tempArray = ro.TObjArray()
		self.trans_grid.histo[nameh].FitSlicesY(0, 0, -1, 0, 'QNR', tempArray)
		if not doNC:
			nameh2 = 'signal_noise_fitted_sigma_Vs_event_{c}'.format(c=suffix) if not isFriend else 'signal_noise_buffer_{o}_Fitted_sigma_Vs_event_{c}'.format(o=self.trans_grid.noise_friend_buffer, c=suffix)
		else:
			nameh2 = 'signal_noise_NC_chs_fitted_sigma_Vs_event_{c}'.format(c=suffix) if not isFriend else 'signal_noise_NC_chs_buffer_{o}_Fitted_sigma_Vs_event_{c}'.format(o=self.trans_grid.noise_friend_buffer, c=suffix)
		self.trans_grid.histo[nameh2] = tempArray[2]
		self.trans_grid.histo[nameh2].SetNameTitle('h_' + nameh2, 'h_' + nameh2)
		self.trans_grid.canvas[nameh2] = ro.TCanvas('c_' + nameh2, 'c_' + nameh2, 1)
		self.trans_grid.histo[nameh2].Draw('e hist')
		ro.gPad.Update()
		SetDefault1DCanvasSettings(self.trans_grid.canvas[nameh2])
		SetDefault1DStats(self.trans_grid.histo[nameh2], y1=0.15, y2=0.45, optstat=1000000001)
		self.trans_grid.FitPol(nameh2, 1)
		self.PosCanvas(nameh2)

	def DoNoiseAnalysis(self, cells='all', typ='adc', isFriend=False):
		if isFriend:
			if not self.trans_grid.trans_tree.GetFriend('pedTree'):
				print 'The transparent tree has no pedTree friend. Add it first in transparent grid class'
				return
		self.DoCommonMode(isFriend)
		self.DoPedestalEventHistograms(isFriend)
		self.DoStrips2DHistograms(typ, isFriend)
		self.PlotNoiseNotInCluster(cells, typ, False, isFriend)
		self.DoNoiseVsEventStudies(cells, isFriend=isFriend)
		self.PlotNoiseNotInCluster('all', typ, True, isFriend)
		self.DoNoiseVsEventStudies('all', doNC=True, isFriend=isFriend)

if __name__ == '__main__':
	c = NoiseAnalysis(None, 0, 0)
