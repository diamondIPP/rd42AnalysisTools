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
		self.phdelta = 0
		self.window_shift = 3
		# self.min_snr_neg, self.max_snr_neg, self.delta_snr = -64.25, 0.25, 0.125
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

		self.suffix = {'all': 'all', 'good': 'selection', 'bad': 'not_selection'}

		self.noise_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.noise_friend_cuts = {t: '' for t in ['all', 'good', 'bad']}

		self.in_transp_cluster = ''

		self.noise_varz = {}
		self.noise_friend_varz = {}

	def PosCanvas(self, canvas_name):
		self.w = PositionCanvas(self.trans_grid, canvas_name, self.w, self.window_shift)

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

	def PlotNoiseNotInCluster(self, cells='all', optending=''):
		temp_cut_noise = self.noise_cuts[cells]
		temph = ro.TH1F('temph0', 'temph0', int(RoundInt((self.max_adc_noise - self.min_adc_noise) / float(self.delta_adc_noise))), self.min_adc_noise, self.max_adc_noise)
		self.trans_grid.trans_tree.Draw('diaChSignal>>temph0', temp_cut_noise, 'goff')
		mean, sigma = temph.GetMean(), temph.GetRMS()
		temph.Delete()
		self.min_snr_noise, self.max_snr_noise, self.delta_snr_noise = (ni / float(sigma) for ni in [self.min_adc_noise, self.max_adc_noise, self.delta_adc_noise])
		suffix = self.suffix[cells]
		self.trans_grid.DrawHisto1D('signal_noise_{c}_snr{s}'.format(c=suffix, s=optending), self.min_snr_noise, self.max_snr_noise, self.delta_snr_noise, self.noise_varz['snr'], varname='Signal not in cluster (SNR)', cuts=temp_cut_noise, option='e hist')
		self.trans_grid.FitGaus('signal_noise_{c}_snr{s}'.format(c=suffix, s=optending))
		self.trans_grid.histo['signal_noise_{c}_snr{s}'.format(c=suffix, s=optending)].GetXaxis().SetRangeUser(-3.2, 3.2)
		self.PosCanvas('signal_noise_{c}_snr{s}'.format(c=suffix, s=optending))
		self.trans_grid.DrawHisto1D('signal_noise_{c}_adc{s}'.format(c=suffix, s=optending), self.min_adc_noise, self.max_adc_noise, self.delta_adc_noise, self.noise_varz['adc'], varname='Signal not in cluster (SNR)', cuts=temp_cut_noise, option='e hist')
		self.trans_grid.FitGaus('signal_noise_{c}_adc{s}'.format(c=suffix, s=optending))
		self.trans_grid.histo['signal_noise_{c}_adc{s}'.format(c=suffix, s=optending)].GetXaxis().SetRangeUser(-32, 32)
		self.PosCanvas('signal_noise_{c}_adc{s}'.format(c=suffix, s=optending))

	def PlotFriendNoiseNotInCluster(self, cells='all'):
		if self.trans_grid.trans_tree.GetFriend('pedTree'):
			optending = 'buffer_{v}'.format(v=int(RoundInt(self.trans_grid.trans_tree.GetMaximum('pedTree.slidingLength'))))
			temp_cut_noise = self.noise_friend_cuts[cells]
			temph = ro.TH1F('temph0', 'temph0', int(RoundInt((self.max_adc_noise - self.min_adc_noise) / float(self.delta_adc_noise))), self.min_adc_noise, self.max_adc_noise)
			self.trans_grid.trans_tree.Draw('pedTree.diaChSignal>>temph0', temp_cut_noise, 'goff')
			mean, sigma = temph.GetMean(), temph.GetRMS()
			temph.Delete()
			self.min_snr_noise, self.max_snr_noise, self.delta_snr_noise = (ni / float(sigma) for ni in [self.min_adc_noise, self.max_adc_noise, self.delta_adc_noise])
			suffix = self.suffix[cells]
			self.trans_grid.DrawHisto1D('signal_noise_{s}_{c}_snr'.format(c=suffix, s=optending), self.min_snr_noise, self.max_snr_noise, self.delta_snr_noise, self.noise_friend_varz['snr'], varname='Signal not in cluster (SNR)', cuts=temp_cut_noise, option='e hist')
			self.trans_grid.FitGaus('signal_noise_{s}_{c}_snr'.format(c=suffix, s=optending))
			self.trans_grid.histo['signal_noise_{s}_{c}_snr'.format(c=suffix, s=optending)].GetXaxis().SetRangeUser(-3.2, 3.2)
			self.PosCanvas('signal_noise_{s}_{c}_snr'.format(c=suffix, s=optending))
			self.trans_grid.DrawHisto1D('signal_noise_{s}_{c}_adc'.format(c=suffix, s=optending), self.min_adc_noise, self.max_adc_noise, self.delta_adc_noise, self.noise_friend_varz['adc'], varname='Signal not in cluster (SNR)', cuts=temp_cut_noise, option='e hist')
			self.trans_grid.FitGaus('signal_noise_{s}_{c}_adc'.format(c=suffix, s=optending))
			self.trans_grid.histo['signal_noise_{s}_{c}_adc'.format(c=suffix, s=optending)].GetXaxis().SetRangeUser(-32, 32)
			self.PosCanvas('signal_noise_{s}_{c}_adc'.format(c=suffix, s=optending))
		else:
			print 'The transparent tree has no pedTree friend. Cannot do these plots'

	def DoProfileMaps(self):
		minev, maxev = self.trans_grid.trans_tree.GetMinimum('event'), self.trans_grid.trans_tree.GetMaximum('event')

		self.trans_grid.DrawProfile1D('cm_event_profile', minev + self.delta_ev / 2.0, maxev - self.delta_ev / 2.0, self.delta_ev, 'event', 'event', 'cmn', 'cm [ADC]', self.trans_grid.cuts_man.transp_ev)
		self.PosCanvas('cm_event_profile')

		self.trans_grid.DrawProfile2D('adc_channel_event_profile', minev + self.delta_ev / 2.0, maxev - self.delta_ev / 2.0, self.delta_ev, 'event', 0, 127, 1, 'VA channel', 'event', 'diaChannels', 'diaChADC', 'ADC', self.trans_grid.cuts_man.transp_ev)
		self.PosCanvas('adc_channel_event_profile')
		self.trans_grid.DrawProfile2D('signal_channel_event_profile', minev + self.delta_ev / 2.0, maxev - self.delta_ev / 2.0, self.delta_ev, 'event', 0, 127, 1, 'VA channel', 'event', 'diaChannels', 'diaChSignal', 'Signal [ADC]', self.trans_grid.cuts_man.transp_ev)
		self.PosCanvas('signal_channel_event_profile')
		self.trans_grid.DrawProfile2D('noise_cmc_channel_event_profile_all', minev + self.delta_ev / 2.0, maxev - self.delta_ev / 2.0, self.delta_ev, 'event', 0, 127, 1, 'VA channel', 'event', 'diaChannels', 'diaChPedSigmaCmc', 'Noise [ADC]', self.trans_grid.cuts_man.transp_ev)
		self.PosCanvas('noise_cmc_channel_event_profile_all')

	def DoFriendProfileMaps(self):
		optending = 'buffer_{v}'.format(v=int(RoundInt(self.trans_grid.trans_tree.GetMaximum('pedTree.slidingLength'))))
		minev, maxev = self.trans_grid.trans_tree.GetMinimum('event'), self.trans_grid.trans_tree.GetMaximum('event')

		self.trans_grid.DrawProfile1D('cm_event_profile_{s}'.format(s=optending), minev + self.delta_ev / 2.0, maxev - self.delta_ev / 2.0, self.delta_ev, 'event', 'event', 'pedTree.cmn', 'cm [ADC]', self.trans_grid.cuts_man.transp_ev)
		self.PosCanvas('cm_event_profile_{s}'.format(s=optending))

		# self.trans_grid.DrawProfile2D('adc_channel_event_profile_{s}'.format(s=optending), minev + self.delta_ev / 2.0, maxev - self.delta_ev / 2.0, self.delta_ev, 'event', 0, 127, 1, 'VA channel', 'event', 'diaChannels', 'diaChADC', 'ADC', self.trans_grid.cuts_man.transp_ev)
		# self.PosCanvas('adc_channel_event_profile_{s}'.format(s=optending))
		self.trans_grid.DrawProfile2D('signal_channel_event_profile_{s}'.format(s=optending), minev + self.delta_ev / 2.0, maxev - self.delta_ev / 2.0, self.delta_ev, 'event', 0, 127, 1, 'VA channel', 'event', 'diaChannels', 'pedTree.diaChSignal', 'Signal [ADC]', self.trans_grid.cuts_man.transp_ev)
		self.PosCanvas('signal_channel_event_profile_{s}'.format(s=optending))
		self.trans_grid.DrawProfile2D('noise_cmc_channel_event_profile_all_{s}'.format(s=optending), minev + self.delta_ev / 2.0, maxev - self.delta_ev / 2.0, self.delta_ev, 'event', 0, 127, 1, 'VA channel', 'event', 'diaChannels', 'pedTree.diaChPedSigmaCmc', 'Noise [ADC]', self.trans_grid.cuts_man.transp_ev)
		self.PosCanvas('noise_cmc_channel_event_profile_all_{s}'.format(s=optending))

	def DoStrips2DHistograms(self):
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

		tempcuts = self.trans_grid.cuts_man.ConcatenateCuts(self.trans_grid.cuts_man.not_in_cluster, self.trans_grid.cuts_man.valid_ped_sigma)
		minz, maxz = -322.5, 322.5
		DrawHistogram('noise_Vs_channel_adc', minz, maxz, 'signal noise [ADC]', self.noise_varz['adc'], tempcuts, 'adc')
		minz, maxz = -32.25, 32.25
		DrawHistogram('noise_Vs_channel_snr', minz, maxz, 'signal noise [SNR]', self.noise_varz['snr'], tempcuts, 'snr')

	def DoFriendStrips2DHistograms(self):
		optending = 'buffer_{v}'.format(v=int(RoundInt(self.trans_grid.trans_tree.GetMaximum('pedTree.slidingLength'))))
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

		tempcuts = self.trans_grid.cuts_man.ConcatenateCuts(self.trans_grid.cuts_man.not_in_cluster, self.trans_grid.cuts_man.valid_ped_new_sigma)
		minz, maxz = -322.5, 322.5
		DrawHistogram('noise_Vs_channel_adc_{s}'.format(s=optending), minz, maxz, 'signal noise [ADC]', self.noise_friend_varz['adc'], tempcuts, 'adc')
		minz, maxz = -32.25, 32.25
		DrawHistogram('noise_Vs_channel_snr_{s}'.format(s=optending), minz, maxz, 'signal noise [SNR]', self.noise_friend_varz['snr'], tempcuts, 'snr')

	def DoPedestalEventHistograms(self, doFriend=False):
		tempcuts = self.in_transp_cluster
		minev, maxev = self.trans_grid.trans_tree.GetMinimum('event'), self.trans_grid.trans_tree.GetMaximum('event')
		deltaev = 100.
		if not doFriend:
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

	def DoNoiseAnalysis(self, cells='all'):
		self.GetCutsFromCutManager(cells)
		self.GetVarzFromTranspGrid()
		self.DoProfileMaps()
		self.DoPedestalEventHistograms(False)
		self.DoStrips2DHistograms()
		self.PlotNoiseNotInCluster(cells)

	def DoFriendNoiseAnalysis(self, cells='all'):
		if self.trans_grid.trans_tree.GetFriend('pedTree'):
			self.GetCutsFromCutManager(cells)
			self.GetVarzFromTranspGrid()
			self.DoFriendProfileMaps()
			self.DoPedestalEventHistograms(True)
			self.DoFriendStrips2DHistograms()
			self.PlotFriendNoiseNotInCluster(cells)
		else:
			print 'The transparent tree has no pedTree friend. Cannot do these plots'

	def GetCutsFromCutManager(self, cells):
		self.noise_cuts[cells] = self.trans_grid.cuts_man.noise_cuts[cells]
		self.noise_friend_cuts[cells] = self.trans_grid.cuts_man.noise_friend_cuts[cells]
		self.in_transp_cluster = self.trans_grid.cuts_man.ConcatenateCuts(cut2=self.trans_grid.cuts_man.in_transp_cluster, cut1=self.trans_grid.cuts_man.transp_ev)

	def GetVarzFromTranspGrid(self):
		self.noise_varz = self.trans_grid.noise_varz
		self.noise_friend_varz = self.trans_grid.noise_friend_varz

if __name__ == '__main__':
	c = NoiseAnalysis(None, 0, 0)
