#!/usr/bin/env python
import numpy as np
import ipdb
from Utils import *

class CutManager:
	def __init__(self, tree, numstrips=5, clust_size=10, sat_adc=4095):
		self.tree = tree
		self.numstrips = numstrips
		self.cluster_size = clust_size
		self.sat_adc = sat_adc
		self.transp_ev = '(transparentEvent)'
		self.tcutg_cell = {}
		self.selected_cells = '{s}'
		self.not_selected_cells = '{ns}'
		self.all_cells = '({s}||{ns})'
		self.selected_cells_centers = '{s}'
		self.not_selected_cells_centers = '{ns}'
		self.no_up_down_borders = '(({l}<=diaChYPred)&&(diaChYPred<={h}))'

		self.good_chs, self.bad_chs, self.only_bad_chs = [], [], []
		self.good_chs_cut, self.only_bad_chs_cut = '', ''

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

		self.in_transp_cluster = '((diaChSeed)||(diaChHits))'
		self.not_in_cluster = '((!diaChSeed)&&(!diaChHits))'
		self.not_in_transp_cluster = '((!diaChSeed)&&(!diaChHits)&&(!diaChsNoisy)&&(!diaChsScreened)&&(!diaChsNC))'
		# self.not_in_transp_cluster = '((!diaChSeed)&&(!diaChHits)&&(!diaChsNoisy)&&(!diaChsNC))'
		# self.not_in_transp_cluster = '((!diaChSeed)&&(!diaChHits)&&(!diaChsNoisy)&&(!diaChsScreened))'
		# self.not_in_transp_cluster = '((!diaChSeed)&&(!diaChHits)&&(!diaChsNoisy))'
		self.any_saturated = '(diaChADC=={s})'.format(s=self.sat_adc)
		self.non_saturated = '(diaChADC!={s})'.format(s=self.sat_adc)
		self.valid_ped_sigma = '(diaChPedSigmaCmc>0)'
		self.valid_ped_friend_sigma = '(pedTree.diaChPedSigmaCmc>0)'
		self.valid_ped_sigma_ch = {}
		self.valid_ped_sigma_h = {}
		self.valid_ped_sigma_N_ch = {}
		self.valid_ped_sigma_N_h = {}

		self.sat_adc_ch = {}
		self.sat_adc_h = {}
		self.not_sat_adc_ch = {}
		self.not_sat_adc_h = {}

		self.sat_adc_N_ch = {}
		self.sat_adc_N_h = {}
		self.not_sat_adc_N_ch = {}
		self.not_sat_adc_N_h = {}
		
		self.ph_adc_ch = {}
		self.ph_snr_ch = {}
		self.ph_adc_h = {}
		self.ph_snr_h = {}
		
		self.ph_adc_N_ch = {}
		self.ph_snr_N_ch = {}
		self.ph_adc_N_h = {}
		self.ph_snr_N_h = {}

		self.sat_evts = None
		self.sat_evts_region = '(satRegion)'
		self.not_sat_evts_region = '(!satRegion)'
		self.nc_chs = '((!diaChSeed)&&(!diaChHits)&&(!diaChsNoisy)&&(!diaChsScreened)&&(diaChsNC))'
		
		self.noise_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.noise_friend_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.noise_nc_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.noise_nc_friend_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.ph_adc_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.ph_snr_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.ph_adc_h_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.ph_snr_h_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_adc_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_adc_h_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_snr_ch_cuts = {t: {} for t in ['all', 'good', 'bad']}
		self.phN_snr_h_cuts = {t: {} for t in ['all', 'good', 'bad']}

		self.in_central_rect_region = {}
		self.out_central_rect_region = {}

	def SetCuts(self, neg_cut_snr, neg_cut_adc):
		for i in xrange(self.cluster_size):
			if FindLeafInTree(self.tree, 'clusterChannel{i}'.format(i=i)):
				self.neg_snr_ph_ch['PH_Ch{i}'.format(i=i)] = '((diaChPedSigmaCmc[clusterChannel{i}]>0)&&(diaChSignal[clusterChannel{i}]/diaChPedSigmaCmc[clusterChannel{i}]<-{th}))'.format(i=i, th=neg_cut_snr)
				self.not_neg_snr_ph_ch['PH_Ch{i}'.format(i=i)] = '((diaChPedSigmaCmc[clusterChannel{i}]>0)&&(diaChSignal[clusterChannel{i}]/diaChPedSigmaCmc[clusterChannel{i}]>=-{th}))'.format(i=i, th=neg_cut_snr)
				self.neg_adc_ph_ch['PH_Ch{i}'.format(i=i)] = '(diaChSignal[clusterChannel{i}]<-{th})'.format(i=i, th=neg_cut_adc)
				self.not_neg_adc_ph_ch['PH_Ch{i}'.format(i=i)] = '(diaChSignal[clusterChannel{i}]>=-{th})'.format(i=i, th=neg_cut_adc)

				self.sat_adc_ch['Ch{i}'.format(i=i)] = '({sr}&&(diaChADC[clusterChannel{i}]>={s}))'.format(i=i, s=self.sat_adc, sr=self.sat_evts_region)
				self.not_sat_adc_ch['Ch{i}'.format(i=i)] = '({nsr}&&(diaChADC[clusterChannel{i}]<{s}))'.format(i=i, s=self.sat_adc, nsr=self.not_sat_evts_region)

				self.valid_ped_sigma_ch['Ch{i}'.format(i=i)] = '(diaChPedSigmaCmc[clusterChannel{i}]>0)'.format(i=i)

				self.ph_adc_ch['PH_Ch{i}'.format(i=i)] = self.transp_ev
				self.ph_snr_ch['PH_Ch{i}'.format(i=i)] = self.ConcatenateCuts(self.transp_ev, self.valid_ped_sigma_ch['Ch{i}'.format(i=i)])

			if FindLeafInTree(self.tree, 'clusterChannelHighest{i}'.format(i=i+1)):
				self.neg_snr_ph_h['PH_H{i}'.format(i=i+1)] = '((diaChPedSigmaCmc[clusterChannelHighest{i}]>0)&&(diaChSignal[clusterChannelHighest{i}]/diaChPedSigmaCmc[clusterChannelHighest{i}]<-{th}))'.format(i=i+1, th=neg_cut_snr)
				self.not_neg_snr_ph_h['PH_H{i}'.format(i=i+1)] = '((diaChPedSigmaCmc[clusterChannelHighest{i}]>0)&&(diaChSignal[clusterChannelHighest{i}]/diaChPedSigmaCmc[clusterChannelHighest{i}]>=-{th}))'.format(i=i+1, th=neg_cut_snr)
				self.neg_adc_ph_h['PH_H{i}'.format(i=i+1)] = '(diaChSignal[clusterChannelHighest{i}]<-{th})'.format(i=i+1, th=neg_cut_adc)
				self.not_neg_adc_ph_h['PH_H{i}'.format(i=i+1)] = '(diaChSignal[clusterChannelHighest{i}]>=-{th})'.format(i=i + 1, th=neg_cut_adc)

				self.sat_adc_h['H{i}'.format(i=i+1)] = '({sr}&&(diaChADC[clusterChannelHighest{i}]>={s}))'.format(i=i+1, s=self.sat_adc, sr=self.sat_evts_region)
				self.not_sat_adc_h['H{i}'.format(i=i+1)] = '({nsr}&&(diaChADC[clusterChannelHighest{i}]<{s}))'.format(i=i+1, s=self.sat_adc, nsr=self.not_sat_evts_region)

				self.valid_ped_sigma_h['H{i}'.format(i=i+1)] = '(diaChPedSigmaCmc[clusterChannelHighest{i}]>0)'.format(i=i+1)

				self.ph_adc_h['PH_H{i}'.format(i=i+1)] = self.transp_ev
				self.ph_snr_h['PH_H{i}'.format(i=i+1)] = self.ConcatenateCuts(self.transp_ev, self.valid_ped_sigma_h['H{i}'.format(i=i + 1)])

		for ch in xrange(self.cluster_size):
			list_neg_snr_phN_ch = []
			list_not_neg_snr_phN_ch = []
			list_neg_adc_phN_ch = []
			list_not_neg_adc_phN_ch = []
			list_neg_snr_phN_h = []
			list_not_neg_snr_phN_h = []
			list_neg_adc_phN_h = []
			list_not_neg_adc_phN_h = []

			list_sat_N_ch = []
			list_not_sat_N_ch = []
			list_sat_N_h = []
			list_not_sat_N_h = []

			list_valid_ped_sigma_N_ch = []
			list_valid_ped_sigma_N_h = []
			
			list_ph_adc_N_ch = []
			list_ph_snr_N_ch = []
			list_ph_adc_N_h = []
			list_ph_snr_N_h = []

			for chi in xrange(ch + 1):
				if 'PH_Ch{i}'.format(i=chi) in self.neg_snr_ph_ch.keys():
					list_neg_snr_phN_ch.append(self.neg_snr_ph_ch['PH_Ch{i}'.format(i=chi)])
					list_neg_adc_phN_ch.append(self.neg_adc_ph_ch['PH_Ch{i}'.format(i=chi)])
				if 'PH_Ch{i}'.format(i=chi) in self.not_neg_snr_ph_ch.keys():
					list_not_neg_snr_phN_ch.append(self.not_neg_snr_ph_ch['PH_Ch{i}'.format(i=chi)])
					list_not_neg_adc_phN_ch.append(self.not_neg_adc_ph_ch['PH_Ch{i}'.format(i=chi)])
				if 'PH_H{i}'.format(i=chi+1) in self.neg_snr_ph_h.keys():
					list_neg_snr_phN_h.append(self.neg_snr_ph_h['PH_H{i}'.format(i=chi+1)])
					list_neg_adc_phN_h.append(self.neg_adc_ph_h['PH_H{i}'.format(i=chi+1)])
				if 'PH_H{i}'.format(i=chi+1) in self.not_neg_snr_ph_h.keys():
					list_not_neg_snr_phN_h.append(self.not_neg_snr_ph_h['PH_H{i}'.format(i=chi+1)])
					list_not_neg_adc_phN_h.append(self.not_neg_adc_ph_h['PH_H{i}'.format(i=chi+1)])

				if 'Ch{i}'.format(i=chi) in self.sat_adc_ch.keys():
					list_sat_N_ch.append(self.sat_adc_ch['Ch{i}'.format(i=chi)])
					list_not_sat_N_ch.append(self.not_sat_adc_ch['Ch{i}'.format(i=chi)])
				if 'H{i}'.format(i=chi+1) in self.sat_adc_h.keys():
					list_sat_N_h.append(self.sat_adc_h['H{i}'.format(i=chi+1)])
					list_not_sat_N_h.append(self.not_sat_adc_h['H{i}'.format(i=chi+1)])

				if 'Ch{i}'.format(i=chi) in self.valid_ped_sigma_ch.keys():
					list_valid_ped_sigma_N_ch.append(self.valid_ped_sigma_ch['Ch{i}'.format(i=chi)])
				if 'H{i}'.format(i=chi+1) in self.valid_ped_sigma_h.keys():
					list_valid_ped_sigma_N_h.append(self.valid_ped_sigma_h['H{i}'.format(i=chi+1)])

				if 'PH_Ch{i}'.format(i=chi) in self.ph_adc_ch.keys():
					list_ph_adc_N_ch.append(self.ph_adc_ch['PH_Ch{i}'.format(i=chi)])
					list_ph_snr_N_ch.append(self.ph_snr_ch['PH_Ch{i}'.format(i=chi)])

				if 'PH_H{i}'.format(i=chi+1) in self.ph_adc_h.keys():
					list_ph_adc_N_h.append(self.ph_adc_h['PH_H{i}'.format(i=chi+1)])
					list_ph_snr_N_h.append(self.ph_snr_h['PH_H{i}'.format(i=chi+1)])


			self.neg_snr_phN_ch['PH{i}_Ch'.format(i=ch + 1)] = '(' + '||'.join(list_neg_snr_phN_ch) + ')' if len(list_neg_snr_phN_ch) > 0 else ''
			self.neg_adc_phN_ch['PH{i}_Ch'.format(i=ch + 1)] = '(' + '||'.join(list_neg_adc_phN_ch) + ')' if len(list_neg_adc_phN_ch) > 0 else ''
			self.not_neg_snr_phN_ch['PH{i}_Ch'.format(i=ch + 1)] = '(' + '&&'.join(list_not_neg_snr_phN_ch) + ')' if len(list_not_neg_snr_phN_ch) > 0 else ''
			self.not_neg_adc_phN_ch['PH{i}_Ch'.format(i=ch + 1)] = '(' + '&&'.join(list_not_neg_adc_phN_ch) + ')' if len(list_not_neg_adc_phN_ch) > 0 else ''
			self.neg_snr_phN_h['PH{i}_H'.format(i=ch + 1)] = '(' + '||'.join(list_neg_snr_phN_h) + ')' if len(list_neg_snr_phN_h) > 0 else ''
			self.neg_adc_phN_h['PH{i}_H'.format(i=ch + 1)] = '(' + '||'.join(list_neg_adc_phN_h) + ')' if len(list_neg_adc_phN_h) > 0 else ''
			self.not_neg_snr_phN_h['PH{i}_H'.format(i=ch + 1)] = '(' + '&&'.join(list_not_neg_snr_phN_h) + ')' if len(list_not_neg_snr_phN_h) > 0 else ''
			self.not_neg_adc_phN_h['PH{i}_H'.format(i=ch + 1)] = '(' + '&&'.join(list_not_neg_adc_phN_h) + ')' if len(list_not_neg_adc_phN_h) > 0 else ''

			self.sat_adc_N_ch['{i}_Ch'.format(i=ch + 1)] = '(' + '||'.join(list_sat_N_ch) + ')' if len(list_sat_N_ch) > 0 else ''
			self.not_sat_adc_N_ch['{i}_Ch'.format(i=ch + 1)] = '(' + '&&'.join(list_not_sat_N_ch) + ')' if len(list_not_sat_N_ch) > 0 else ''
			self.sat_adc_N_h['{i}_H'.format(i=ch + 1)] = '(' + '||'.join(list_sat_N_h) + ')' if len(list_sat_N_h) > 0 else ''
			self.not_sat_adc_N_h['{i}_H'.format(i=ch + 1)] = '(' + '&&'.join(list_not_sat_N_h) + ')' if len(list_not_sat_N_h) > 0 else ''

			self.valid_ped_sigma_N_ch['{i}_Ch'.format(i=ch + 1)] = '(' + '&&'.join(list_valid_ped_sigma_N_ch) + ')' if len(list_valid_ped_sigma_N_ch) > 0 else ''
			self.valid_ped_sigma_N_h['{i}_H'.format(i=ch + 1)] = '(' + '&&'.join(list_valid_ped_sigma_N_h) + ')' if len(list_valid_ped_sigma_N_h) > 0 else ''
			
			self.ph_adc_N_ch['PH{i}_Ch'.format(i=ch+1)] = '(' + '&&'.join(list_ph_adc_N_ch) + ')' if len(list_ph_adc_N_ch) > 0 else ''
			self.ph_snr_N_ch['PH{i}_Ch'.format(i=ch+1)] = '(' + '&&'.join(list_ph_snr_N_ch) + ')' if len(list_ph_snr_N_ch) > 0 else ''

			self.ph_adc_N_h['PH{i}_H'.format(i=ch+1)] = '(' + '&&'.join(list_ph_adc_N_h) + ')' if len(list_ph_adc_N_h) > 0 else ''
			self.ph_snr_N_h['PH{i}_H'.format(i=ch+1)] = '(' + '&&'.join(list_ph_snr_N_h) + ')' if len(list_ph_snr_N_h) > 0 else ''

	def SetChs(self, goodchs, badchs):
		self.good_chs = goodchs
		self.bad_chs = badchs
		self.only_bad_chs = [bi for bi in self.bad_chs if bi not in self.good_chs]

		self.good_chs_cut = '(' + '||'.join(['(diaChannels=={c})'.format(c=ch) for ch in self.good_chs]) + ')'
		self.only_bad_chs_cut = '(!' + self.good_chs_cut + ')'

	def SetCells(self, selection, not_selection):
		self.selected_cells = self.selected_cells.format(s=selection)
		self.not_selected_cells = self.not_selected_cells.format(ns=not_selection)
		self.all_cells = self.all_cells.format(s=selection, ns=not_selection)

	def SetCellsCenters(self, selection, not_selection):
		self.selected_cells_centers = self.selected_cells_centers.format(s=selection)
		self.not_selected_cells_centers = self.not_selected_cells_centers.format(ns=not_selection)

	def SetNoiseCuts(self):
		# self.noise_cuts = {cells: self.ConcatenateCutWithCells(cut=self.ConcatenateCuts(cut1=self.not_in_transp_cluster, cut2=self.valid_ped_sigma), cells=cells) for cells in ['all', 'good', 'bad']}
		self.noise_cuts = {'all': self.ConcatenateCuts(self.valid_ped_sigma, self.not_in_transp_cluster), 'good': self.ConcatenateCuts(self.valid_ped_sigma, self.ConcatenateCuts(self.not_in_transp_cluster, self.good_chs_cut)), 'bad': self.ConcatenateCuts(self.valid_ped_sigma, self.ConcatenateCuts(self.not_in_transp_cluster, self.only_bad_chs_cut))}
		# self.noise_friend_cuts = {cells: self.ConcatenateCutWithCells(cut=self.ConcatenateCuts(cut1=self.not_in_transp_cluster, cut2=self.valid_ped_friend_sigma), cells=cells) for cells in ['all', 'good', 'bad']}
		self.noise_friend_cuts = {'all': self.ConcatenateCuts(self.valid_ped_friend_sigma, self.not_in_transp_cluster), 'good': self.ConcatenateCuts(self.valid_ped_friend_sigma, self.ConcatenateCuts(self.not_in_transp_cluster, self.good_chs_cut)), 'bad': self.ConcatenateCuts(self.valid_ped_friend_sigma, self.ConcatenateCuts(self.not_in_transp_cluster, self.only_bad_chs_cut))}
		self.noise_nc_cuts = {cells: self.ConcatenateCutWithCells(cut=self.ConcatenateCuts(cut1=self.nc_chs, cut2=self.valid_ped_sigma), cells=cells) for cells in ['all', 'good', 'bad']}
		self.noise_nc_friend_cuts = {cells: self.ConcatenateCutWithCells(cut=self.ConcatenateCuts(cut1=self.nc_chs, cut2=self.valid_ped_friend_sigma), cells=cells) for cells in ['all', 'good', 'bad']}

	def SetPHCuts(self):
		for cells in ['all', 'good', 'bad']:
			for ch in xrange(self.cluster_size):
				if FindLeafInTree(self.tree, 'clusterChannel{c}'.format(c=ch)):
					self.ph_adc_ch_cuts[cells]['PH_Ch{c}'.format(c=ch)] = self.ConcatenateCutWithCells(cut=self.ph_adc_ch['PH_Ch{c}'.format(c=ch)], cells=cells)
					self.ph_snr_ch_cuts[cells]['PH_Ch{c}'.format(c=ch)] = self.ConcatenateCutWithCells(cut=self.ph_snr_ch['PH_Ch{c}'.format(c=ch)], cells=cells)
				if FindLeafInTree(self.tree, 'clusterChannelHighest{c}'.format(c=ch+1)):
					self.ph_adc_h_cuts[cells]['PH_H{c}'.format(c=ch+1)] = self.ConcatenateCutWithCells(cut=self.ph_adc_h['PH_H{c}'.format(c=ch+1)], cells=cells)
					self.ph_snr_h_cuts[cells]['PH_H{c}'.format(c=ch+1)] = self.ConcatenateCutWithCells(cut=self.ph_snr_h['PH_H{c}'.format(c=ch+1)], cells=cells)
				self.phN_adc_ch_cuts[cells]['PH{c}_Ch'.format(c=ch+1)] = self.ConcatenateCutWithCells(cut=self.ph_adc_N_ch['PH{c}_Ch'.format(c=ch+1)], cells=cells)
				self.phN_snr_ch_cuts[cells]['PH{c}_Ch'.format(c=ch+1)] = self.ConcatenateCutWithCells(cut=self.ph_snr_N_ch['PH{c}_Ch'.format(c=ch+1)], cells=cells)
				self.phN_adc_h_cuts[cells]['PH{c}_H'.format(c=ch+1)] = self.ConcatenateCutWithCells(cut=self.ph_adc_N_h['PH{c}_H'.format(c=ch+1)], cells=cells)
				self.phN_snr_h_cuts[cells]['PH{c}_H'.format(c=ch+1)] = self.ConcatenateCutWithCells(cut=self.ph_snr_N_h['PH{c}_H'.format(c=ch+1)], cells=cells)

	def SetUpDownBorderCuts(self, lower, upper):
		self.no_up_down_borders = self.no_up_down_borders.format(l=lower, h=upper)

	def GetThCut(self, var, th, cells, cuts='', op='>='):
		temp_cuts = '({v}{o}{th})'.format(v=var, o=op, th=th)
		temp_cuts = temp_cuts if cuts == '' else self.ConcatenateCuts(cuts, temp_cuts, '&&')
		return self.ConcatenateCutWithCells(temp_cuts, cells, '&&')

	def FindSaturationEvents(self):
		tempsat = self.ConcatenateCuts(self.transp_ev, self.any_saturated)
		leng = self.tree.Draw('event', tempsat, 'goff')
		if leng == -1:
			print 'Error, could not get the branch event. Check tree structure.'
		while leng > self.tree.GetEstimate():
			self.tree.SetEstimate(leng)
			leng = self.tree.Draw('event', tempsat, 'goff')
		buffev = self.tree.GetVal(0)
		self.sat_evts = [int(buffev[ev]) for ev in xrange(leng)]

	def IsEventInSaturationRegion(self, ev, skipAfter=0, skipBefore=0):
		for satev in self.sat_evts:
			if satev - skipBefore <= ev < satev + skipAfter:
				return True
		return False

	def ConcatenateCuts(self, cut1, cut2, operator='&&'):
		return '(' + operator.join([cut1, cut2]) + ')'

	def ConcatenateCutWithCells(self, cut, cells='all', operator='&&'):
		return self.ConcatenateCuts(self.selected_cells, cut, operator) if cells == 'good' else self.ConcatenateCuts(self.not_selected_cells, cut, operator) if cells == 'bad' else self.ConcatenateCuts(self.all_cells, cut, operator)


if __name__ == '__main__':
	cm = CutManager(None)
