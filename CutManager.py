#!/usr/bin/env python
import numpy as np
import ipdb
from Utils import *

class CutManager:
	def __init__(self, tree, numstrips=5, clust_size=10, sat_adc=4095):
		self.tree = tree
		self.numstrips = numstrips
		self.clust_size = clust_size
		self.sat_adc = sat_adc
		self.transp_ev = '(transparentEvent)'
		self.tcutg_cell = {}
		self.selected_cells = '{s}'
		self.not_selected_cells = '{ns}'
		self.all_cells = '({s}||{ns})'
		self.selected_cells_centers = '{s}'
		self.not_selected_cells_centers = '{ns}'
		self.no_up_down_borders = '(({l}<=diaChYPred)&&(diaChYPred<={h}))'

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

		self.not_in_transp_cluster = '((!diaChSeed)&&(!diaChHits)&&(!diaChsNoisy)&&(!diaChsScreened)&&(!diaChsNC))'
		self.any_saturated = '(diaChADC=={s})'.format(s=self.sat_adc)
		self.non_saturated = '(diaChADC!={s})'.format(s=self.sat_adc)
		self.valid_ped_sigma = '(diaChPedSigmaCmc>0)'
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

		self.sat_evts = None
		self.sat_evts_region = '(satRegion)'
		self.not_sat_evts_region = '(!satRegion)'

	def SetNegAndSatAndPedSigmaCuts(self, neg_cut_snr, neg_cut_adc):
		for i in xrange(self.clust_size):
			if FindLeafInTree(self.tree, 'clusterChannel{i}'.format(i=i)):
				self.neg_snr_ph_ch['PH_Ch{i}'.format(i=i)] = '((diaChPedSigmaCmc[clusterChannel{i}]>0)&&(diaChSignal[clusterChannel{i}]/diaChPedSigmaCmc[clusterChannel{i}]<-{th}))'.format(i=i, th=neg_cut_snr)
				self.not_neg_snr_ph_ch['PH_Ch{i}'.format(i=i)] = '((diaChPedSigmaCmc[clusterChannel{i}]>0)&&(diaChSignal[clusterChannel{i}]/diaChPedSigmaCmc[clusterChannel{i}]>=-{th}))'.format(i=i, th=neg_cut_snr)
				self.neg_adc_ph_ch['PH_Ch{i}'.format(i=i)] = '(diaChSignal[clusterChannel{i}]<-{th})'.format(i=i, th=neg_cut_adc)
				self.not_neg_adc_ph_ch['PH_Ch{i}'.format(i=i)] = '(diaChSignal[clusterChannel{i}]>=-{th})'.format(i=i, th=neg_cut_adc)

				self.sat_adc_ch['Ch{i}'.format(i=i)] = '(diaChADC[clusterChannel{i}]>={s})'.format(i=i, s=self.sat_adc)
				self.not_sat_adc_ch['Ch{i}'.format(i=i)] = '(diaChADC[clusterChannel{i}]<{s})'.format(i=i, s=self.sat_adc)

				self.valid_ped_sigma_ch['Ch{i}'.format(i=i)] = '(diaChPedSigmaCmc[clusterChannel{i}]>0)'.format(i=i)
			if FindLeafInTree(self.tree, 'clusterChannelHighest{i}'.format(i=i+1)):
				self.neg_snr_ph_h['PH_H{i}'.format(i=i+1)] = '((diaChPedSigmaCmc[clusterChannelHighest{i}]>0)&&(diaChSignal[clusterChannelHighest{i}]/diaChPedSigmaCmc[clusterChannelHighest{i}]<-{th}))'.format(i=i+1, th=neg_cut_snr)
				self.not_neg_snr_ph_h['PH_H{i}'.format(i=i+1)] = '((diaChPedSigmaCmc[clusterChannelHighest{i}]>0)&&(diaChSignal[clusterChannelHighest{i}]/diaChPedSigmaCmc[clusterChannelHighest{i}]>=-{th}))'.format(i=i+1, th=neg_cut_snr)
				self.neg_adc_ph_h['PH_H{i}'.format(i=i+1)] = '(diaChSignal[clusterChannelHighest{i}]<-{th})'.format(i=i+1, th=neg_cut_adc)
				self.not_neg_adc_ph_h['PH_H{i}'.format(i=i+1)] = '(diaChSignal[clusterChannelHighest{i}]>=-{th})'.format(i=i+1, th=neg_cut_adc)

				self.sat_adc_h['H{i}'.format(i=i+1)] = '(diaChADC[clusterChannelHighest{i}]>={s})'.format(i=i+1, s=self.sat_adc)
				self.not_sat_adc_h['H{i}'.format(i=i+1)] = '(diaChADC[clusterChannelHighest{i}]<{s})'.format(i=i+1, s=self.sat_adc)

				self.valid_ped_sigma_h['H{i}'.format(i=i+1)] = '(diaChPedSigmaCmc[clusterChannelHighest{i}]>0)'.format(i=i+1)

		for ch in xrange(self.numstrips):
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

			self.neg_snr_phN_ch['PH{i}_Ch'.format(i=ch + 1)] = '(' + '||'.join(list_neg_snr_phN_ch) + ')' if len(list_neg_snr_phN_ch) > 0 else ''
			self.neg_adc_phN_ch['PH{i}_Ch'.format(i=ch + 1)] = '(' + '||'.join(list_neg_adc_phN_ch) + ')' if len(list_neg_adc_phN_ch) > 0 else ''
			self.not_neg_snr_phN_ch['PH{i}_Ch'.format(i=ch + 1)] = '(' + '&&'.join(list_not_neg_snr_phN_ch) + ')' if len(list_not_neg_snr_phN_ch) > 0 else ''
			self.not_neg_adc_phN_ch['PH{i}_Ch'.format(i=ch + 1)] = '(' + '&&'.join(list_not_neg_adc_phN_ch) + ')' if len(list_not_neg_adc_phN_ch) > 0 else ''
			self.neg_snr_phN_h['PH{i}_H'.format(i=ch + 1)] = '(' + '||'.join(list_neg_snr_phN_h) + ')' if len(list_neg_snr_phN_h) > 0 else ''
			self.neg_adc_phN_h['PH{i}_H'.format(i=ch + 1)] = '(' + '||'.join(list_neg_adc_phN_h) + ')' if len(list_neg_adc_phN_h) > 0 else ''
			self.not_neg_snr_phN_h['PH{i}_H'.format(i=ch + 1)] = '(' + '&&'.join(list_not_neg_snr_phN_h) + ')' if len(list_not_neg_snr_phN_h) > 0 else ''
			self.not_neg_adc_phN_h['PH{i}_H'.format(i=ch + 1)] = '(' + '&&'.join(list_not_neg_adc_phN_h) + ')' if len(list_not_neg_adc_phN_h) > 0 else ''

			self.sat_adc_N_ch['{i}_Ch'.format(i=ch + 1)] = '(' + '||'.join(list_sat_N_ch) + ')' if len(list_sat_N_ch) > 0 else ''
			self.not_sat_adc_N_ch['{i}_Ch'.format(i=ch + 1)] = '(' + '&&'.join(list_sat_N_ch) + ')' if len(list_not_sat_N_ch) > 0 else ''
			self.sat_adc_N_h['{i}_H'.format(i=ch + 1)] = '(' + '||'.join(list_sat_N_h) + ')' if len(list_sat_N_h) > 0 else ''
			self.not_sat_adc_N_h['{i}_H'.format(i=ch + 1)] = '(' + '&&'.join(list_sat_N_h) + ')' if len(list_not_sat_N_h) > 0 else ''

			self.valid_ped_sigma_N_ch['{i}_Ch'.format(i=ch + 1)] = '(' + '&&'.join(list_valid_ped_sigma_N_ch) + ')' if len(list_valid_ped_sigma_N_ch) > 0 else ''
			self.valid_ped_sigma_N_h['{i}_H'.format(i=ch + 1)] = '(' + '&&'.join(list_valid_ped_sigma_N_h) + ')' if len(list_valid_ped_sigma_N_h) > 0 else ''

	def SetCells(self, selection, not_selection):
		self.selected_cells = self.selected_cells.format(s=selection)
		self.not_selected_cells = self.not_selected_cells.format(ns=not_selection)
		self.all_cells = self.all_cells.format(s=selection, ns=not_selection)

	def SetCellsCenters(self, selection, not_selection):
		self.selected_cells_centers = self.selected_cells_centers.format(s=selection)
		self.not_selected_cells_centers = self.not_selected_cells_centers.format(ns=not_selection)

	def SetUpDownBorderCuts(self, lower, upper):
		self.no_up_down_borders = self.no_up_down_borders.format(l=lower, h=upper)

	def GetThCut(self, var='clusterChargeN', th=100):
		return '({v}>={t})'.format(v=var, t=th)

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
