#!/usr/bin/env python
import numpy as np
import ipdb
from Utils import *

class CutManager:
	def __init__(self, tree, numstrips=5, clust_size=10, sat_adc=4095, neg_adc_th=4100, neg_snr_th=410):
		self.tree = tree
		self.numstrips = numstrips
		self.cluster_size = clust_size
		self.sat_adc = sat_adc
		self.neg_adc_th = neg_adc_th
		self.neg_snr_th = neg_snr_th
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

		self.noisy_chs = '(diaChsNoisy)'
		self.nc_chs = '(diaChsNC)'
		self.screened_chs = '(diaChsScreened)'
		self.masked_chs = self.OrCuts([self.noisy_chs, self.screened_chs, self.nc_chs])
		self.not_masked_chs = self.NotCut(self.masked_chs)
		self.hit_chs = '(diaChHits)'
		self.seed_chs = '(diaChSeed)'
		self.in_transp_cluster = self.OrCuts([self.hit_chs, self.seed_chs])
		self.not_in_transp_cluster = self.NotCut(self.in_transp_cluster)
		self.not_in_transp_cluster_not_masked = self.AndCuts([self.not_in_transp_cluster, self.not_masked_chs])

		self.any_saturated = '(diaChADC=={s})'.format(s=self.sat_adc)
		self.non_saturated = '(diaChADC!={s})'.format(s=self.sat_adc)
		self.valid_ped_sigma = '(diaChPedSigmaCmc>0)'
		self.valid_ped_friend_sigma = '(pedTree.diaChPedSigmaCmc>0)'
		self.even_pred_chs = '(TMath::Floor(diaChXPred+0.5)%2==0)'
		self.odd_pred_chs = '(TMath::Floor(diaChXPred+0.5)%2!=0)'
		self.even_cols = '(col%2==0)'
		self.odd_cols = '(col%2==1)'

		self.sat_evts = None
		self.sat_evts_region = '(satRegion)'
		self.not_sat_evts_region = '(!satRegion)'

		self.noise_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.noise_friend_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.noise_nc_cuts = {t: '' for t in ['all', 'good', 'bad']}
		self.noise_nc_friend_cuts = {t: '' for t in ['all', 'good', 'bad']}

		self.in_central_rect_region = {}
		self.out_central_rect_region = {}

	def GetValidPedChCut(self, ch, typ='Ch', isFriend=False):
		varch = 'clusterChannel{c}'.format(c=ch) if typ == 'Ch' else 'clusterChannelHighest{c}'.format(c=ch) if typ == 'H' else str(ch)
		return '(diaChPedSigmaCmc[{c}]>0)'.format(c=varch) if not isFriend else '(pedTree.diaChPedSigmaCmc[{c}]>0)'.format(c=varch)

	def GetValidPedNChsCut(self, n, typ='Ch', isFriend=False, op='&&'):
		temp = []
		for chi in xrange(1, n + 1):
			temp.append(self.GetValidPedChCut(chi if typ == 'H' else chi - 1, typ, isFriend))
		return self.AndCuts(temp) if op == '&&' else self.OrCuts(temp) if op == '||' else '(0)'

	def GetNoiseCut(self, cells, isFriend=False):
		cut1 = self.valid_ped_sigma if not isFriend else self.valid_ped_friend_sigma
		if cells == 'bad':
			return self.AndCuts([cut1, self.not_in_transp_cluster_not_masked, self.only_bad_chs_cut])
		elif cells == 'good':
			return self.AndCuts([cut1, self.not_in_transp_cluster_not_masked, self.good_chs_cut])
		else:
			return self.AndCuts([cut1, self.not_in_transp_cluster_not_masked])

	def GetSatChCut(self, ch, typ='Ch'):
		varch = 'clusterChannel' + str(ch) if typ == 'Ch' else 'clusterChannelHighest' + str(ch) if typ == 'H' else str(ch)
		return self.AndCuts([self.sat_evts_region, '(diaChADC[{c}]>={t})'.format(c=varch, t=self.sat_adc)])

	def GetSatNChsCut(self, n, typ='H', op='||'):
		temp = []
		for chi in xrange(1, n + 1):
			temp.append(self.GetSatChCut(chi if typ == 'H' else chi - 1, typ))
		return self.AndCuts(temp) if op == '&&' else self.OrCuts(temp) if op == '||' else '(0)'

	def GetNotSatChCut(self, ch, typ='Ch'):
		varch = 'clusterChannel' + str(ch) if typ == 'Ch' else 'clusterChannelHighest' + str(ch) if typ == 'H' else str(ch)
		return self.AndCuts([self.not_sat_evts_region, '(diaChADC[{c}]<{t})'.format(c=varch, t=self.sat_adc)])

	def GetNotSatNChsCut(self, n, typ='H', op='&&'):
		temp = []
		for chi in xrange(1, n + 1):
			temp.append(self.GetNotSatChCut(chi if typ == 'H' else chi - 1, typ))
		return self.AndCuts(temp) if op == '&&' else self.OrCuts(temp) if op == '||' else '(0)'

	def GetNegPHChCut(self, ch, typ='Ch', isSNR=False, isFriend=False):
		varch = 'clusterChannel' + str(ch) if typ == 'Ch' else 'clusterChannelHighest' + str(ch) if typ == 'H' else str(ch)
		cut1 = self.GetValidPedChCut(ch, typ, isFriend)
		if not isFriend:
			cut2 = '(diaChSignal[{c}]<-{t})'.format(c=varch, t=self.neg_adc_th) if not isSNR else '(diaChSignal[{c}]<-{t}*diaChPedSigmaCmc[{c}])'.format(c=varch, t=self.neg_snr_th)
		else:
			cut2 = '(pedTree.diaChSignal[{c}]<-{t})'.format(c=varch, t=self.neg_adc_th) if not isSNR else '(pedTree.diaChSignal[{c}]<-{t}*pedTree.diaChPedSigmaCmc[{c}])'.format(c=varch, t=self.neg_snr_th)
		return self.AndCuts([cut1, cut2])

	def GetNegPHNChsCut(self, n, typ='Ch', isSNR=False, isFriend=False, op='||'):
		temp = []
		for chi in xrange(1, n + 1):
			temp.append(self.GetNegPHChCut(chi if typ == 'H' else chi - 1, typ, isSNR, isFriend))
		return self.AndCuts(temp) if op == '&&' else self.OrCuts(temp) if op == '||' else '(0)'

	def GetNotNegPHChCut(self, ch, typ='Ch', isSNR=False, isFriend=False):
		varch = 'clusterChannel' + str(ch) if typ == 'Ch' else 'clusterChannelHighest' + str(ch) if typ == 'H' else str(ch)
		cut1 = self.GetValidPedChCut(ch, typ, isFriend)
		if not isFriend:
			cut2 = '(diaChSignal[{c}]>=-{t})'.format(c=varch, t=self.neg_adc_th) if not isSNR else '(diaChSignal[{c}]>=-{t}*diaChPedSigmaCmc[{c}])'.format(c=varch, t=self.neg_snr_th)
		else:
			cut2 = '(pedTree.diaChSignal[{c}]>=-{t})'.format(c=varch, t=self.neg_adc_th) if not isSNR else '(pedTree.diaChSignal[{c}]>=-{t}*pedTree.diaChPedSigmaCmc[{c}])'.format(c=varch, t=self.neg_snr_th)
		return self.AndCuts([cut1, cut2])

	def GetNotNegPHNChsCut(self, n, typ='Ch', isSNR=False, isFriend=False, op='&&'):
		temp = []
		for chi in xrange(1, n + 1):
			temp.append(self.GetNotNegPHChCut(chi if typ == 'H' else chi - 1, typ, isSNR, isFriend))
		return self.AndCuts(temp) if op == '&&' else self.OrCuts(temp) if op == '||' else '(0)'

	def GetPHCuts(self, varKey, isFriend=False):
		if varKey.startswith('PH_Ch'):
			return self.GetValidPedChCut(varKey.strip('PH_Ch'), 'Ch',isFriend)
		elif varKey.startswith('PH_H'):
			return self.GetValidPedChCut(varKey.strip('PH_H'), 'H', isFriend)
		elif varKey.endswith('Ch'):
			return self.GetValidPedNChsCut(int(varKey.strip('PH').strip('_Ch')), 'Ch', isFriend, '&&')
		elif varKey.endswith('H'):
			return self.GetValidPedNChsCut(int(varKey.strip('PH').strip('_H')), 'H', isFriend, '&&')


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
		self.noise_cuts = {'all': self.AndCuts([self.valid_ped_sigma, self.not_in_transp_cluster_not_masked]), 'good': self.AndCuts([self.valid_ped_sigma, self.not_in_transp_cluster_not_masked, self.good_chs_cut]), 'bad': self.AndCuts([self.valid_ped_sigma, self.not_in_transp_cluster_not_masked, self.only_bad_chs_cut])}
		self.noise_friend_cuts = {'all': self.AndCuts([self.valid_ped_friend_sigma, self.not_in_transp_cluster_not_masked]), 'good': self.AndCuts([self.valid_ped_friend_sigma, self.not_in_transp_cluster_not_masked, self.good_chs_cut]), 'bad': self.AndCuts([self.valid_ped_friend_sigma, self.not_in_transp_cluster_not_masked, self.only_bad_chs_cut])}
		self.noise_nc_cuts = {cells: self.ConcatenateCutWithCells(cut=self.AndCuts([self.nc_chs, self.valid_ped_sigma]), cells=cells) for cells in ['all', 'good', 'bad']}
		self.noise_nc_friend_cuts = {cells: self.ConcatenateCutWithCells(cut=self.AndCuts([self.nc_chs, self.valid_ped_friend_sigma]), cells=cells) for cells in ['all', 'good', 'bad']}

	def SetUpDownBorderCuts(self, lower, upper):
		self.no_up_down_borders = self.no_up_down_borders.format(l=lower, h=upper)

	def GetThCut(self, var, th, cells, cuts='', op='>='):
		temp_cuts = '({v}{o}{th})'.format(v=var, o=op, th=th)
		temp_cuts = temp_cuts if cuts == '' else self.AndCuts([cuts, temp_cuts])
		return self.ConcatenateCutWithCells(temp_cuts, cells, '&&')

	def FindSaturationEvents(self):
		tempsat = self.AndCuts([self.transp_ev, self.any_saturated])
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

	def AndCuts(self, cuts=[]):
		return self.ConcatenateCuts(cuts, '&&')

	def OrCuts(self, cuts=[]):
		return self.ConcatenateCuts(cuts, '||')

	def ConcatenateCuts(self, cuts=[], operator='&&'):
		if len(cuts) == 0:
			return '(0)'
		elif len(cuts) == 1:
			return cuts[0]
		else:
			temp = []
			for cut in cuts:
				if cut != '':
					temp.append(cut)
			return '(' + operator.join(temp) + ')' if len(temp) > 0 else '(0)'

	def ConcatenateCutWithCells(self, cut, cells='all', operator='&&'):
		return self.ConcatenateCuts([self.selected_cells, cut], operator) if cells == 'good' else self.ConcatenateCuts([self.not_selected_cells, cut], operator) if cells == 'bad' else self.ConcatenateCuts([self.all_cells, cut], operator) if cells == 'all' else cut

	def NotCut(self, cut):
		return '(!{c})'.format(c=cut)

if __name__ == '__main__':
	cm = CutManager(None)
