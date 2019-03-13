import numpy as np
import ipdb
from Utils import *

class CutManager:
	def __init__(self, tree, numstrips=5, clust_size=10, sat_adc=4095):
		self.tree = tree
		self.numstrips = numstrips
		self.clust_size = clust_size
		self.sat_adc = 4095
		self.transp_ev = '(transparentEvent)'
		self.tcutg_cell = {}
		self.selected_cells = ''
		self.not_selected_cells = ''
		self.all_cells = '(({s})||({ns}))'
		self.selected_cells_centers = ''
		self.not_selected_cells_centers = ''
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

		self.InitializeNegCuts()

	def InitializeNegCuts(self):
		for i in xrange(self.clust_size):
			if FindLeafInTree(self.tree, 'clusterChannel{i}'.format(i=i)):
				self.neg_snr_ph_ch['PH_Ch{i}'.format(i=i)] = '((diaChPedSigmaCmc[clusterChannel{i}]>0)&&(diaChSignal[clusterChannel{i}]/diaChPedSigmaCmc[clusterChannel{i}]<-{th}))'.format(i=i, th='{th}')
				self.not_neg_snr_ph_ch['PH_Ch{i}'.format(i=i)] = '((diaChPedSigmaCmc[clusterChannel{i}]>0)&&(diaChSignal[clusterChannel{i}]/diaChPedSigmaCmc[clusterChannel{i}]>=-{th}))'.format(i=i, th='{th}')
				self.neg_adc_ph_ch['PH_Ch{i}'.format(i=i)] = '(diaChSignal[clusterChannel{i}]<-{th})'.format(i=i, th='{th}')
				self.not_neg_adc_ph_ch['PH_Ch{i}'.format(i=i)] = '(diaChSignal[clusterChannel{i}]>=-{th})'.format(i=i, th='{th}')
			if FindLeafInTree(self.tree, 'clusterChannelHighest{i}'.format(i=i+1)):
				self.neg_snr_ph_h['PH_H{i}'.format(i=i+1)] = '((diaChPedSigmaCmc[clusterChannelHighest{i}]>0)&&(diaChSignal[clusterChannelHighest{i}]/diaChPedSigmaCmc[clusterChannelHighest{i}]<-{th}))'.format(i=i+1, th='{th}')
				self.not_neg_snr_ph_h['PH_H{i}'.format(i=i+1)] = '((diaChPedSigmaCmc[clusterChannelHighest{i}]>0)&&(diaChSignal[clusterChannelHighest{i}]/diaChPedSigmaCmc[clusterChannelHighest{i}]>=-{th}))'.format(i=i+1, th='{th}')
				self.neg_adc_ph_h['PH_H{i}'.format(i=i+1)] = '(diaChSignal[clusterChannelHighest{i}]<-{th})'.format(i=i+1, th='{th}')
				self.not_neg_adc_ph_h['PH_H{i}'.format(i=i+1)] = '(diaChSignal[clusterChannelHighest{i}]>=-{th})'.format(i=i+1, th='{th}')

		for ch in xrange(self.numstrips):
			list_neg_snr_phN_ch = []
			list_not_neg_snr_phN_ch = []
			list_neg_adc_phN_ch = []
			list_not_neg_adc_phN_ch = []
			list_neg_snr_phN_h = []
			list_not_neg_snr_phN_h = []
			list_neg_adc_phN_h = []
			list_not_neg_adc_phN_h = []
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
					list_not_neg_snr_phN_h.append(self.not_neg_snr_ph_h['PH_Ch{i}'.format(i=chi+1)])
					list_not_neg_adc_phN_h.append(self.not_neg_adc_ph_h['PH_Ch{i}'.format(i=chi+1)])
				self.neg_snr_phN_ch['PH{i}_Ch'.format(i=ch + 1)] = '(' + '||'.join(list_neg_snr_phN_ch) + ')' if len(list_neg_snr_phN_ch) > 0 else ''
				self.neg_adc_phN_ch['PH{i}_Ch'.format(i=ch + 1)] = '(' + '||'.join(list_neg_adc_phN_ch) + ')' if len(list_neg_adc_phN_ch) > 0 else ''
				self.not_neg_snr_phN_ch['PH{i}_Ch'.format(i=ch + 1)] = '(' + '&&'.join(list_not_neg_snr_phN_ch) + ')' if len(list_not_neg_snr_phN_ch) > 0 else ''
				self.not_neg_adc_phN_ch['PH{i}_Ch'.format(i=ch + 1)] = '(' + '&&'.join(list_not_neg_adc_phN_ch) + ')' if len(list_not_neg_adc_phN_ch) > 0 else ''
				self.neg_snr_phN_h['PH{i}_H'.format(i=ch + 1)] = '(' + '||'.join(list_neg_snr_phN_h) + ')' if len(list_neg_snr_phN_h) > 0 else ''
				self.neg_adc_phN_h['PH{i}_H'.format(i=ch + 1)] = '(' + '||'.join(list_neg_adc_phN_h) + ')' if len(list_neg_adc_phN_h) > 0 else ''
				self.not_neg_snr_phN_h['PH{i}_H'.format(i=ch + 1)] = '(' + '&&'.join(list_not_neg_snr_phN_h) + ')' if len(list_not_neg_snr_phN_h) > 0 else ''
				self.not_neg_adc_phN_h['PH{i}_H'.format(i=ch + 1)] = '(' + '&&'.join(list_not_neg_adc_phN_h) + ')' if len(list_not_neg_adc_phN_h) > 0 else ''


	def GetThCut(self, var='clusterChargeN', th=100):
		return '({v}>={t})'.format(v=var, t=th)

