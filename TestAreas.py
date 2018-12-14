#!/usr/bin/env python
import ROOT as ro
import numpy as np
import ipdb, os
from TransparentGrid import TransparentGrid
from optparse import OptionParser

tests = range(0, 10) + [100]

class TestAreas:
	def __init__(self, num=0, clust_size=1, dir='.', run=0, cellsize=50):
		self.num = num
		self.clust_size = clust_size
		self.dir = dir
		self.run = run
		self.cellsize = cellsize
		self.trans_grid = TransparentGrid(dir, run, cellsize)
		self.threshold = 800 if cellsize == 50 else 200
		self.window_shift = 10

	def SetTest(self):
		self.trans_grid.cell_resolution = 50.0 / 13.0 if self.trans_grid.col_pitch == 50 else 100.0 / 51
		if self.num == 0:
			self.trans_grid.cell_resolution = 50.0 / 25.0
			self.trans_grid.ResetAreas()
			self.trans_grid.AddGoodAreasCol(2 + 1 - self.clust_size, 17, 20)
			self.trans_grid.AddGoodAreasCol(3 + 1 - self.clust_size, 7, 9)
			self.trans_grid.AddGoodAreasCol(3 + 1 - self.clust_size, 17, 19)
			self.trans_grid.AddGoodAreasCol(4 + 1 - self.clust_size, 6, 13)
			self.trans_grid.AddGoodAreasCol(4 + 1 - self.clust_size, 16, 20)
			# self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 8, 19)
			self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 8, 16)
			self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 7, 9)
			# self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 11, 18)
			self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 11, 17)
			# self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 8, 8)
			# self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 12, 12)
			self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 14, 19)
			self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 21, 21)
			self.trans_grid.AddGoodAreasCol(8 + 1 - self.clust_size, 15, 21)
			self.trans_grid.AddGoodAreasCol(8 + 1 - self.clust_size, 24, 24)
			self.trans_grid.AddGoodAreasCol(9 + 1 - self.clust_size, 12, 12)
			self.trans_grid.AddGoodAreasCol(9 + 1 - self.clust_size, 14, 17)
			self.trans_grid.AddGoodAreasCol(9 + 1 - self.clust_size, 19, 25)
			self.trans_grid.AddGoodAreasCol(10 + 1 - self.clust_size, 12, 17)
			self.trans_grid.AddGoodAreasCol(10 + 1 - self.clust_size, 19, 19)
			self.trans_grid.AddGoodAreasCol(11 + 1 - self.clust_size, 13, 17)
			self.trans_grid.AddGoodAreasCol(11 + 1 - self.clust_size, 19, 20)
			self.trans_grid.RemoveFromGoodArea(3 + 1 - self.clust_size, 17)
			self.trans_grid.RemoveFromGoodArea(5 + 1 - self.clust_size, 17)
			self.trans_grid.RemoveFromGoodArea(5 + 1 - self.clust_size, 18)
			self.trans_grid.RemoveFromGoodArea(6 + 1 - self.clust_size, 7)
			self.trans_grid.RemoveFromGoodArea(7 + 1 - self.clust_size, 14)
			self.trans_grid.RemoveFromGoodArea(9 + 1 - self.clust_size, 16)
			self.trans_grid.RemoveFromGoodArea(9 + 1 - self.clust_size, 17)
			self.trans_grid.RemoveFromGoodArea(10 + 1 - self.clust_size, 17)
			self.trans_grid.RemoveFromGoodArea(10 + 1 - self.clust_size, 19)
			self.trans_grid.AddRemainingToBadAreas()
		elif self.num == 1:
			self.trans_grid.cell_resolution = 50.0 / 11.0
			# self.trans_grid.cell_resolution = 50.0 / 17.0
			self.trans_grid.ResetAreas()
			# self.trans_grid.AddGoodAreas(2 + 1 - self.clust_size, 8)
			self.trans_grid.AddGoodAreasCol(3 + 1 - self.clust_size, 7, 9)
			self.trans_grid.AddGoodAreasCol(4 + 1 - self.clust_size, 6, 13)
			self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 8, 19)
			# self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 11, 18)
			self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 11, 16)
			# self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 14, 19)
			self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 14, 17)
			# self.trans_grid.AddGoodAreasCol(8 + 1 - self.clust_size, 15, 21)
			self.trans_grid.AddGoodAreasCol(8 + 1 - self.clust_size, 15, 17)
			self.trans_grid.AddGoodAreasCol(9 + 1 - self.clust_size, 14, 17)
			self.trans_grid.AddGoodAreasCol(10 + 1 - self.clust_size, 12, 17)
			self.trans_grid.AddGoodAreasCol(11 + 1 - self.clust_size, 13, 17)
			self.trans_grid.RemoveFromGoodArea(5 + 1 - self.clust_size, 17)
			self.trans_grid.RemoveFromGoodArea(5 + 1 - self.clust_size, 18)
			self.trans_grid.RemoveFromGoodArea(5 + 1 - self.clust_size, 19)
			self.trans_grid.RemoveFromGoodArea(6 + 1 - self.clust_size, 18)
			self.trans_grid.RemoveFromGoodArea(7 + 1 - self.clust_size, 14)
			self.trans_grid.RemoveFromGoodArea(9 + 1 - self.clust_size, 16)
			self.trans_grid.RemoveFromGoodArea(9 + 1 - self.clust_size, 17)
			self.trans_grid.RemoveFromGoodArea(10 + 1 - self.clust_size, 16)
			self.trans_grid.RemoveFromGoodArea(10 + 1 - self.clust_size, 17)
			self.trans_grid.AddRemainingToBadAreas()
		elif self.num == 2:
			self.trans_grid.cell_resolution = 50.0 / 15.0
			self.trans_grid.ResetAreas()
			self.trans_grid.AddGoodAreasCol(4 + 1 - self.clust_size, 6, 13)
			# self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 8, 19)
			self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 8, 16)
			# self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 11, 18)
			self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 11, 17)
			# self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 14, 19)
			self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 15, 19)
			self.trans_grid.AddGoodAreasCol(8 + 1 - self.clust_size, 15, 21)
			self.trans_grid.AddGoodAreasCol(8 + 1 - self.clust_size, 15, 21)
			self.trans_grid.AddRemainingToBadAreas()
		elif self.num == 3:
			self.trans_grid.cell_resolution = 50.0 / 15.0
			self.trans_grid.ResetAreas()
			self.trans_grid.AddGoodAreasCol(4 + 1 - self.clust_size, 8, 19)
			self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 8, 19)
			self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 8, 19)
			self.trans_grid.AddRemainingToBadAreas()
		elif self.num == 4:
			self.trans_grid.ResetAreas()
			# self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 8, 19)
			self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 8, 16)
			# self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 8, 19)
			self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 8, 16)
			self.trans_grid.AddRemainingToBadAreas()
		elif self.num == 5:
			self.trans_grid.ResetAreas()
			self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 14, 19)
			self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 14, 19)
			self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 14, 19)
			self.trans_grid.AddGoodAreasCol(8 + 1 - self.clust_size, 14, 19)
			self.trans_grid.AddGoodAreasCol(9 + 1 - self.clust_size, 14, 19)
			self.trans_grid.RemoveFromGoodArea(5 + 1 - self.clust_size, 17)
			self.trans_grid.RemoveFromGoodArea(5 + 1 - self.clust_size, 18)
			self.trans_grid.RemoveFromGoodArea(6 + 1 - self.clust_size, 18)
			self.trans_grid.RemoveFromGoodArea(7 + 1 - self.clust_size, 14)
			self.trans_grid.RemoveFromGoodArea(9 + 1 - self.clust_size, 16)
			self.trans_grid.RemoveFromGoodArea(9 + 1 - self.clust_size, 17)
			self.trans_grid.RemoveFromGoodArea(9 + 1 - self.clust_size, 18)
			self.trans_grid.AddRemainingToBadAreas()
		elif self.num == 6:
			self.trans_grid.cell_resolution = 50.0 / 11.0
			self.trans_grid.ResetAreas()
			self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 14, 18)
			self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 14, 18)
			self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 14, 18)
			self.trans_grid.RemoveFromGoodArea(5 + 1 - self.clust_size, 17)
			self.trans_grid.RemoveFromGoodArea(5 + 1 - self.clust_size, 18)
			self.trans_grid.RemoveFromGoodArea(6 + 1 - self.clust_size, 18)
			self.trans_grid.RemoveFromGoodArea(7 + 1 - self.clust_size, 14)
			self.trans_grid.AddRemainingToBadAreas()
		elif self.num == 7:
			self.trans_grid.cell_resolution = 50.0 / 11.0
			self.trans_grid.ResetAreas()
			# self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 15, 18)
			self.trans_grid.AddGoodAreasCol(5 + 1 - self.clust_size, 15, 16)
			# self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 15, 18)
			self.trans_grid.AddGoodAreasCol(6 + 1 - self.clust_size, 15, 17)
			self.trans_grid.AddGoodAreasCol(7 + 1 - self.clust_size, 15, 18)
			self.trans_grid.AddGoodAreasCol(8 + 1 - self.clust_size, 15, 18)
			self.trans_grid.AddRemainingToBadAreas()
		elif self.num == 8:
			self.trans_grid.cell_resolution = 50.0 / 11.0
			self.trans_grid.ResetAreas()
			self.trans_grid.AddGoodAreasRow(15, 5 + 1 - self.clust_size, 11 + 2 - 2 * self.clust_size)
			self.trans_grid.AddGoodAreasRow(16, 5 + 1 - self.clust_size, 11 + 2 - 2 * self.clust_size)
			self.trans_grid.AddGoodAreasRow(17, 5 + 1 - self.clust_size, 11 + 2 - 2 * self.clust_size)
			self.trans_grid.AddRemainingToBadAreas()
		elif self.num == 9:
			self.trans_grid.cell_resolution = 2
			self.trans_grid.ResetAreas()
			self.trans_grid.SelectGoodAndBadByThreshold(self.threshold)
		elif self.num == 100:
			self.trans_grid.cell_resolution = 2
			self.trans_grid.ResetAreas()
			self.trans_grid.SelectGoodAndBadByThreshold(self.threshold)
			
	def PlotTest(self):
		w = 0
		num = self.num
		y0, rowpitch, numrows, xoff, yoff, colpitch, numcols, yup = self.trans_grid.row_info_diamond['0'], self.trans_grid.row_info_diamond['pitch'], self.trans_grid.row_info_diamond['num'], self.trans_grid.row_info_diamond['x_off'], self.trans_grid.row_info_diamond['y_off'], self.trans_grid.col_pitch, self.trans_grid.num_cols, self.trans_grid.row_info_diamond['up']
		for ch in xrange(1, self.clust_size + 1):
			#  plot map
			self.trans_grid.DrawProfile2DDiamond('ph{c}_map_test{n}'.format(c=ch, n=num), varz='clusterCharge' + str(ch), cuts='({l}<=diaChYPred)&&(diaChYPred<={h})'.format(l=y0, h=yup))
			# ro.gPad.Update()
			self.trans_grid.canvas['ph{c}_map_test{n}'.format(c=ch, n=num)].SetWindowPosition(w, w)
			w += self.window_shift
			self.trans_grid.gridTextDiamond.Draw('same TEXT0')
			self.trans_grid.profile['ph{c}_map_test'.format(c=ch) + str(num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile['ph{c}_map_test'.format(c=ch) + str(num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
			#  draw selected areas
			self.trans_grid.DrawGoodAreasDiamond('ph{c}_map_test{n}'.format(c=ch, n=num))
			self.trans_grid.DrawGoodAreasDiamondCenters('ph{c}_map_test{n}'.format(c=ch, n=num))
			#  plot cell overlay
			self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{c}_cells_test{n}'.format(c=ch, n=num), var='clusterCharge' + str(ch), cells='good')
			self.trans_grid.canvas['ph{c}_cells_test{n}'.format(c=ch, n=num)].SetWindowPosition(w, w)
			w += self.window_shift
			#  show center region in cell
			self.trans_grid.DrawTCutCentersInCellOverlay('ph{c}_cells_test{n}'.format(c=ch, n=num))
			#  draw ph of selected areas
			self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}'.format(c=ch, n=num), var='clusterCharge' + str(ch))
			self.trans_grid.canvas['ph{c}_test{n}'.format(c=ch, n=num)].SetWindowPosition(w, w)
			w += self.window_shift
			self.trans_grid.DrawPHCentralRegion('ph{c}_test{n}_centers'.format(c=ch, n=num), cells='good', var='clusterCharge' + str(ch))
			self.trans_grid.canvas['ph{c}_test{n}_centers'.format(c=ch, n=num)].SetWindowPosition(w, w)
			w += self.window_shift
			#  fit distribution for central region
			self.trans_grid.FitLanGaus('ph{c}_test{n}_centers'.format(c=ch, n=num), color=ro.kRed)
			#  get difference between cell and center
			self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}_periphery'.format(c=ch, n=num))
			self.trans_grid.histo['ph{c}_test{n}_periphery'.format(c=ch, n=num)].Reset('ICES')
			self.trans_grid.histo['ph{c}_test{n}_periphery'.format(c=ch, n=num)].Add(self.trans_grid.histo['ph{c}_test{n}'.format(c=ch, n=num)], self.trans_grid.histo['ph{c}_test{n}_centers'.format(c=ch, n=num)], 1, -1)
			self.trans_grid.canvas['ph{c}_test{n}_periphery'.format(c=ch, n=num)].SetWindowPosition(w, w)
			w += self.window_shift
			self.trans_grid.FitLanGaus('ph{c}_test{n}_periphery'.format(c=ch, n=num), color=ro.kBlue)
			self.trans_grid.canvas['ph{c}_test{n}'.format(c=ch, n=num)].cd()
			self.trans_grid.langaus['ph{c}_test{n}_centers'.format(c=ch, n=num)].fit.Draw('same')
			self.trans_grid.langaus['ph{c}_test{n}_periphery'.format(c=ch, n=num)].fit.Draw('same')
			self.trans_grid.DrawDoubleLangaus('ph{c}_test{n}'.format(c=ch, n=num), 'ph{c}_test{n}_centers'.format(c=ch, n=num), 'ph{c}_test{n}_periphery'.format(c=ch, n=num), color=ro.kBlack)
			# ro.gPad.Update()
			#  position of negative clusters
			self.trans_grid.Draw2DHistoDiamond('hm{c}_test{n}_negative'.format(c=ch, n=num))
			self.trans_grid.histo['hm{c}_test{n}_negative'.format(c=ch, n=num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.histo['hm{c}_test{n}_negative'.format(c=ch, n=num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
			self.trans_grid.canvas['hm{c}_test{n}_negative'.format(c=ch, n=num)].SetWindowPosition(w, w)
			w += self.window_shift
			self.trans_grid.DrawGoodAreasDiamond('hm{c}_test{n}_negative'.format(c=ch, n=num))
			self.trans_grid.DrawBadAreasDiamond('hm{c}_test{n}_negative'.format(c=ch, n=num))
			self.trans_grid.DrawGoodAreasDiamondCenters('hm{c}_test{n}_negative'.format(c=ch, n=num))
			self.trans_grid.DrawProfile2DDiamond('ph{c}_map_test{n}_negative'.format(c=ch, n=num), varz='clusterCharge' + str(ch),
			                          cuts='((diaChSignal[int(TMath::Floor(diaChXPred+0.5))]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))-1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))+1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c}))'.format(
				                          c=self.trans_grid.neg_cut))
			self.trans_grid.profile['ph{c}_map_test{n}_negative'.format(c=ch, n=num)].GetXaxis().SetRangeUser(self.trans_grid.ch_ini - 1, self.trans_grid.ch_end + 1)
			self.trans_grid.profile['ph{c}_map_test{n}_negative'.format(c=ch, n=num)].GetYaxis().SetRangeUser(y0 - int(rowpitch / self.trans_grid.bins_per_ch_y), yup + int(rowpitch / self.trans_grid.bins_per_ch_y))
			self.trans_grid.canvas['ph{c}_map_test{n}_negative'.format(c=ch, n=num)].SetWindowPosition(w, w)
			w += self.window_shift
			self.trans_grid.DrawGoodAreasDiamond('ph{c}_map_test{n}_negative'.format(c=ch, n=num))
			self.trans_grid.DrawBadAreasDiamond('ph{c}_map_test{n}_negative'.format(c=ch, n=num))
			self.trans_grid.DrawGoodAreasDiamondCenters('ph{c}_map_test{n}_negative'.format(c=ch, n=num))
			self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{c}_cells_test{n}_negative'.format(c=ch, n=num), var='clusterCharge' + str(ch), cells='good',
			                                     cuts='((diaChSignal[int(TMath::Floor(diaChXPred+0.5))]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))-1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))-1]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))+1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))+1]*{c}))'.format(
				                                     c=self.trans_grid.neg_cut))
			self.trans_grid.canvas['ph{c}_cells_test{n}_negative'.format(c=ch, n=num)].SetWindowPosition(w, w)
			w += self.window_shift
			self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}_negative'.format(c=ch, n=num), var='clusterCharge' + str(ch),
			                     cuts='((diaChSignal[int(TMath::Floor(diaChXPred+0.5))]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))-1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))-1]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))+1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))+1]*{c}))'.format(
				                     c=self.trans_grid.neg_cut))
			self.trans_grid.canvas['ph{c}_test{n}_negative'.format(c=ch, n=num)].SetWindowPosition(w, w)
			w += self.window_shift
			self.trans_grid.FitLanGaus('ph{c}_test{n}_negative'.format(c=ch, n=num), color=ro.kBlue)
			self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{c}_cells_test{n}_no_negative'.format(c=ch, n=num), var='clusterCharge' + str(ch), cells='good',
			                                     cuts='((diaChSignal[int(TMath::Floor(diaChXPred+0.5))]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})&&(diaChSignal[int(TMath::Floor(diaChXPred+0.5))-1]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))-1]*{c})&&(diaChSignal[int(TMath::Floor(diaChXPred+0.5))+1]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))+1]*{c}))'.format(
				                                     c=self.trans_grid.neg_cut))
			self.trans_grid.canvas['ph{c}_cells_test{n}_no_negative'.format(c=ch, n=num)].SetWindowPosition(w, w)
			w += self.window_shift
			self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}_no_negative'.format(c=ch, n=num), var='clusterCharge' + str(ch),
			                     cuts='((diaChSignal[int(TMath::Floor(diaChXPred+0.5))]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})&&(diaChSignal[int(TMath::Floor(diaChXPred+0.5))-1]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))-1]*{c})&&(diaChSignal[int(TMath::Floor(diaChXPred+0.5))+1]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))+1]*{c}))'.format(
				                     c=self.trans_grid.neg_cut))
			self.trans_grid.canvas['ph{c}_test{n}_no_negative'.format(c=ch, n=num)].SetWindowPosition(w, w)
			w += self.window_shift
			self.trans_grid.FitLanGaus('ph{c}_test{n}_no_negative'.format(c=ch, n=num), color=ro.kRed)
		if self.clust_size >= 2:
			self.trans_grid.DrawProfile2DDiamondCellOverlay('ph2_minus_ph1_map_test{n}'.format(n=num), 'clusterCharge2-clusterCharge1', cells='good')
			self.trans_grid.canvas['ph2_minus_ph1_map_test{n}'.format(n=num)].SetWindowPosition(w, w)
			w += self.window_shift
			self.trans_grid.DrawTCutCentersInCellOverlay('ph2_minus_ph1_map_test{n}'.format(n=num))
			tempmin, tempmax, tempbins = self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins
			self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins = self.trans_grid.phmin_neg, self.trans_grid.phmax_neg, self.trans_grid.phbins_neg
			self.trans_grid.DrawPHGoodAreas('ph2_minus_ph1_test{n}'.format(n=num), 'clusterCharge2-clusterCharge1')
			self.trans_grid.canvas['ph2_minus_ph1_test{n}'.format(n=num)].SetWindowPosition(w, w)
			w += self.window_shift

			self.trans_grid.phmin, self.trans_grid.phmax, self.trans_grid.phbins = tempmin, tempmax, tempbins
			self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}'.format(c='N', n=num), var='clusterCharge' + str('N'))
			self.trans_grid.canvas['ph{c}_test{n}'.format(c='N', n=num)].SetWindowPosition(w, w)
			w += self.window_shift
			self.trans_grid.DrawPHCentralRegion('ph{c}_test{n}_centers'.format(c='N', n=num), cells='good', var='clusterCharge' + str('N'))
			self.trans_grid.canvas['ph{c}_test{n}_centers'.format(c='N', n=num)].SetWindowPosition(w, w)
			w += self.window_shift
			#  fit distribution for central region
			self.trans_grid.FitLanGaus('ph{c}_test{n}_centers'.format(c='N', n=num), color=ro.kRed)
			#  get difference between cell and center
			self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}_periphery'.format(c='N', n=num))
			self.trans_grid.histo['ph{c}_test{n}_periphery'.format(c='N', n=num)].Reset('ICES')
			self.trans_grid.histo['ph{c}_test{n}_periphery'.format(c='N', n=num)].Add(self.trans_grid.histo['ph{c}_test{n}'.format(c='N', n=num)], self.trans_grid.histo['ph{c}_test{n}_centers'.format(c='N', n=num)], 1, -1)
			self.trans_grid.canvas['ph{c}_test{n}_periphery'.format(c='N', n=num)].SetWindowPosition(w, w)
			w += self.window_shift
			self.trans_grid.FitLanGaus('ph{c}_test{n}_periphery'.format(c='N', n=num), color=ro.kBlue)
			self.trans_grid.canvas['ph{c}_test{n}'.format(c='N', n=num)].cd()
			self.trans_grid.langaus['ph{c}_test{n}_centers'.format(c='N', n=num)].fit.Draw('same')
			self.trans_grid.langaus['ph{c}_test{n}_periphery'.format(c='N', n=num)].fit.Draw('same')
			self.trans_grid.DrawDoubleLangaus('ph{c}_test{n}'.format(c='N', n=num), 'ph{c}_test{n}_centers'.format(c='N', n=num), 'ph{c}_test{n}_periphery'.format(c='N', n=num), color=ro.kBlack)

			# self.trans_grid.FitLanGaus('ph2_minus_ph1_test{n}'.format(n=num), color=ro.kViolet)
			self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{c}_cells_test{n}_negative'.format(c=3, n=num), var='clusterChargeN', cells='good',
			                                     cuts='((diaChSignal[int(TMath::Floor(diaChXPred+0.5))]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))-1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))-1]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))+1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))+1]*{c}))'.format(
				                                     c=self.trans_grid.neg_cut))
			self.trans_grid.canvas['ph{c}_cells_test{n}_negative'.format(c=3, n=num)].SetWindowPosition(w, w)
			w += self.window_shift
			self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}_negative'.format(c=3, n=num), var='clusterChargeN',
			                     cuts='((diaChSignal[int(TMath::Floor(diaChXPred+0.5))]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))-1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))-1]*{c})||(diaChSignal[int(TMath::Floor(diaChXPred+0.5))+1]<-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))+1]*{c}))'.format(
				                     c=self.trans_grid.neg_cut))
			self.trans_grid.canvas['ph{c}_test{n}_negative'.format(c=3, n=num)].SetWindowPosition(w, w)
			self.trans_grid.FitLanGaus('ph{c}_test{n}_negative'.format(c=3, n=num), color=ro.kBlue)
			w += self.window_shift
			self.trans_grid.DrawProfile2DDiamondCellOverlay('ph{c}_cells_test{n}_no_negative'.format(c=3, n=num), var='clusterChargeN', cells='good',
			                                     cuts='((diaChSignal[int(TMath::Floor(diaChXPred+0.5))]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})&&(diaChSignal[int(TMath::Floor(diaChXPred+0.5))-1]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))-1]*{c})&&(diaChSignal[int(TMath::Floor(diaChXPred+0.5))+1]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))+1]*{c}))'.format(
				                                     c=self.trans_grid.neg_cut))
			self.trans_grid.canvas['ph{c}_cells_test{n}_no_negative'.format(c=3, n=num)].SetWindowPosition(w, w)
			w += self.window_shift
			self.trans_grid.DrawPHGoodAreas('ph{c}_test{n}_no_negative'.format(c=3, n=num), var='clusterChargeN',
			                     cuts='((diaChSignal[int(TMath::Floor(diaChXPred+0.5))]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))]*{c})&&(diaChSignal[int(TMath::Floor(diaChXPred+0.5))-1]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))-1]*{c})&&(diaChSignal[int(TMath::Floor(diaChXPred+0.5))+1]>-diaChPedSigmaCmc[int(TMath::Floor(diaChXPred+0.5))+1]*{c}))'.format(
				                     c=self.trans_grid.neg_cut))
			self.trans_grid.canvas['ph{c}_test{n}_no_negative'.format(c=3, n=num)].SetWindowPosition(w, w)
			w += self.window_shift
			self.trans_grid.FitLanGaus('ph{c}_test{n}_no_negative'.format(c=3, n=num), color=ro.kRed)

	def SaveCanvas(self):
		if not os.path.isdir('{d}/{r}/test{n}'.format(d=self.trans_grid.dir, r=self.trans_grid.run, n=self.num)):
			os.makedirs('{d}/{r}/test{n}'.format(d=self.trans_grid.dir, r=self.trans_grid.run, n=self.num))
		self.trans_grid.SaveCanvasInlist(self.trans_grid.canvas.keys(), 'test{n}'.format(n=self.num))

if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-d', '--dir', dest='dir', type='string', help='Path to the subdirectory that contains the output of different runs')
	parser.add_option('-r', '--run', dest='run', type='int', help='run number to be analysed (e.g. 25209)')
	parser.add_option('-c', '--cellsize', dest='cellsize', type='int', default=50, help='cell size of the square 3D device')
	parser.add_option('-t', '--test', dest='testnumber', type='int', default=-1, help='Run a automatically one of the predefined tests')
	parser.add_option('-n', '--numstrips', dest='numstrips', type='int', default=2, help='Number of strips to use')
	parser.add_option('-a', '--auto', dest='auto', default=False, action='store_true', help='Sets up test, creates plots and saves them automatically if toggled')

	(options, args) = parser.parse_args()
	run = int(options.run)
	dir = str(options.dir)
	testnum = int(options.testnumber)
	numstrips = int(options.numstrips)
	cellsize = int(options.cellsize) if testnum < 100 else 100
	autom = bool(options.auto)

	t = TestAreas(testnum, numstrips, dir, run, cellsize)
	t.trans_grid.SetLines()
	t.trans_grid.CreateTCutGs()
	if testnum in tests:
		t.SetTest()
		if autom:
			t.PlotTest()
			t.SaveCanvas()

