#!/usr/bin/env python
import numpy as np
import ROOT as ro
from Cell3D import Cell3D

class HexagonalCell3D(Cell3D):
	def __init__(self, col_num=0, row_num=0, height=0, width_pitch_ratio=0, run=0):
		Cell3D.__init__(self, col_num, row_num, sides=6, height=height, width_pitch_ratio=width_pitch_ratio, run=run)

	def CreateTCutG(self):
		tempx = np.add(self.xcenter, np.divide(np.array([-(2 * self.p - self.w), -self.w, -(2 * self.p - self.w), 2 * self.p - self.w, self.w, 2 * self.p - self.w, -(2 * self.p - self.w)], 'f8'), 2.0, dtype='f8'), dtype='f8')
		tempy = np.add(self.ycenter, np.divide(np.array([-self.h, 0, self.h, self.h, 0, -self.h, -self.h], 'f8'), 2.0, dtype='f8'), dtype='f8')
		tempname = 'cutg_dia_' + str(self.run) + '_{c}_{r}_hex'.format(c=self.col_num, r=self.row_num)
		self.cutg = ro.TCutG(tempname, 7, tempx, tempy)
		self.cutg.SetNameTitle(tempname, tempname)
		self.cutg.SetVarX('diaChXPred')
		self.cutg.SetVarY('diaChYPred')
		self.cutg.SetLineColor(ro.kBlack)

	def CreateTCutGCenter(self):
		tempx = np.add(self.xcenter, np.divide(np.array([-(2 * self.p - self.w), -self.w, -(2 * self.p - self.w), 2 * self.p - self.w, self.w, 2 * self.p - self.w, -(2 * self.p - self.w)], 'f8'), 4.0, dtype='f8'), dtype='f8')
		tempy = np.add(self.ycenter, np.divide(np.array([-self.h, 0, self.h, self.h, 0, -self.h, -self.h], 'f8'), 4.0, dtype='f8'), dtype='f8')
		tempname = 'cutg_dia_' + str(self.run) + '_center_{c}_{r}_hex'.format(c=self.col_num, r=self.row_num)
		self.cutg_center = ro.TCutG(tempname, 7, tempx, tempy)
		self.cutg_center.SetNameTitle(tempname, tempname)
		self.cutg_center.SetVarX('diaChXPred')
		self.cutg_center.SetVarY('diaChYPred')
		self.cutg_center.SetLineColor(ro.kViolet)


if __name__ == '__main__':
	z = HexagonalCell3D()