#!/usr/bin/env python
import numpy as np
import ROOT as ro
from Cell3D import Cell3D

class RectangularCell3D(Cell3D):
	def __init__(self, col_num=0, row_num=0, height=0, run=0):
		Cell3D.__init__(self, col_num, row_num, sides=4, run=run)
		self.p = 1.0
		self.h = height
		self.cutg = None

	def CreateTCutG(self):
		tempx = np.add(self.xcenter, np.divide(np.array([-self.p, -self.p, self.p, self.p, -self.p], 'f8'), 2.0, dtype='f8'), dtype='f8')
		tempy = np.add(self.ycenter, np.divide(np.array([-self.h, self.h, self.h, -self.h, -self.h], 'f8'), 2.0, dtype='f8'), dtype='f8')
		tempname = 'cutg_dia_' + str(self.run) + '_{c}_{r}_rect'.format(c=self.col_num, r=self.row_num)
		self.cutg = ro.TCutG(tempname, 5, tempx, tempy)
		self.cutg.SetNameTitle(tempname, tempname)
		self.cutg.SetVarX('diaChXPred')
		self.cutg.SetVarY('diaChYPred')
		self.cutg.SetLineColor(ro.kBlack)

if __name__ == '__main__':
	z = RectangularCell3D()