#!/usr/bin/env python

import ROOT as ro
import numpy as np
import shapely.geometry as geom
import shapely.ops as ops
import ipdb
import sys
from Utils import *

class GridAreas:
	def __init__(self, numcols=0, numrows=0, run=0):
		self.num_cols = numcols
		self.num_rows = numrows
		self.run = run
		self.goodAreas_diamond = []
		self.goodAreas_simplified_diamond = []
		self.goodAreas_diamond_centers = []
		self.goodAreas_index = []
		self.badAreas_diamond = []
		self.badAreas_simplified_diamond = []
		self.badAreas_diamond_centers = []
		self.badAreas_index = []
		self.goodAreasCutNames_diamond = ''
		self.goodAreasCutNames_simplified_diamond = ''
		self.goodAreasCutNames_diamond_centers = ''
		self.badAreasCutNames_diamond = ''
		self.badAreasCutNames_simplified_diamond = ''
		self.badAreasCutNames_diamond_centers = ''

		self.center_rectangles = {}

	def CreateRectangleSymmetricCentralRegion(self, percentage=80, pitchCol=50, pitchRow=50, xvar='(((diaChXPred-0.5)*(50*10000))%(50*10000))/10000', yvar='(((diaChYPred-25)*10000)%(50*10000))/10000'):
		xinf = np.divide(pitchCol * (1 - np.sqrt(percentage / 100., dtype='float128')), 2.0, dtype='float128')
		xsup = np.subtract(pitchCol, xinf, dtype='float128')
		yinf = np.divide(pitchRow * (1 - np.sqrt(percentage / 100., dtype='float128')), 2.0, dtype='float128')
		ysup = np.subtract(pitchRow, yinf, dtype='float128')
		xarray = np.array([xinf, xinf, xsup, xsup, xinf], dtype='float64')
		yarray = np.array([yinf, ysup, ysup, yinf, yinf], dtype='float64')
		tempname = 'cutg_center_rect_{r}_{p}'.format(r=self.run, p=percentage)
		tempCutG = ro.TCutG(tempname, 5, xarray, yarray)
		tempCutG.SetNameTitle(tempname, tempname)
		tempCutG.SetVarX(xvar)
		tempCutG.SetVarY(yvar)
		tempCutG.SetLineColor(ro.kRed)
		self.center_rectangles[percentage] = {'name': tempname, 'tcutg': tempCutG}

	def AddGoodAreas(self, col, row, tcutgs_diamond, tcutgs_diamond_center):
		if 0 <= col < self.num_cols and 0 <= row < self.num_rows:
			self.goodAreas_index.append((col, row))
			tcutgs_diamond[col][row].SetLineColor(ro.kRed)
			tcutgs_diamond_center[col][row].SetLineColor(ro.kViolet)
			self.goodAreas_diamond.append(tcutgs_diamond[col][row])
			self.goodAreas_diamond_centers.append(tcutgs_diamond_center[col][row])
			tempgood = [cut.GetName() for cut in self.goodAreas_diamond]
			self.goodAreasCutNames_diamond = '((' + ')||('.join(tempgood) + '))'
			tempgood = [cut.GetName() for cut in self.goodAreas_diamond_centers]
			self.goodAreasCutNames_diamond_centers = '((' + ')||('.join(tempgood) + '))'

	def AddBadAreas(self, col, row, tcutgs_diamond, tcutgs_diamond_center):
		if 0 <= col < self.num_cols and 0 <= row < self.num_rows:
			self.badAreas_index.append((col, row))
			tcutgs_diamond[col][row].SetLineColor(ro.kBlue)
			tcutgs_diamond_center[col][row].SetLineColor(ro.kViolet)
			self.badAreas_diamond.append(tcutgs_diamond[col][row])
			self.badAreas_diamond_centers.append(tcutgs_diamond_center[col][row])
			tempbad = [cut.GetName() for cut in self.badAreas_diamond]
			self.badAreasCutNames_diamond = '((' + ')||('.join(tempbad) + '))'
			tempbad = [cut.GetName() for cut in self.badAreas_diamond_centers]
			self.badAreasCutNames_diamond_centers = '((' + ')||('.join(tempbad) + '))'

	def AddGoodAreasRow(self, row, coli=0, colf=0, tcutgs_diamond=None, tcutgs_diamond_center=None):
		(colii, colff) = (0, self.num_cols) if coli == 0 and colf == 0 else (coli, colf)
		if 0 <= colii <= colff < self.num_cols and 0 <= row < self.num_rows:
			for col in xrange(colii, colff + 1):
				self.goodAreas_index.append((col, row))
				tcutgs_diamond[col][row].SetLineColor(ro.kRed)
				tcutgs_diamond_center[col][row].SetLineColor(ro.kViolet)
				self.goodAreas_diamond.append(tcutgs_diamond[col][row])
				self.goodAreas_diamond_centers.append(tcutgs_diamond_center[col][row])

			tempgood = [cut.GetName() for cut in self.goodAreas_diamond]
			self.goodAreasCutNames_diamond = '((' + ')||('.join(tempgood) + '))'
			tempgood = [cut.GetName() for cut in self.goodAreas_diamond_centers]
			self.goodAreasCutNames_diamond_centers = '((' + ')||('.join(tempgood) + '))'

	def AddGoodAreasCol(self, col, rowi=0, rowf=0, tcutgs_diamond=None, tcutgs_diamond_center=None):
		(rowii, rowff) = (0, self.num_rows) if rowi == 0 and rowf == 0 else (rowi, rowf)
		if 0 <= col < self.num_cols and 0 <= rowii <= rowff < self.num_rows:
			for row in xrange(rowii, rowff + 1):
				self.goodAreas_index.append((col, row))
				tcutgs_diamond[col][row].SetLineColor(ro.kRed)
				tcutgs_diamond_center[col][row].SetLineColor(ro.kViolet)
				self.goodAreas_diamond.append(tcutgs_diamond[col][row])
				self.goodAreas_diamond_centers.append(tcutgs_diamond_center[col][row])

			tempgood = [cut.GetName() for cut in self.goodAreas_diamond]
			self.goodAreasCutNames_diamond = '((' + ')||('.join(tempgood) + '))'
			tempgood = [cut.GetName() for cut in self.goodAreas_diamond_centers]
			self.goodAreasCutNames_diamond_centers = '((' + ')||('.join(tempgood) + '))'

	def AddRemainingToBadAreas(self, tcutgs_diamond=None, tcutgs_diamond_center=None):
		for col in xrange(0, self.num_cols):
			for row in xrange(0, self.num_rows):
				if (col, row) not in self.goodAreas_index and (col, row) not in self.badAreas_index:
					self.AddBadAreas(col, row, tcutgs_diamond, tcutgs_diamond_center)

	def RemoveFromGoodArea(self, col, row, tcutgs_diamond=None, tcutgs_diamond_center=None):
		if (col, row) in self.goodAreas_index:
			index_g = self.goodAreas_index.index((col, row))
			self.goodAreas_diamond.pop(index_g)
			self.goodAreas_diamond_centers.pop(index_g)
			self.goodAreas_index.pop(index_g)
			self.AddBadAreas(col, row, tcutgs_diamond, tcutgs_diamond_center)
			tempgood = [cut.GetName() for cut in self.goodAreas_diamond]
			self.goodAreasCutNames_diamond = '((' + ')||('.join(tempgood) + '))'
			tempgood = [cut.GetName() for cut in self.goodAreas_diamond_centers]
			self.goodAreasCutNames_diamond_centers = '((' + ')||('.join(tempgood) + '))'

	def SimplifyGoodAndBadAreas(self):
		def CreateTCutGDic(polygon_dic, color0=ro.kRed, color1=ro.kBlue):
			polygons = polygon_dic['polygons']
			names = polygon_dic['names']
			list_simplified = []
			for i, polyg in enumerate(polygons):
				xi, yi = polyg.exterior.xy
				temptcg = ro.TCutG(names[i], len(xi), xi, yi)
				temptcg.SetNameTitle(names[i], names[i])
				temptcg.SetVarX(self.goodAreas_diamond[0].GetVarX())
				temptcg.SetVarY(self.goodAreas_diamond[0].GetVarY())
				temptcg.SetLineColor(color0)
				list_holes = []
				for hi in xrange(len(polyg.interiors)):
					xhi, yhi = polyg.interiors[hi].xy
					tempthi = ro.TCutG(names[i] + '_h' + str(hi), len(xhi), xhi, yhi)
					tempthi.SetNameTitle(names[i] + '_h' + str(hi), names[i] + '_h' + str(hi))
					tempthi.SetVarX(self.goodAreas_diamond[0].GetVarX())
					tempthi.SetVarY(self.goodAreas_diamond[0].GetVarY())
					tempthi.SetLineColor(color1)
					list_holes.append(tempthi)
				list_simplified.append({'polygon': temptcg, 'holes': list_holes})
			return list_simplified

		def ReturnNameCutGWithHoles(polygon_dic_list):
			list_names = []
			for polygon_dic in polygon_dic_list:
				holes = polygon_dic['holes']
				polygon = polygon_dic['polygon']
				not_hole_names = ['(!' + hole.GetName() + ')' for hole in holes]
				not_holes_part = '(' + '&&'.join(not_hole_names) + ')'
				polygon_part_with_holes = '(' + polygon.GetName() + '&&' + not_holes_part + ')' if len(holes) != 0 else '(' + polygon.GetName() + ')'
				list_names.append(polygon_part_with_holes)
			return '(' + '||'.join(list_names) + ')'

		good_polygons_dic = self.SimplifyAreas(self.goodAreas_diamond, [], [], prefix='g')
		self.goodAreas_simplified_diamond = CreateTCutGDic(good_polygons_dic, ro.kRed, ro.kBlue)
		self.goodAreasCutNames_simplified_diamond = ReturnNameCutGWithHoles(self.goodAreas_simplified_diamond)

		bad_polygons_dic = self.SimplifyAreas(self.badAreas_diamond, [], [], prefix='b')
		self.badAreas_simplified_diamond = CreateTCutGDic(bad_polygons_dic, ro.kBlue, ro.kRed)
		self.badAreasCutNames_simplified_diamond = ReturnNameCutGWithHoles(self.badAreas_simplified_diamond)

	def SimplifyAreas(self, area_list, polygon_list=[], name_list=[], prefix='g', num_merged=0, finished=False):
		def CheckAreas(ai, aj):
			def DistancePointLine(xa, ya, xb, yb, x0, y0):
				return np.divide(np.abs((yb - ya) * x0 - (xb - xa) * y0 + xb * ya - yb * xa, dtype='float64'), np.sqrt(np.power(yb - ya, 2) + np.power(xb - xa, 2), dtype='float64'), dtype='float64')

			def CheckPointBetweenPoints(xa, ya, xb, yb, x0, y0):
				d1 = np.sqrt(np.power(x0 - xa, 2.0, dtype='f8') + np.power(y0 - ya, 2.0, dtype='f8'), dtype='f8')
				d2 = np.sqrt(np.power(x0 - xb, 2.0, dtype='f8') + np.power(y0 - yb, 2.0, dtype='f8'), dtype='f8')
				d3 = np.sqrt(np.power(xb - xa, 2.0, dtype='f8') + np.power(yb - ya, 2.0, dtype='f8'), dtype='f8')
				# return True if d3 == d1 + d2 else False
				return np.abs(d1 + d2 - d3) < 1e-12
			nai, naj = ai.GetN(), aj.GetN()
			xai, xaj, yai, yaj = ai.GetX(), aj.GetX(), ai.GetY(), aj.GetY()
			xai, xaj, yai, yaj = [xai[i] for i in xrange(nai)], [xaj[i] for i in xrange(naj)], [yai[i] for i in xrange(nai)], [yaj[i] for i in xrange(naj)]
			for i in xrange(nai - 1):
				for j in xrange(naj - 1):
					# if DistancePointLine(xai[i], yai[i], xai[i+1], yai[i+1], xaj[j], yaj[j]) + DistancePointLine(xai[i], yai[i], xai[i+1], yai[i+1], xaj[j+1], yaj[j+1]) == 0:
					if CheckPointBetweenPoints(xai[i], yai[i], xai[i+1], yai[i+1], xaj[j], yaj[j]) and CheckPointBetweenPoints(xai[i], yai[i], xai[i+1], yai[i+1], xaj[j+1], yaj[j+1]):
						return i, i+1, j, j+1
			return 0, 0, 0, 0

		def MergeAreas(ai, aj, prefix='g', num=0):
			# nai, naj = ai.GetN(), aj.GetN()
			# xai, xaj, yai, yaj = ai.GetX(), aj.GetX(), ai.GetY(), aj.GetY()
			# xai, xaj, yai, yaj = [xai[i] for i in xrange(nai)], [xaj[i] for i in xrange(naj)], [yai[i] for i in xrange(nai)], [yaj[i] for i in xrange(naj)]
			# poli = geom.Polygon([[xai[i], yai[i]] for i in xrange(nai)])
			# polj = geom.Polygon([[xaj[i], yaj[i]] for i in xrange(naj)])
			# polij = ops.cascaded_union([poli, polj])
			polij = ops.cascaded_union([ai, aj])
			polij = polij.simplify(0.1)
			if polij.type == 'Polygon':
			# try:
				xnew, ynew = polij.exterior.xy
				# colrowai = ai.GetName().split('cutg_dia_' + str(self.run) + '_')[-1].split('_')
				# colrowai = namei.split('cutg_dia_' + str(self.run) + '_')[-1].split('_')
				# colrowaj = aj.GetName().split('cutg_dia_' + str(self.run) + '_')[-1].split('_')
				# colrowaj = namej.split('cutg_dia_' + str(self.run) + '_')[-1].split('_')
				# colrownew = colrowai + colrowaj
				newname = 'cutg_dia_' + str(self.run) + '_' + prefix + '_' + str(num)
				# newname = 'cutg_dia_' + str(self.run) + '_' + '_'.join(colrownew)
				# newtcutg = ro.TCutG(newname, len(xnew), xnew, ynew)
				# newtcutg.SetNameTitle(newname, newname)
				# newtcutg.SetVarX(ai.GetVarX())
				# newtcutg.SetVarY(ai.GetVarY())
				# newtcutg.SetLineColor(ro.kRed)
				dic_newpolyg = {'name': newname, 'polygon': polij}
			# except Exception:
			else:
				# ipdb.set_trace()
				dic_newpolyg = None
			return dic_newpolyg

		def CreateNewList(area_listij, name_list, posi, posj, dic_areaij):
			new_area_list = [area_listij[i] for i in xrange(len(area_listij)) if i not in [posi, posj]]
			new_name_list = [name_list[i] for i in xrange(len(name_list)) if i not in [posi, posj]]
			new_area_list.insert(0, dic_areaij['polygon'])
			new_name_list.insert(0, dic_areaij['name'])
			return {'names': new_name_list, 'polygons': new_area_list}

		if finished:
			# self.goodAreas_simplified_diamond = area_list
			return {'polygons': polygon_list, 'names': name_list}
		if len(polygon_list) == 0:
			lentemp0 = len(area_list)
			maxperm = int(lentemp0 * (lentemp0 - 1))
			# maxperm = int(lentemp0 * (lentemp0 - 1) / 2.0)
			sys.setrecursionlimit(maxperm)
			for areai in area_list:
				ni = areai.GetN()
				xi, yi = areai.GetX(), areai.GetY()
				xi, yi = [xi[i] for i in xrange(ni)], [yi[i] for i in xrange(ni)]
				poli = geom.Polygon([[xi[i], yi[i]] for i in xrange(ni)])
				polygon_list.append(poli)
				namei = areai.GetName()
				namei = 'cutg_dia_' + str(self.run) + '_' + prefix + '_' + '_'.join(namei.split('_')[-2:])
				name_list.append(namei)
		finito = True
		new_list = polygon_list[:]
		new_names = name_list[:]
		exitflag = False
		temp0 = polygon_list[:]
		tempnames0 = name_list[:]
		for i in xrange(len(temp0)):
			for j in xrange(i+1, len(temp0)):
				ai, aj = temp0[i], temp0[j]
				# i0, i1, j0, j1 = CheckAreas(ai, aj)
				# if i0 + i1 + j0 + j1 != 0:
				dic_areaij = MergeAreas(ai, aj, prefix, num_merged)
				if dic_areaij:
					num_merged += 1
					dic_new_list = CreateNewList(temp0, tempnames0, i, j, dic_areaij)
					new_list = dic_new_list['polygons']
					new_names = dic_new_list['names']
					finito = False
					exitflag = True
					break
			if exitflag:
				break
		# ipdb.set_trace()
		# if finito:
		# 	ipdb.set_trace()
		return self.SimplifyAreas(area_list, new_list, new_names, prefix, num_merged, finito)

	def ResetAreas(self):
		self.goodAreas_diamond = []
		self.goodAreas_diamond_centers = []
		self.badAreas_diamond = []
		self.badAreas_diamond_centers = []
		self.goodAreas_index = []
		self.badAreas_index = []
		self.goodAreasCutNames_diamond = ''
		self.goodAreasCutNames_diamond_centers = ''
		self.badAreasCutNames_diamond = ''
		self.badAreasCutNames_diamond_centers = ''

if __name__ == '__main__':
	ga = GridAreas(0, 0)
