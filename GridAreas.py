#!/usr/bin/env python

import ROOT as ro
import numpy as np
import shapely.geometry as geom
import shapely.ops as ops
import ipdb
import sys
from Utils import *

class GridAreas:
	def __init__(self, numcols=0, numrows_even=0, numrows_odd=0, run=0):
		self.num_cols = numcols
		self.num_rows_even = numrows_even
		self.num_rows_odd = numrows_odd
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

		self.center_polygons = {}

		self.good_channels = []
		self.bad_channels = []

	def CreateProportionalCellCentralRegion(self, percentage, diacols, width, varx, vary):
		"""
		Creates a dictionary with a polygon scaled by a percentage with the shape of the cell for inner cell studies. The polygon as a TCutG and its name are added in the center_polygons dictionary as another dictionary
		:param percentage: percentage (0<=percentage<=100) to scale the polygon of the cell
		:param diacols: DiamondColumns object that has the information of the columns and the cells in the rows of the transparent grid
		:return:
		"""
		fraction = np.sqrt(np.divide(percentage, 100., dtype='f8'), dtype='f8') if percentage != 0 else np.sqrt(np.divide(0.0001, 100., dtype='f8'), dtype='f8')
		xarray = np.multiply(diacols.cols[0].cells[0].GetXCoordinatesPolygon(xcenter=0, fraction=fraction), width, dtype='f8')
		yarray = diacols.cols[0].cells[0].GetYCoordinatesPolygon(ycenter=0, fraction=fraction)

		tempname = 'cutg_center_polygon_{r}_{p}'.format(r=self.run, p=percentage)
		tempCutG = ro.TCutG(tempname, xarray.size, xarray, yarray)
		tempCutG.SetNameTitle(tempname, tempname)
		tempCutG.SetVarX(varx)
		tempCutG.SetVarY(vary)
		tempCutG.SetLineColor(ro.kBlack)
		tempCutG.SetLineWidth(2)
		self.center_polygons[percentage] = {'name': tempname, 'tcutg': tempCutG}

	def AddGoodAreas(self, col, row, dia_cols):
		"""
		This method tags a cell as 'good' (selected). It updates the list goodAreas_diamond, goodAreas_diamond_centers, goodAreas_index and updates the strings goodAreasCutNames_diamond, goodAreasCutNames_diamond_centers
		:param col: is the column number of the cell to be added
		:param row: is the row number of the cell to be added
		:param dia_cols: is the DiamondColumns object with the information of the columns in the transparent grid
		:return:
		"""
		if 0 <= col < self.num_cols:
			if (col % 2 == 0) and (0 <= row < self.num_rows_even) or (col % 2 != 0) and (0 <= row < self.num_rows_odd):
				self.goodAreas_index.append((col, row))
				dia_cols.cols[col].cells[row].cutg.SetLineColor(ro.kRed)
				self.goodAreas_diamond.append(dia_cols.cols[col].cells[row].cutg)
				self.goodAreas_diamond_centers.append(dia_cols.cols[col].cells[row].cutg_center)
				tempgood = [cut.GetName() for cut in self.goodAreas_diamond]
				self.goodAreasCutNames_diamond = '((' + ')||('.join(tempgood) + '))'
				tempgood = [cut.GetName() for cut in self.goodAreas_diamond_centers]
				self.goodAreasCutNames_diamond_centers = '((' + ')||('.join(tempgood) + '))'

	def AddBadAreas(self, col, row, dia_cols):
		"""
		This method tags a cell as 'bad' (not selected). It updates the list badAreas_diamond, badAreas_diamond_centers, badAreas_index and updates the strings badAreasCutNames_diamond, badAreasCutNames_diamond_centers
		:param col: is the column number of the cell to be added
		:param row: is the row number of the cell to be added
		:param dia_cols: is the DiamondColumns object with the information of the columns in the transparent grid
		:return:
		"""
		if 0 <= col < self.num_cols:
			if (col % 2 == 0) and (0 <= row < self.num_rows_even) or (col % 2 != 0) and (0 <= row < self.num_rows_odd):
				self.badAreas_index.append((col, row))
				dia_cols.cols[col].cells[row].cutg.SetLineColor(ro.kBlack)
				self.badAreas_diamond.append(dia_cols.cols[col].cells[row].cutg)
				self.badAreas_diamond_centers.append(dia_cols.cols[col].cells[row].cutg_center)
				tempbad = [cut.GetName() for cut in self.badAreas_diamond]
				self.badAreasCutNames_diamond = '((' + ')||('.join(tempbad) + '))'
				tempbad = [cut.GetName() for cut in self.badAreas_diamond_centers]
				self.badAreasCutNames_diamond_centers = '((' + ')||('.join(tempbad) + '))'

	def AddGoodAreasRow(self, row, coli=0, colf=0, dia_cols=None):
		"""
		Adds the row 'row' starting from column 'coli' until 'colf' as good (in selection)
		:param row: number of the row
		:param coli: starting column
		:param colf: finishing column
		:param dia_cols: DiamondColumns object with the information of the columns in the transparent grid
		:return:
		"""
		(colii, colff) = (0, self.num_cols) if coli == 0 and colf == 0 else (coli, colf)
		if 0 <= colii <= colff < self.num_cols and 0 <= row < min(self.num_rows_even, self.num_rows_odd):
			for col in xrange(colii, colff + 1):
				self.AddGoodAreas(col, row, dia_cols)

	def AddGoodAreasCol(self, col, rowi=0, rowf=0, dia_cols=None):
		"""
		Adds the column 'col' starting from row 'rowi' until 'rowf' as good (in selection)
		:param col: number of the row
		:param rowi: starting column
		:param rowf: finishing column
		:param dia_cols: DiamondColumns object with the information of the columns in the transparent grid
		:return:
		"""
		if col % 2 == 0:
			(rowii, rowff) = (0, self.num_rows_even) if rowi == 0 and rowf == 0 else (rowi, rowf)
		else:
			(rowii, rowff) = (0, self.num_rows_odd) if rowi == 0 and rowf == 0 else (rowi, rowf)
		if 0 <= col < self.num_cols and ((self.num_rows_even > rowff >= rowii >= 0 == col % 2) or (self.num_rows_odd > rowff >= rowii >= 0 != col % 2)):
			for row in xrange(rowii, rowff + 1):
				self.AddGoodAreas(col, row, dia_cols)

	def AddRemainingToBadAreas(self, dia_cols=None):
		"""
		Tags the remaining cells not in goodAreas_index as bad (not in selection)
		:param dia_cols: is the DiamondColumns object that has the information of the columns in the transparent grid
		:return:
		"""
		if dia_cols:
			for col in xrange(dia_cols.num_cols):
				for row in xrange(dia_cols.cols[col].num_rows):
					if (col, row) not in self.goodAreas_index and (col, row) not in self.badAreas_index:
						self.AddBadAreas(col, row, dia_cols)

	def RemoveFromGoodArea(self, col, row, dia_cols=None):
		"""
		If a cell identified by column 'col' and row 'row' was tagged as good (in selection), it is untagged. The cell is then tagged as bad (not in selection)
		:param col: cell's column number
		:param row: cell's row number
		:param dia_cols: DiamondColumns object that has the information of the columns in the transparent grid
		:return:
		"""
		if (col, row) in self.goodAreas_index:
			index_g = self.goodAreas_index.index((col, row))
			self.goodAreas_diamond.pop(index_g)
			self.goodAreas_diamond_centers.pop(index_g)
			self.goodAreas_index.pop(index_g)
			tempgood = [cut.GetName() for cut in self.goodAreas_diamond]
			self.goodAreasCutNames_diamond = '((' + ')||('.join(tempgood) + '))'
			tempgood = [cut.GetName() for cut in self.goodAreas_diamond_centers]
			self.goodAreasCutNames_diamond_centers = '((' + ')||('.join(tempgood) + '))'
		if (col, row) not in self.badAreas_index:
			self.AddBadAreas(col, row, dia_cols)

	def SimplifyGoodAndBadAreas(self):
		"""
		This method simplifies connected cells into a bit TCutG with TCutG holes if necessary to reduce the computational time if the number of cells is big. This is done by creating polygons and merging them if they are connected
		The result is stored in {good,bad}Areas_simplified_diamond which contains a list with dictionaries, each of which contains a 'polygon' with the TCutG of the polygon, and 'holes' with a list of TCutGs of the holes inside the polygon.
		It also stores the {good,bad}AreasCutNames_simplified_diamond with the cut string used to select these areas in the analysis
		:return:
		"""
		def CreateTCutGDic(polygon_dic, color0=ro.kRed, color1=ro.kBlack):
			"""
			This method creates new TCutGs for the polygon_dic
			:param polygon_dic: Its a dictionary that contains the information of the polygons under 'polygons' and the names of the polygons under 'names'
			:param color0: color used for the polygons
			:param color1: color used for the holes inside each of the polygons
			:return: returns a list of dictionaries. Each dictionary contains the TCutG of the polygon under 'polygon' and a list of TCutGs of the holes inside the polygon under the key 'holes'
			"""
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
			"""
			This method creates the string of the cut that will be used to select the polygons and exclude the holes inside the polygons
			:param polygon_dic_list: is a list of dictionaries. Each dictionary has the TCutG of the polygon and a list of TCutGs of the holes inside the polygon
			:return: returns a string with of the form '((polygon0_name&&((!polygon0_hole0_name)&&(!polygon0_hole1_name)&&...&&(!polygon0_holeM0_name)))||(polygon1_name&&((!polygon1_hole0_name)&&...&&(!polygon1_holeM1_name)))||...||(polygonN_name&&((!polygonN_hole0_name)&&...&&(!polygonN_holeMN_name))))'
			"""
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
		self.goodAreas_simplified_diamond = CreateTCutGDic(good_polygons_dic, ro.kRed, ro.kBlack)
		self.goodAreasCutNames_simplified_diamond = ReturnNameCutGWithHoles(self.goodAreas_simplified_diamond)

		bad_polygons_dic = self.SimplifyAreas(self.badAreas_diamond, [], [], prefix='b')
		self.badAreas_simplified_diamond = CreateTCutGDic(bad_polygons_dic, ro.kBlack, ro.kRed)
		self.badAreasCutNames_simplified_diamond = ReturnNameCutGWithHoles(self.badAreas_simplified_diamond)

		self.GetGoodAndBadChannels()

	def SimplifyAreas(self, area_list, polygon_list=[], name_list=[], prefix='g', num_merged=0, finished=False):
		"""
		This Method simplifies a list of areas if they are connected
		:param area_list: is a list of TCutGs
		:param polygon_list: is a list of polygons
		:param name_list: is a list of names
		:param prefix: is either 'g' for good areas (selected) or 'b' for bad areas (not selected)
		:param num_merged: its the identifier of the number of merged polygons. This is used to have unique ids on the merged polygons
		:param finished: parameter that ends the recursion to return the result
		:return: returns a dictionary with a polygon list under 'polygons' and a name list under 'names'
		"""
		def MergeAreas(ai, aj, prefix='g', num=0):
			"""
			Merges areas ai, aj with the prefix 'g' or 'b' for 'good' or 'bad' tags respectively
			:param ai: TCutG ai to be merged
			:param aj: TCutG aj to be merged
			:param prefix: 'g' or 'b' which corresponds to 'good' or 'bad' tags respectively
			:param num: identifier of merged area
			:return: If merging is successful, it returns a dictionary with the new name of the merged polygon under 'name' and the merged polygon under 'polygon'. If the merge fails, it returns None
			"""
			polij = ops.cascaded_union([ai, aj])
			polij = polij.simplify(0.1)
			if polij.type == 'Polygon':
				newname = 'cutg_dia_' + str(self.run) + '_' + prefix + '_' + str(num)
				dic_newpolyg = {'name': newname, 'polygon': polij}
			else:
				dic_newpolyg = None
			return dic_newpolyg

		def CreateNewList(area_listij, name_list, posi, posj, dic_areaij):
			"""
			Creates a new list of with the merged polygon by popping out the polygons in positions 'posi' and 'posj' and inserting the new merged polygon in the polygon list and in the names list.
			:param area_listij: is the polygon list to be modified with the new merged polygon
			:param name_list: is the name list to be modified with the name of the new merged polygon
			:param posi: position posi of the first polygon used for the merging
			:param posj: position posj of the second polygon used for the merging
			:param dic_areaij: dictrionary with the name and polygon of the new merged polygon
			:return: returns a dictionary with the new modified name list 'names' and the new modified polygon list 'polygons'
			"""
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
			maxperm = max(int(lentemp0 * (lentemp0 - 1)), 1)  # maximum number of permutations
			sys.setrecursionlimit(maxperm)
			for areai in area_list:
				ni = areai.GetN()
				xi, yi = areai.GetX(), areai.GetY()
				xi, yi = [xi[i] for i in xrange(ni)], [yi[i] for i in xrange(ni)]
				poli = geom.Polygon([[xi[i], yi[i]] for i in xrange(ni)])
				polygon_list.append(poli)
				namei = areai.GetName()
				namei = 'cutg_dia_' + str(self.run) + '_' + prefix + '_' + '_'.join(namei.split('_')[-3:])
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
		return self.SimplifyAreas(area_list, new_list, new_names, prefix, num_merged, finito)

	def GetGoodAndBadChannels(self):
		#  starting point should be the leftmost of the bottomest point of the region!
		self.good_channels = list({RoundInt((area.GetX()[1] + area.GetX()[area.GetN() - 3]) / 2.0) for area in self.goodAreas_diamond})
		self.bad_channels = list({RoundInt((area.GetX()[1] + area.GetX()[area.GetN() - 3]) / 2.0) for area in self.badAreas_diamond})

	def ResetAreas(self):
		"""
		Resets all lists and strings
		:return:
		"""
		self.goodAreas_diamond = []
		self.goodAreas_diamond_centers = []
		self.badAreas_diamond = []
		self.badAreas_diamond_centers = []
		self.goodAreas_index = []
		self.badAreas_index = []
		self.goodAreasCutNames_diamond = ''
		self.goodAreasCutNames_simplified_diamond = ''
		self.goodAreasCutNames_diamond_centers = ''
		self.badAreasCutNames_diamond = ''
		self.badAreasCutNames_simplified_diamond = ''
		self.badAreasCutNames_diamond_centers = ''

if __name__ == '__main__':
	ga = GridAreas(0, 0)
