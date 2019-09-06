#!/usr/bin/env python
from ROOT import TFile, TH2F, TH3F, TH1F, TCanvas, TCutG, kRed, gStyle, TBrowser, Long, gROOT
from optparse import OptionParser
from numpy import array, floor
from copy import deepcopy
import os

__author__ = 'DA'

class AnalyseSelectionArea:
    def __init__(self, run=22011, sourceDir='./', outputDir=''):
        print 'Creating AnalyseSelectionArea instance for run:', run
        self.run = run
        self.sourceDir = sourceDir
        self.outputDir = outputDir if outputDir != '' else '{s}/{r}/selectionAnalysis/'.format(r=self.run, s=self.sourceDir)
        self.rootFile = TFile(sourceDir + '/{r}/selectionAnalysis/root/histograms.{r}.{r}.root'.format(r=self.run))
        self.histo3D = TH3F(self.rootFile.Get('hChargeVsFidCut'))
        self.histo3D.GetXaxis().SetTitle('Silicon X/ch')
        self.histo3D.GetYaxis().SetTitle('Silicon Y/ch')
        self.histo1D = 0
        self.map = 0
        self.map_fid = 0
        self.sel_old = {'x_low': self.histo3D.GetXaxis().GetXmin(), 'x_high': self.histo3D.GetXaxis().GetXmax(), 'y_low': self.histo3D.GetYaxis().GetXmin(), 'y_high': self.histo3D.GetYaxis().GetXmax()}
        self.fidcut = 0
        self.fidpoints = []
        self.nameFid = ''
        if not os.path.isdir('{dir}/Plots'.format(dir=self.outputDir)):
            os.makedirs('{dir}/Plots'.format(dir=self.outputDir))
        if not os.path.isdir('{dir}/root'.format(dir=self.outputDir)):
            os.makedirs('{dir}/root'.format(dir=self.outputDir))
        gStyle.SetPalette(55)
        gStyle.SetNumberContours(999)
        self.bla = []

    def Get1DHisto(self, bat=False):
        nbins = self.histo3D.GetNbinsZ()
        nbins_merge = int(floor(nbins / 64))
        print 'Merge {n} bins'.format(n=nbins_merge)
        self.histo3D.RebinZ(nbins_merge)
        if self.fidcut == 0:
            self.CreateFidCut()
        if bat: gROOT.SetBatch(True)
        if len(self.fidpoints) == 5:
            self.histo3D.GetXaxis().SetRangeUser(self.fidpoints[0]['x'], self.fidpoints[2]['x'])
            self.histo3D.GetYaxis().SetRangeUser(self.fidpoints[0]['y'], self.fidpoints[2]['y'])
            self.histo1D = self.histo3D.Project3D('z')
        else:
            self.histo1D = TH1F('hChargeVsFidCut_z', 'hChargeVsFidCut_z', 64, 0, 4096)
            nbins3D = (self.histo3D.GetNbinsY()*self.histo3D.GetNbinsX()*self.histo3D.GetNbinsZ()+1)
            for bin in xrange(1, nbins3D):
                x, y, z = Long(0), Long(0), Long(0)
                if self.histo3D.GetBinContent(bin) >= 1:
                    self.histo3D.GetBinXYZ(bin, x, y, z)
                    point = {'x': self.histo3D.GetXaxis().GetBinCenter(x), 'y': self.histo3D.GetYaxis().GetBinCenter(y)}
                    # if self.WindingIsPointInPoly(point):
                    if bool(int(self.fidcut.IsInside(point['x'], point['y']))):
                        self.histo1D.Fill(self.histo3D.GetZaxis().GetBinCenter(z))
        canvas_name = 'c_s_1D_{r}_{n}'.format(r=self.run, n=self.nameFid)
        canvas = self.CreateCanvas(canvas_name)
        gStyle.SetOptStat('nemr')
        canvas.cd()
        self.histo1D.SetLineWidth(3)
        self.histo1D.GetYaxis().SetTitle('entries')
        self.histo1D.Draw()
        name = 'histo1D_charge_sel_{r}_{n}'.format(r=self.run, n=self.nameFid)
        self.SaveCanvas(canvas, name)
        self.bla.append(canvas)
        if bat: gROOT.SetBatch(False)

    def Get2DMap(self, bat=False):
        if self.fidcut == 0:
            self.CreateFidCut()
        self.histo3D.GetXaxis().SetRangeUser(self.sel_old['x_low'], self.sel_old['x_high'])
        self.histo3D.GetYaxis().SetRangeUser(self.sel_old['y_low'], self.sel_old['y_high'])
        if bat: gROOT.SetBatch(True)
        if self.map == 0:
            self.map = self.histo3D.Project3DProfile('yx')
        self.map.SetTitle(self.map.GetName())
        self.map.GetZaxis().SetTitle('charge/ADC')
        self.map.GetZaxis().SetRangeUser(0,3000)
        canvas_name = 'c_s_2D_{r}_{n}'.format(r=self.run, n=self.nameFid)
        canvas = self.CreateCanvas(canvas_name)
        gStyle.SetOptStat('ne')
        canvas.cd()
        self.map.GetXaxis().SetTitle('Silicon X/ch')
        self.map.GetYaxis().SetTitle('Silicon Y/ch')
        self.map.Draw('colz')
        self.fidcut.Draw('same')
        self.bla.append(canvas)
        name = 'histo2D_charge_sel_{r}_{n}'.format(r=self.run, n=self.nameFid)
        self.SaveCanvas(canvas, name)
        if bat: gROOT.SetBatch(False)

    def Get2DMapFiducial(self):
        if self.fidcut == 0:
            self.CreateFidCut()
        self.histo3D.GetXaxis().SetRangeUser(self.sel_old['x_low'], self.sel_old['x_high'])
        self.histo3D.GetYaxis().SetRangeUser(self.sel_old['y_low'], self.sel_old['y_high'])

        if self.map == 0:
            self.map = self.histo3D.Project3DProfile('yx')
        self.map.GetZaxis().SetTitle('charge/ADC')
        self.map.GetZaxis().SetRangeUser(0, 3000)
        gStyle.SetOptStat('n')
        self.map_fid = self.map.Clone('{n}_FidRegion'.format(n=self.map.GetName()))
        self.map_fid.SetTitle(self.map_fid.GetName())
        nbinsMap = (self.map.GetNbinsY() * self.map.GetNbinsX() + 1)
        for bin in xrange(1, nbinsMap):
            x, y, z = Long(0), Long(0), Long(0)
            if self.map.GetBinContent(bin) >= 0:
                self.map.GetBinXYZ(bin, x, y, z)
                point = {'x': self.map.GetXaxis().GetBinCenter(x), 'y': self.map.GetYaxis().GetBinCenter(y)}
                # if not self.WindingIsPointInPoly(point):
                if not bool(int(self.fidcut.IsInside(point['x'], point['y']))):
                    self.map_fid.SetBinContent(bin, 0)
        canvas_name = 'c_s_2D_{r}_fid_{n}'.format(r=self.run, n=self.nameFid)
        canvas = self.CreateCanvas(canvas_name)
        gStyle.SetOptStat('n')
        canvas.cd()
        self.map_fid.GetXaxis().SetTitle('Silicon X/ch')
        self.map_fid.GetYaxis().SetTitle('Silicon Y/ch')
        self.map_fid.Draw('colz')
        self.fidcut.Draw('same')
        self.bla.append(canvas)
        name = 'histo2D_charge_sel_{r}_fid_{n}'.format(r=self.run, n=self.nameFid)
        self.SaveCanvas(canvas, name)

    def CreateFidCut(self, sel=0):
        if self.fidcut == 0:
            if len(self.fidpoints) > 3:
                self.fidcut = TCutG('fidcut', len(self.fidpoints))
                self.fidcut.SetVarX(self.histo3D.GetXaxis().GetName())
                self.fidcut.SetVarY(self.histo3D.GetYaxis().GetName())
                npoints = 0
                for ipoint in self.fidpoints:
                    self.fidcut.SetPoint(npoints, ipoint['x'], ipoint['y'])
                    npoints += 1
            elif len(self.fidpoints) == 0:
                self.SetFidPoints(0)
            else:
                print 'The fiducial region is not specified correctly. Enter the region again.'
                self.fidpoints = []
                self.SetFidPoints()

        else:
            print 'The fiducial region has already been defined'
            return
        if self.fidcut != 0:
            self.fidcut.SetLineColor(kRed)
            self.fidcut.SetLineWidth(3)

    def SetFidPoints(self, npoint=0):
        if npoint == 0:
            if not bool(int(raw_input('Change default Fiducial Region? (yes: 1 / no: 0'))):
                self.fidpoints.append({'x': self.sel_old['x_low'], 'y': self.sel_old['y_low']})
                self.fidpoints.append({'x': self.sel_old['x_low'], 'y': self.sel_old['y_high']})
                self.fidpoints.append({'x': self.sel_old['x_high'], 'y': self.sel_old['y_high']})
                self.fidpoints.append({'x': self.sel_old['x_high'], 'y': self.sel_old['y_low']})
                self.fidpoints.append(deepcopy(self.fidpoints[0]))
                self.nameFid = ''
                self.CreateFidCut()
                return
            else:
                self.nameFid = raw_input('Enter the name of this region: ')
        point = {'x': 0, 'y': 0}
        stop = True
        while stop:
            point['x'] = float(raw_input('Enter value of X for point {n}: '.format(n=npoint)))
            point['y'] = float(raw_input('Enter value of Y for point {n}: '.format(n=npoint)))
            if self.sel_old['x_low'] <= float(point['x']) <= self.sel_old['x_high'] and self.sel_old['y_low'] <= float(point['y']) <= self.sel_old['y_high']:
                stop = False
            else:
                print 'The specified point is outside of the possible range. Try again...'
        self.fidpoints.append(point)
        if npoint < 2:
            self.SetFidPoints(npoint + 1)
        else:
            stop = True
            while stop:
                cont = raw_input('Insert another point (yes: 1 / no: 0) ? ')
                if cont == '0' or cont == '1':
                    stop = False
                else:
                    print 'Did not entered 1 or 0. Try again...'
            if bool(int(cont)):
                self.SetFidPoints(npoint + 1)
            else:
                self.fidpoints.append(deepcopy(self.fidpoints[0]))
        self.CreateFidCut()

    def SetDefualtFidPoints(self, geom):
        if geom == 4:
            self.nameFid = 'Bla_Square'
            self.fidpoints.append({'x': 95.2, 'y': 65.3})
            self.fidpoints.append({'x': 95.2, 'y': 67.3})
            self.fidpoints.append({'x': 97.2, 'y': 67.3})
            self.fidpoints.append({'x': 97.2, 'y': 65.3})
            self.fidpoints.append(deepcopy(self.fidpoints[0]))
            self.CreateFidCut()
        elif geom == 6:
            self.nameFid = 'Bla_Hexa'
            self.fidpoints.append({'x': 86.8, 'y': 101.8})
            self.fidpoints.append({'x': 85.5, 'y': 101.8})
            self.fidpoints.append({'x': 84.8, 'y': 102.9})
            self.fidpoints.append({'x': 85.5, 'y': 104.1})
            self.fidpoints.append({'x': 86.8, 'y': 104.1})
            self.fidpoints.append({'x': 87.5, 'y': 102.9})
            self.fidpoints.append(deepcopy(self.fidpoints[0]))
            self.CreateFidCut()

    def SetDefaultFidPointsStrip(self):
        self.nameFid = 'Bla_Strip'
        self.fidpoints.append({'x': 127.6, 'y': 90.5})
        self.fidpoints.append({'x': 127.6, 'y': 108.5})
        self.fidpoints.append({'x': 132.6, 'y': 108.5})
        self.fidpoints.append({'x': 132.6, 'y': 90.5})
        self.fidpoints.append(deepcopy(self.fidpoints[0]))
        self.CreateFidCut()
        
    def SetDefaultFidPointsSingleC(self):
        self.nameFid = 'Bla_SingleC'
        self.fidpoints.append({'x': 110, 'y': 65})
        self.fidpoints.append({'x': 110, 'y': 95})
        self.fidpoints.append({'x': 130, 'y': 95})
        self.fidpoints.append({'x': 130, 'y': 65})
        self.fidpoints.append(deepcopy(self.fidpoints[0]))
        self.CreateFidCut()
        
    def SetDefaultFidPointsFull3D(self):
        self.nameFid = 'Bla_Full3D'
        self.fidpoints.append({'x': 165.66, 'y': 42})
        self.fidpoints.append({'x': 165.66, 'y': 60})
        self.fidpoints.append({'x': 189.66, 'y': 60})
        self.fidpoints.append({'x': 189.66, 'y': 42})
        self.fidpoints.append(deepcopy(self.fidpoints[0]))
        self.CreateFidCut()

    def SetDefaultFidPointsBigHex(self):
        self.nameFid = 'Bla_BigHex'
        self.fidpoints.append({'x': 82, 'y': 98})
        self.fidpoints.append({'x': 82, 'y': 103})
        self.fidpoints.append({'x': 99, 'y': 103})
        self.fidpoints.append({'x': 99, 'y': 98})
        self.fidpoints.append(deepcopy(self.fidpoints[0]))
        self.CreateFidCut()

    def SetDefaultFidPointsBigSquare(self):
        self.nameFid = 'Bla_BigSquare'
        self.fidpoints.append({'x': 85, 'y': 63})
        self.fidpoints.append({'x': 85, 'y': 67})
        self.fidpoints.append({'x': 99, 'y': 67})
        self.fidpoints.append({'x': 99, 'y': 63})
        self.fidpoints.append(deepcopy(self.fidpoints[0]))
        self.CreateFidCut()

    def CreateCanvas(self, name='c0', w=1024, h=1024):
        c = TCanvas(name, name, w, h)
        c.SetWindowSize(w + (w - c.GetWw()), h + (h - c.GetWh()))
        return deepcopy(c)

    def SaveCanvas(self, canvas, name):
        canvas.SaveAs(self.outputDir + '/root/{n}.root'.format(n=name))
        canvas.SaveAs(self.outputDir + '/Plots/{n}.png'.format(n=name))

    def GetHistograms(self, bat=False):
        self.Get2DMap(bat)
        self.Get1DHisto(bat)

    # ------------------------------------------------------------------------------------------------------------------
    # From http://geomalgorithms.com/a03-_inclusion.html#wn_PnPoly()
    def IsLeft(self,p0, p1, p2):
        return ((p1['x']-p0['x'])*(p2['y']-p0['y'])-(p2['x']-p0['x'])*(p1['y']-p0['y']))

    def WindingIsPointInPoly(self, point):
        wn = 0
        n = len(self.fidpoints) - 1
        for i in xrange(n):
            if self.fidpoints[i]['y'] <= point['y']:
                if self.fidpoints[i+1]['y'] > point['y']:
                    if self.IsLeft(self.fidpoints[i], self.fidpoints[i+1], point) > 0:
                        wn += 1
            elif self.fidpoints[i+1]['y'] <= point['y']:
                if self.IsLeft(self.fidpoints[i], self.fidpoints[i+1], point) < 0:
                    wn -= 1
        return bool(int(wn))
    # ------------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option('-r', '--run', dest='run', default=22011, type='int', help='Run to be analysed (e.g. 22011)')
    parser.add_option('-s', '--sourceDir', dest='source', default='./', type='string', help='source folder containing processed data of different runs')
    parser.add_option('-o', '--outputDir', dest='output', default='', type='string', help='output folder containing the analysed results')

    (options, args) = parser.parse_args()
    run = int(options.run)
    source = str(options.source)
    output = str(options.output)

    z = AnalyseSelectionArea(run, source, output)
