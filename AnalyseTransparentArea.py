#!/usr/bin/env python
from ROOT import TFile, TH2F, TH3F, TH1F, TCanvas, TCutG, kRed, gStyle, TBrowser, Long, TProfile2D
from optparse import OptionParser
from numpy import array, floor
from copy import deepcopy
import os

__author__ = 'DA'

class AnalyseTransparentArea:
    def __init__(self, run=22011, sourceDir='./', outputDir=''):
        print 'Creating AnalyseTransparentArea instance for run:', run
        self.run = run
        self.sourceDir = sourceDir
        self.outputDir = outputDir if outputDir != '' else '{s}/{r}/transparentAnalysis/'.format(r=self.run, s=self.sourceDir)
        self.rootFile1 = TFile(sourceDir + '/{r}/transparentAnalysis/root/hLandau1HighestHitProfile_1OutOf10.{r}.root'.format(r=self.run))
        self.rootFile2 = TFile(sourceDir + '/{r}/transparentAnalysis/root/hLandau2HighestHitProfile_2OutOf10.{r}.root'.format(r=self.run))
        self.histo2D_1 = TProfile2D(self.rootFile1.Get('cRoot_hLandau1HighestHitProfile_1OutOf10').GetPrimitive('hLandau1HighestHitProfile_1OutOf10'))
        self.histo2D_2 = TProfile2D(self.rootFile2.Get('cRoot_hLandau2HighestHitProfile_2OutOf10').GetPrimitive('hLandau2HighestHitProfile_2OutOf10'))
        self.histo2D_1_fid = 0
        self.histo2D_2_fid = 0
        self.histo2D_1.GetXaxis().SetTitle('X/\mu m')
        self.histo2D_1.GetYaxis().SetTitle('Y/\mu m')
        self.histo2D_2.GetXaxis().SetTitle('X/\mu m')
        self.histo2D_2.GetYaxis().SetTitle('Y/\mu m')
        self.histo1D_1 = 0
        self.histo1D_2 = 0
        self.sel_old = {'x_low': self.histo2D_1.GetXaxis().GetXmin(), 'x_high': self.histo2D_1.GetXaxis().GetXmax(), 'y_low': self.histo2D_1.GetYaxis().GetXmin(), 'y_high': self.histo2D_1.GetYaxis().GetXmax()}
        self.fidcut_1 = 0
        self.fidcut_2 = 0
        self.fidpoints = []
        self.nameFid = ''
        if not os.path.isdir('{dir}/Plots'.format(dir=self.outputDir)):
            os.makedirs('{dir}/Plots'.format(dir=self.outputDir))
        if not os.path.isdir('{dir}/root'.format(dir=self.outputDir)):
            os.makedirs('{dir}/root'.format(dir=self.outputDir))
        gStyle.SetPalette(55)
        gStyle.SetNumberContours(999)
        self.bla = []

    def Get1DHisto(self):
        if self.fidcut_1 == 0 or self.fidcut_2 == 0:
            self.fidcut_1, self.fidcut_2 = 0, 0
            self.CreateFidCut()

        self.histo1D_1 = TH1F('hChargeTransp1OutOf10VsFidCut_z', 'hChargeTransp1OutOf10VsFidCut_z', 64, 0, 4096)
        nbins2D_1 = (self.histo2D_1.GetNbinsY()*self.histo2D_1.GetNbinsX()+1)
        for bin in xrange(1, nbins2D_1):
            x, y, z = Long(0), Long(0), Long(0)
            self.histo2D_1.GetBinXYZ(bin, x, y, z)
            point = {'x': self.histo2D_1.GetXaxis().GetBinCenter(x), 'y': self.histo2D_1.GetYaxis().GetBinCenter(y)}
            # if self.WindingIsPointInPoly(point):
            if bool(int(self.fidcut_1.IsInside(point['x'], point['y']))):
                self.histo1D_1.Fill(self.histo2D_1.GetBinContent(bin))
        canvas_name_1 = 'c_t_1D_1_out_of_10_{r}_{n}'.format(r=self.run, n=self.nameFid)
        canvas_1 = self.CreateCanvas(canvas_name_1)
        gStyle.SetOptStat('nemr')
        canvas_1.cd()
        self.histo1D_1.SetLineWidth(3)
        self.histo1D_1.GetYaxis().SetTitle('entries')
        self.histo1D_1.Draw()
        name_1 = 'histo1D_charge_transp_1_out_of_10_{r}_{n}'.format(r=self.run, n=self.nameFid)
        self.SaveCanvas(canvas_1, name_1)
        self.bla.append(canvas_1)

        self.histo1D_2 = TH1F('hChargeTransp2OutOf10VsFidCut_z', 'hChargeTransp2OutOf10VsFidCut_z', 64, 0, 4096)
        nbins2D_2 = (self.histo2D_2.GetNbinsY() * self.histo2D_2.GetNbinsX() + 1)
        for bin in xrange(1, nbins2D_2):
            x, y, z = Long(0), Long(0), Long(0)
            self.histo2D_2.GetBinXYZ(bin, x, y, z)
            point = {'x': self.histo2D_2.GetXaxis().GetBinCenter(x), 'y': self.histo2D_2.GetYaxis().GetBinCenter(y)}
            # if self.WindingIsPointInPoly(point):
            if bool(int(self.fidcut_2.IsInside(point['x'], point['y']))):
                self.histo1D_2.Fill(self.histo2D_2.GetBinContent(bin))
        canvas_name_2 = 'c_t_1D_2_out_of_10_{r}_{n}'.format(r=self.run, n=self.nameFid)
        canvas_2 = self.CreateCanvas(canvas_name_2)
        gStyle.SetOptStat('nemr')
        canvas_2.cd()
        self.histo1D_2.SetLineWidth(3)
        self.histo1D_2.GetYaxis().SetTitle('entries')
        self.histo1D_2.Draw()
        name_2 = 'histo1D_charge_transp_2_out_of_10_{r}_{n}'.format(r=self.run, n=self.nameFid)
        self.SaveCanvas(canvas_2, name_2)
        self.bla.append(canvas_2)

    def Get2DMap(self):
        if self.fidcut_1 == 0 or self.fidcut_2 == 0:
            self.fidcut_1, self.fidcut_2 = 0, 0
            self.CreateFidCut()
        self.histo2D_1.GetXaxis().SetRangeUser(self.sel_old['x_low'], self.sel_old['x_high'])
        self.histo2D_1.GetYaxis().SetRangeUser(self.sel_old['y_low'], self.sel_old['y_high'])
        self.histo2D_2.GetXaxis().SetRangeUser(self.sel_old['x_low'], self.sel_old['x_high'])
        self.histo2D_2.GetYaxis().SetRangeUser(self.sel_old['y_low'], self.sel_old['y_high'])

        self.histo2D_1.SetTitle(self.histo2D_1.GetName())
        self.histo2D_1.GetZaxis().SetTitle('charge/ADC')
        self.histo2D_1.GetZaxis().SetRangeUser(0,3000)
        self.histo2D_2.SetTitle(self.histo2D_2.GetName())
        self.histo2D_2.GetZaxis().SetTitle('charge/ADC')
        self.histo2D_2.GetZaxis().SetRangeUser(0,3000)

        canvas_name_1 = 'c_t_2D_1_out_of_10_{r}_{n}'.format(r=self.run, n=self.nameFid)
        canvas_1 = self.CreateCanvas(canvas_name_1)
        gStyle.SetOptStat('ne')
        canvas_1.cd()
        self.histo2D_1.Draw('colz')
        self.fidcut_1.Draw('same')
        self.bla.append(canvas_1)
        name_1 = 'histo2D_charge_transp_1_out_of_10_{r}_{n}'.format(r=self.run, n=self.nameFid)
        self.SaveCanvas(canvas_1, name_1)

        canvas_name_2 = 'c_t_2D_2_out_of_10_{r}_{n}'.format(r=self.run, n=self.nameFid)
        canvas_2 = self.CreateCanvas(canvas_name_2)
        gStyle.SetOptStat('ne')
        canvas_2.cd()
        self.histo2D_2.Draw('colz')
        self.fidcut_2.Draw('same')
        self.bla.append(canvas_2)
        name_2 = 'histo2D_charge_transp_2_out_of_10_{r}_{n}'.format(r=self.run, n=self.nameFid)
        self.SaveCanvas(canvas_2, name_2)

    def Get2DMapFiducial(self):
        if self.fidcut_1 == 0 or self.fidcut_2 == 0:
            self.fidcut_1, self.fidcut_2 = 0, 0
            self.CreateFidCut()
        self.histo2D_1.GetXaxis().SetRangeUser(self.sel_old['x_low'], self.sel_old['x_high'])
        self.histo2D_1.GetYaxis().SetRangeUser(self.sel_old['y_low'], self.sel_old['y_high'])
        self.histo2D_2.GetXaxis().SetRangeUser(self.sel_old['x_low'], self.sel_old['x_high'])
        self.histo2D_2.GetYaxis().SetRangeUser(self.sel_old['y_low'], self.sel_old['y_high'])

        self.histo2D_1.GetZaxis().SetTitle('charge/ADC')
        self.histo2D_1.GetZaxis().SetRangeUser(0, 3000)

        gStyle.SetOptStat('n')
        self.histo2D_1_fid = self.histo2D_1.Clone('{n}_FidRegion'.format(n=self.histo2D_1.GetName()))
        self.histo2D_1_fid.SetTitle(self.histo2D_1_fid.GetName())
        nbinshisto2D_1 = (self.histo2D_1.GetNbinsY() * self.histo2D_1.GetNbinsX() + 1)
        for bin in xrange(1, nbinshisto2D_1):
            x, y, z = Long(0), Long(0), Long(0)
            self.histo2D_1.GetBinXYZ(bin, x, y, z)
            point = {'x': self.histo2D_1.GetXaxis().GetBinCenter(x), 'y': self.histo2D_1.GetYaxis().GetBinCenter(y)}
            # if not self.WindingIsPointInPoly(point):
            if not (bool(int(self.fidcut_1.IsInside(point['x'], point['y'])))):
                self.histo2D_1_fid.SetBinContent(bin, 0)
        canvas_name_1 = 'c_t_2D_1_out_of_10_{r}_fid_{n}'.format(r=self.run, n=self.nameFid)
        canvas_1 = self.CreateCanvas(canvas_name_1)
        gStyle.SetOptStat('n')
        canvas_1.cd()
        self.histo2D_1_fid.GetXaxis().SetTitle('X/\mu m')
        self.histo2D_1_fid.GetYaxis().SetTitle('Y/\mu m')
        self.histo2D_1_fid.Draw('colz')
        self.fidcut_1.Draw('same')
        self.bla.append(canvas_1)
        name_1 = 'histo2D_charge_transp_1_out_of_10_{r}_fid_{n}'.format(r=self.run, n=self.nameFid)
        self.SaveCanvas(canvas_1, name_1)

        gStyle.SetOptStat('n')
        self.histo2D_2_fid = self.histo2D_2.Clone('{n}_FidRegion'.format(n=self.histo2D_2.GetName()))
        self.histo2D_2_fid.SetTitle(self.histo2D_2_fid.GetName())
        nbinshisto2D_2 = (self.histo2D_2.GetNbinsY() * self.histo2D_2.GetNbinsX() + 1)
        for bin in xrange(1, nbinshisto2D_2):
            x, y, z = Long(0), Long(0), Long(0)
            self.histo2D_2.GetBinXYZ(bin, x, y, z)
            point = {'x': self.histo2D_2.GetXaxis().GetBinCenter(x), 'y': self.histo2D_2.GetYaxis().GetBinCenter(y)}
            # if not self.WindingIsPointInPoly(point):
            if not bool(int(self.fidcut_2.IsInside(point['x'], point['y']))):
                self.histo2D_2_fid.SetBinContent(bin, 0)
        canvas_name_2 = 'c_t_2D_2_out_of_10_{r}_fid_{n}'.format(r=self.run, n=self.nameFid)
        canvas_2 = self.CreateCanvas(canvas_name_2)
        gStyle.SetOptStat('n')
        canvas_2.cd()
        self.histo2D_2_fid.GetXaxis().SetTitle('X/\mu m')
        self.histo2D_2_fid.GetYaxis().SetTitle('Y/\mu m')
        self.histo2D_2_fid.Draw('colz')
        self.fidcut_2.Draw('same')
        self.bla.append(canvas_2)
        name_2 = 'histo2D_charge_transp_2_out_of_10_{r}_fid_{n}'.format(r=self.run, n=self.nameFid)
        self.SaveCanvas(canvas_2, name_2)

    def CreateFidCut(self):
        if self.fidcut_1 == 0 and self.fidcut_2 == 0:
            if len(self.fidpoints) > 3:
                self.ConvertFidSilToTranspSpatial()
                self.fidcut_1 = TCutG('fidcut_1', len(self.fidpoints))
                self.fidcut_1.SetVarX(self.histo2D_1.GetXaxis().GetName())
                self.fidcut_1.SetVarY(self.histo2D_1.GetYaxis().GetName())
                self.fidcut_2 = TCutG('fidcut_2', len(self.fidpoints))
                self.fidcut_2.SetVarX(self.histo2D_2.GetXaxis().GetName())
                self.fidcut_2.SetVarY(self.histo2D_2.GetYaxis().GetName())
                npoints = 0
                for ipoint in self.fidpoints:
                    self.fidcut_1.SetPoint(npoints, ipoint['x'], ipoint['y'])
                    self.fidcut_2.SetPoint(npoints, ipoint['x'], ipoint['y'])
                    npoints += 1
            else:
                print 'The fiducial region is not specified correctly. Enter the region again.'
                self.fidpoints = []
                self.SetFidPoints()
        else:
            print 'The fiducial region has already been defined'
            return
        if self.fidcut_1 != 0 and self.fidcut_2 != 0:
            self.fidcut_1.SetLineColor(kRed)
            self.fidcut_1.SetLineWidth(3)
            self.fidcut_2.SetLineColor(kRed)
            self.fidcut_2.SetLineWidth(3)

    def ConvertFidSilToTranspSpatial(self):
        print 'Converting fiducial region in silicon space to transparent physical space...'
        for ipoint in xrange(len(self.fidpoints)):
            self.fidpoints[ipoint]['x'] = self.fidpoints[ipoint]['x']*50+179.06
            self.fidpoints[ipoint]['y'] = self.fidpoints[ipoint]['y']*50-9.55

    def SetFidPoints(self, npoint=0):
        if npoint == 0:
            if not bool(int(raw_input('Change default Fiducial Region? (yes: 1 / no: 0)'))):
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


    def CreateCanvas(self, name='c0', w=1024, h=1024):
        c = TCanvas(name, name, w, h)
        c.SetWindowSize(w + (w - c.GetWw()), h + (h - c.GetWh()))
        return deepcopy(c)

    def SaveCanvas(self, canvas, name):
        canvas.SaveAs(self.outputDir + '/root/{n}.root'.format(n=name))
        canvas.SaveAs(self.outputDir + '/Plots/{n}.png'.format(n=name))

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

    z = AnalyseTransparentArea(run, source, output)
