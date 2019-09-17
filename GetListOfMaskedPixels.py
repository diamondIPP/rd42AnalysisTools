#!/usr/bin/env python
import os, sys, shutil
sys.path.append('/home/sandiego/rd42AnalysisTools')  # TODO: HARDCODED!!! NEEDED TO RUN IN BATCH!!! CHANGE ACCORDINGLY
from optparse import OptionParser
import ROOT as ro
from collections import OrderedDict

__author__ = 'DA'
columns = 52
rows = 80



def Correct_Path(path, times=2):
	abs_path = ''
	if path[0] == '~':
		abs_path += os.path.expanduser('~')
		abs_path += path[1:]
	elif os.path.isabs(path):
		abs_path += path
	else:
		abs_path += os.path.abspath(path)
	if times != 1:
		return Correct_Path(abs_path, 1)
	return abs_path

def IsA2DMap(obj):
	if obj.InheritsFrom(ro.TProfile2D.Class().GetName()):
		return True
	elif obj.InheritsFrom(ro.TH2D.Class().GetName()):
		return True
	elif obj.InheritsFrom(ro.TH2F.Class().GetName()):
		return True
	elif obj.InheritsFrom(ro.TH2I.Class().GetName()):
		return True
	else:
		return False

def CheckIfRepeatedLines(filef):
	with open(filef, 'r') as f0:
		with open(filef + '.tmp', 'w') as ftemp:
			lines = f0.readlines()
			lines_params = []
			for line in lines:
				if line.startswith('pix'):
					lines_params.append('{i} {c} {r}'.format(i=line.split(' ')[1], c=line.split(' ')[2], r=line.split(' ')[3]))
				else:
					ftemp.write(line)
			# lines_params2 = list(set(lines_params))
			lines_params2 = list(OrderedDict.fromkeys(lines_params))
			for line2 in lines_params2:
				ftemp.write('pix {i} {c} {r}'.format(i=line2.split(' ')[0], c=line2.split(' ')[1], r=line2.split(' ')[2]))
	shutil.move(filef + '.tmp', filef)

if __name__ == '__main__':
	parser = OptionParser()
	parser.add_option('-f', '--file', dest='file', type='string', help='Path to the root file containing the histogram')
	parser.add_option('-t', '--threshold', dest='threshold', type='float', default=9, help='Threshold to use for each pixel')
	parser.add_option('-g', '--greater', dest='isgreaterthan', default=False, action='store_true', help='the trheshold is applied to greater than the value instead of default lower than the value')
	parser.add_option('-i', '--i2c', dest='i2c', type='int', default=0, help='i2c of device')

	(options, args) = parser.parse_args()
	filef = Correct_Path(str(options.file))
	dirf = '/'.join(filef.split('/')[:-1])
	th = float(options.threshold)
	i2c = int(options.i2c)
	dirt = 1 if bool(options.isgreaterthan) else -1

	print 'Taking file {f}. Marking any pixel in the device with i2c {i} with a threshold {dt} or equal than {v} in the defaultMaskfile.dat'.format(f=filef, i=i2c, dt='greater' if dirt > 0 else 'lower', v=th)

	blaf = ro.TFile(filef)
	blac = blaf.Get(blaf.GetListOfKeys()[0].GetName())
	blah = None

	for it in xrange(blac.GetListOfPrimitives().GetSize()):
		if IsA2DMap(blac.GetListOfPrimitives()[it]):
			blah = blac.GetListOfPrimitives()[it]
			break

	if blah:
		with open(dirf + '/defaultMaskFile.dat', 'a') as f:
			for col in xrange(columns):
				for row in xrange(rows):
					if dirt < 0:
						if blah.GetBinContent(col + 1, row + 1) <= th:
							print 'pix', i2c, col, row
							f.write('pix {i} {c} {r}\n'.format(i=i2c, c=col, r=row))
					else:
						if blah.GetBinContent(col + 1, row + 1) >= th:
							print 'pix', i2c, col, row
							f.write('pix {i} {c} {r}\n'.format(i=i2c, c=col, r=row))

		CheckIfRepeatedLines(dirf + '/defaultMaskFile.dat')

