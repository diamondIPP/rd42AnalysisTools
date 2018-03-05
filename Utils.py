import os, shutil

def IsInt(i):
	try:
		int(i)
		return True
	except ValueError:
		return False

def IsFloat(f):
	try:
		float(f)
		return True
	except ValueError:
		return False

def CreateDefaultSettingsFile(dir, run_no, events):
	if not os.path.isdir(dir):
		os.makedirs(dir)
	with open('settings.{r}.ini'.format(r=run_no), 'w') as f:
		f.write('runNo = {r}\n'.format(r=run_no))
		f.write('Events = {e}\n'.format(e=events))
		f.write('repeaterCardNo = 2\n')
		f.write('voltage = 0\n')
		f.write('diamondName = default\n')
		f.write('currentBegin = 0\n')
		f.write('currentEnd = 0\n\n')
		f.write('diamondPattern = {0,50,0,127}\n\n')
		f.write('Iter_Size = 500\n\n')
		f.write('dia_input = 0\n\n')
		f.write('Dia_channel_screen_channels = {}\n\n')
		f.write('Dia_channel_not_connected = {}\n\n')
		f.write('Si_Pedestal_Hit_Factor = 5\n')
		f.write('Di_Pedestal_Hit_Factor = 3\n\n')
		f.write('DO_CMC = 1\n\n')
		f.write('clusterSeedFactors = {16,20,20,20,18,18,14,10,5}\n')
		f.write('clusterHitFactors = {10,16,14,14,14,14,12,10,3}\n\n')
		f.write('UseAutoFidCut = 0\n')
		f.write('nDiamonds = 1\n')
		f.write('si_avg_fidcut_xlow = 0\n')
		f.write('si_avg_fidcut_xhigh = 255\n')
		f.write('si_avg_fidcut_ylow = 0\n')
		f.write('si_avg_fidcut_yhigh = 255\n')
		f.write('selectionfidCut0 = {0-255,0-255}\n')
		f.write('RerunSelection = 0\n\n')
		f.write('DetectorsToAlign = 2\n')
		f.write('alignment_training_track_number = {e}\n'.format(e=int(round(events)/10.0)))
		f.write('alignment_training_method = 1\n')
		f.write('alignment_chi2 = 5\n\n')
		f.write('Double_t detectorD0Z = 0.725\n')
		f.write('Double_t detectorD1Z = 1.625\n')
		f.write('Double_t detectorD2Z = 18.725\n')
		f.write('Double_t detectorD3Z = 19.625\n')
		f.write('Double_t detectorDiaZ = 10.2\n\n')
		f.write('3dShortAnalysis = 0;\n')
		f.write('3dLongAnalysis = 0;\n')
		f.write('3dTransparentAnalysis = 0;\n')

def CloseSubprocess(p, stdin=False, stdout=False):
	pid = p.pid
	if stdin:
		p.stdin.close()
	if stdout:
		p.stdout.close()
	if p.wait() is None:
		print 'Could not terminate subprocess... forcing termination'
		p.kill()
		if p.wait() is None:
			print 'Could not kill subprocess... quitting'
			exit()
	try:
		os.kill(pid, 0)
	except OSError:
		pass
	else:
		print 'The subprocess is still running. Killing it with os.kill'
		os.kill(pid, 15)
		try:
			os.kill(pid, 0)
		except OSError:
			pass
		else:
			print 'The process does not die... quitting program'
			exit()
	del p, pid
	p = None
