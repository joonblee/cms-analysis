# Check : You should have 4 muon selection result files of all signal samples

import sys, os, glob, time, datetime

if __name__ == "__main__":
  os.system('rm *.root')
  os.system('rm *.png')

  M_A = [1]

  for i in M_A:

    os.system('''root -l 'DCBfit.cc(200,%s)' &>log.H200A%s&''' %(i,i))
    os.system('''root -l 'DCBfit.cc(300,%s)' &>log.H300A%s&''' %(i,i))
    os.system('''root -l 'DCBfit.cc(400,%s)' &>log.H400A%s&''' %(i,i))
    os.system('''root -l 'DCBfit.cc(500,%s)' &>log.H500A%s&''' %(i,i))
    os.system('''root -l 'DCBfit.cc(750,%s)' &>log.H750A%s&''' %(i,i))
    os.system('''root -l 'DCBfit.cc(1000,%s)' &>log.H1000A%s&''' %(i,i))
    os.system('''root -l 'DCBfit.cc(1250,%s)' &>log.H1250A%s&''' %(i,i))
    os.system('''root -l 'DCBfit.cc(1500,%s)' &>log.H1500A%s&''' %(i,i))
    os.system('''root -l 'DCBfit.cc(1750,%s)' &>log.H1750A%s&''' %(i,i))
    os.system('''root -l 'DCBfit.cc(2000,%s)' &>log.H2000A%s&''' %(i,i))
    time.sleep(5)

  files = []
  path_ = './'

  for i in range(0,10):
    files = []
    for iDir in glob.glob(os.path.join(path_,'*.root')):
      files.append(iDir)

    print '- # of root files : ',len(files)
    if len(files) == 10*len(M_A):
      break

    time.sleep(5)

  print '\n +++ root file list +++ '
  for iDir in files:
    print iDir
  print 

  for i in M_A:
    os.system('hadd sigma_phi%s.root sigma_*_%s.root' %(i,i))
    os.system('''root -l 'hist_res.cc("phi%s")' &''' %i)






