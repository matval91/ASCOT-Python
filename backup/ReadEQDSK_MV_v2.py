"""
NAME: ReadEQDSK

AUTHOR: Leonardo Pigatto

DESCRIPTION
Python function to read EQDSK files

CALLING SEQUENCE
out = ReadEQDSK(in_filename)

CHANGE HISTORY:
-started on 29/09/2015 - porting from Matlab

NOTES:


"""

import numpy as np
from pylab import *
import re

class eqdsk():
	def __init__(self,comment,switch,nrbox,nzbox,rboxlength,zboxlength,R0EXP,rboxleft,Raxis,Zaxis,psiaxis,psiedge,B0EXP,Ip,T,p,TTprime,pprime,psi,q,nLCFS,nlimits,R,Z,R_limits,Z_limits,R_grid,Z_grid,psi_grid,rhopsi):
		self.comment=comment
		self.switch=switch
		self.nrbox=nrbox
		self.nzbox=nzbox
		self.rboxlength=rboxlength
		self.zboxlength=zboxlength
		self.R0EXP=R0EXP
		self.rboxleft=rboxleft
		self.Raxis=Raxis
		self.Zaxis=Zaxis
		self.psiaxis=psiaxis
		self.psiedge=psiedge
		self.B0EXP=B0EXP
		self.Ip=Ip
		self.T=T
		self.p=p
		self.TTprime=TTprime
		self.pprime=pprime
		self.psi=psi
		self.q=q
		self.nLCFS=nLCFS
		self.nlimits=nlimits
		self.R=R
		self.Z=Z
		self.R_limits=R_limits
		self.Z_limits=Z_limits
		self.R_grid=R_grid
		self.Z_grid=Z_grid
		self.psi_grid=psi_grid
		self.rhopsi=rhopsi

def ReadEQDSK(in_filename):

	#in_filename = 'g032138.03100'
	#with open(in_filename) as fin:
	#	data = [ii.strip().split() for ii in fin.readlines()]
	data = []
	fin = open(in_filename,"r")
        el1=''
        ll=len(el1)
        d='[0-9]-[0-9]'
	for ii in fin.readlines():
		for jj in ii.split():
                        if re.finditer(d, jj):
                                for ind_s,s in enumerate(re.finditer(d,jj)):
                                        #length=s.span()[0]+1-ll*ind_s
                                        length=15
                                        if jj[0]=='-':
                                                length+=1
                                        el1=jj[0:length]
                                        if el1[-4]!='E' or el1[-4]!='e':
                                                el1=jj[0:length]
                                        data.append(el1)
                                        jj=jj[len(el1):]
                                        ll=len(el1)
                                        el1=''
                                ll=0
                        
		        data.append(jj)
	fin.close()
	
	#Comments, date, switch and grid dimensions
	comment = data[0]+' '+data[1]+' '+data[2]+' '+data[3]
	switch = data[4]
	nrbox = int(data[5])
	nzbox = int(data[6])
	print('nrbox', nrbox)
        print('nzbox', nzbox)
	#First line
	rboxlength = float(data[7])
	zboxlength = float(data[8])
	R0EXP = float(data[9])
	rboxleft = float(data[10])
	
	#Second line
	Raxis = float(data[12])
	Zaxis = float(data[13])
	psiaxis = float(data[14]) # psi_axis-psi_edge
	psiedge = float(data[15]) # psi_edge-psi_edge (=0)
	B0EXP = float(data[16]) # Normalizing magnetic field in CHEASE
	
	#Third line
	#Ip is first element, all others are already stored
	Ip = float(data[17])
	
	#Fourth line - nothing or already stored
	
	#T (or F - poloidal flux function)
	Nn = 27
	T = []
	for ii in range(0,nrbox):
		T.append(float(data[Nn+ii]))
	Nn = Nn+nrbox
	#p (pressure)
	p = []
	for ii in range(0,nrbox):
		p.append(float(data[Nn+ii]))
	Nn = Nn+nrbox
	#TT'
	TTprime = []
	for ii in range(0,nrbox):
		TTprime.append(float(data[Nn+ii]))
	Nn = Nn+nrbox
	#p'
	pprime = []
	for ii in range(0,nrbox):
		pprime.append(float(data[Nn+ii]))
	Nn = Nn+nrbox
	#psi
	psi = []
	for ii in range(0,nrbox*nzbox):
		psi.append(float(data[Nn+ii]))
	psi = np.array(psi)
	psi = np.reshape(psi,(nrbox,nzbox))
	Nn = Nn+(nrbox*nzbox)
	#q safety factor
	q = []
	for ii in range(0,nrbox):
		q.append(float(data[Nn+ii]))
        q=np.absolute(q)
	Nn = Nn+nrbox
	
	#n of points for the LCFS and limiter boundary
	nLCFS = int(data[Nn])
	Nn = Nn+1
	nlimits = int(data[Nn])
	Nn = Nn+1
	#RZ LCFS coordinates
	RZtmp = []
	for ii in range(0,nLCFS*2):
		RZtmp.append(float(data[Nn+ii]))
	R = RZtmp[0:-1:2]
	Z = RZtmp[1::2]
	Nn = Nn+nLCFS*2
	#RZ limits
	RZtmp = []
	for ii in range(0,nlimits*2):
		RZtmp.append(float(data[Nn+ii]))
	R_limits = RZtmp[0:-1:2]
	Z_limits = RZtmp[1::2]
	Nn = Nn+nlimits*2
	#RZ grid for psi
	R_grid = np.linspace(rboxleft,rboxleft+rboxlength,nrbox)
	Z_grid = np.linspace(-zboxlength/2.,zboxlength/2.,nzbox)
	#psi grid for radial prifiles
	psi_grid = np.linspace(psiaxis,psiedge,nrbox)
	#corresponding rhopsi
	rhopsi = sqrt(abs(psi_grid-psiaxis)/abs(psiedge-psiaxis))
	
	out = eqdsk(comment,switch,nrbox,nzbox,rboxlength,zboxlength,R0EXP,rboxleft,Raxis,Zaxis,psiaxis,psiedge,B0EXP,Ip,T,p,TTprime,pprime,psi,q,nLCFS,nlimits,R,Z,R_limits,Z_limits,R_grid,Z_grid,psi_grid,rhopsi)
	return out
