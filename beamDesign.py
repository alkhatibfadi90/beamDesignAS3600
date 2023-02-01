"""
Analysis and Design of Reinforced Concrete Beam to AS3600:2018
By: Fadi Alkhatib
"""

import AS3600_2018 as AS
from sympy import*
from math import*
import numpy as np
import beamSketch
import matplotlib.pyplot as plt
import matplotlib.patches as pat


class Rect:
	"""docstring for """
	def __init__(self,name="B1",fc=32, fy=420):
		self.name = name
		self.code = 'AS3600:2018'
		self.fc = fc
		self.fy = fy
		self.Ec = AS.concreteModulus(fc)
		self.Es = AS.steelModulus
		self.unitLength = 'mm'
		self.unitStress = 'MPa'
		self.barloc = 'bottom'
		self.wc = AS.wc
		self.f_ct, self.f_ctf = AS.concUnixialTensile(fc)


	def dimension(self,breadth=300,height=600, cover = 30, length= None, lengthSpan = None, ):
		self.bw = breadth
		self.B = breadth
		self.bv = breadth
		self.beff = breadth
		self.ds = cover
		self.h = height
		self.L = length
		self.Lsp = lengthSpan
		self.Ag = self.bw * self.h
		self.Ig = self.bw * self.h ** 3 / 12


	def rebars (self, barClass = 'N', tensionBarDia = 16, tensionBarNo= 3,compBarDia =10, compBarNo=2, linkDia = 10, linkNo= 2, linkSpacing = 200):
		self.barClass = barClass
		self.diaT = tensionBarDia
		self.noT = tensionBarNo
		self.diaC = compBarDia
		self.noC = compBarNo
		self.diaV = linkDia
		self.legsV = linkNo
		self.Vs = linkSpacing
		self.d = round(self.h - (self.ds + self.diaV + self.noT/2),2)
		self.d1 = round(self.h - self.d - self.diaC / 2, 2)
		self.d_v = round(max(0.72 * self.h, 0.9 * self.d))
		self.HsT = round((self.bw - (2 *(self.ds) + self.diaT)) / (self.noT - 1),2)
		self.HsTc = round((self.bw - (2 * (self.ds) + self.noT * self.noT)) / (self.noT - 1),2)
		self.HsC = round((self.bw - (2 *(self.ds) + self.diaC)) / (self.noC - 1),2)


	def loading (self, SDL = 0, LL = 0, Gk = 1.5, Qk = 1.25):
		self.SDL = SDL * self.Lsp
		self.LL = LL * self.Lsp
		self.Gk = Gk
		self.Qk = Qk
		self.selfweight = self.bw/1000 * self.h/1000 * self.wc
		self.DL = self.SDL + self.selfweight
		self.comboSLS = '1.0 DL  + 1.0 LL'
		self.comboULS = '1.25 DL + 1.5 LL'
		self.SLS = self.DL + self.LL
		self.ULS= self.Gk * self.DL + self.Qk * self.L


	def analysis (self, plot=False):
		interval = 0.1
		xlist = np.arange(0, self.L + interval, interval).tolist()
		Mx = [- ((self.ULS * x) / 2) * (self.L - x) for x in xlist]
		Vx = [self.ULS * (self.L / 2 - x) for x in xlist]
		self.Moment = max(Mx, key = abs)
		self.Shear = max(Vx, key = abs)
		if plot == True:
			beamSketch.analysisPlot(interval, self.L, self.ULS, Vx, Mx)


	def design (self):
		#Flexural
		self.gamma = AS.gamma(self.fc)
		self.alpha_2 = AS.alpha_2(self.fc)
		self.alpha_b = AS.alpha_b
		self.Ast = AS.barArea(self.diaT) * self.noT
		self.rhoT = self.Ast / (self.bw * self.h)
		self.Ast_min = ((self.alpha_b * (self.h / self.d) ** 2 * self.f_ctf / self.fy) * self.B * self.d)
		self.Asc = AS.barArea(self.diaC) * self.noC
		self.k_u = (self.Ast * self.fy) / (self.alpha_2 * self.fc * self.gamma * self.bw * self.d + self.Asc * self.fy)
		self.phi_moment = AS.momentReductionFactor(self.barClass, self.k_u)

		#Shear
		self.alpha_v = 90 #angle of shear fitment: 90 for perpendicular
		self.theta_v = 0.36 #for perpendicular fitment
		self.Asv = AS.barArea(self.diaV) * self.legsV
		self.Asv_min = ((0.08 * sqrt(self.fc) * self.bv * self.Vs) / self.fy)  # Cl 8.2.1.7
		self.k_v = AS.shearKvalue(self.Asv_min, self.Asv, self.Vs, self.d_v)
		self.phi_shear = AS.shearReductionFactor(self.barClass, self.Asv_min, self.Asv)

		#Deflection
		beta = self.beff / self.bw
		rhoFactor = (0.001 * pow(self.fc, (1/3)))/ pow(beta, (2/3))

		I_eff1 = ((5 - 0.04 * self.fc) * self.rhoT + 0.002) * self.beff * pow(self.d, 3)
		I_eff2 = (0.1 / pow(beta, 3)) * self.beff * pow(self.d, 3)
		I_eff3 = ((0.055 * (pow(self.fc, (1 / 3)))) / pow(beta, (2 / 3))) * self.beff * pow(self.d, 3)
		I_eff4 = (0.06 / pow(beta, (2 / 3))) * self.beff * pow(self.d, 3)

		if self.rhoT >= rhoFactor:
			if I_eff1 <= I_eff2:
				self.I_eff = I_eff1
			else:
				self.I_eff = I_eff2
		else:
			if I_eff3 <= I_eff4:
				self.I_eff = I_eff3
			else:
				self.I_eff = I_eff4

		self.deflect = (5 * self.SLS * pow(self.L, 4)) / (384 * self.Ec * self.I_eff)

		#Long-term deflection
		self.K_cs = (2 - 1.2 * (self.Asc / self.Ast)) #long term factor must >= 0.8  Cl.8.5.3.2
		if self.K_cs >= 0.8:
			self.longtermdeflect = 'Long term deflection is satisfied'
		else:
			self.longtermdeflect = ' long-term deflection is NOT satisfied'

		#Cracking
		if self.Ast < self.Ast_min or self.ds > 100 or self.HsT > 300:
			self.cracking = 'Cracking fail: increase reinforcement or reduce spacing and cover'
		else:
			self.cracking = 'Cracking: Satisfied'

		#Stability
		self.stabilityFactor = self.L / self.beff
		if self.stabilityFactor <= min((180 * self.beff / 60), 60):
			self.stability = 'Stable'
		else:
			self.stability = 'Not Stable: Increase beam breadth'


	def plotting(self, size = (10,20), type='info'):
		self.fig, self.axs = plt.subplots(figsize=size)
		fsize = 10
		loc = (0, 0)
		infoSpace = 600
		pad = 100
		mg = 100
		xLeft = - pad * 2
		xRight = self.bw + infoSpace + pad
		yBottom = - pad * 2
		yTop = self.h + pad  * 2
		x0 = loc[0] + pad
		y0 = loc[1] + pad

		#Plot title
		self.axs.vlines(xLeft - pad, yBottom, yTop + pad, linewidth=0.1)
		self.axs.hlines(yBottom, xLeft, xRight, linewidth=0.1)
		if type =='info':
			title = 'Design Info'
		else:
			title = 'Design Detials'
		self.axs.text(x0+self.bw/2, y0+self.h+pad, f"{title}", va='center', fontweight='bold', fontsize=fsize+2)
		xtxt = self.bw + pad * 2
		ytxt = self.h + pad
		self.axs.text(xtxt, ytxt - (1 * mg), f"{self.name}({self.bw}{self.unitLength} x {self.h}{self.unitLength}) ", va='center', fontweight='bold',
				 fontsize=fsize)

		#Plot beam section and dimensions
		beamSketch.beamPlot(x0, y0, self.bw, self.h)
		beamSketch.dimPlot(x0, y0, self.bw, self.h, pad)

		if type == 'info':
			plt.text(xtxt, ytxt - (2 * mg), f"Design Code: {self.code}", va='center', fontsize=fsize)
			plt.text(xtxt, ytxt - (3 * mg), f"$f_c: {self.fc} \ {self.unitStress}$", va='center',fontsize=fsize)
			plt.text(xtxt, ytxt - (4 * mg), f"$f_sy$: {self.fc}  {self.unitStress}$", va='center', fontsize=fsize)
			plt.text(xtxt, ytxt - (5 * mg), f"$f_sy$: {self.fy} {self.unitStress}", va='center', fontsize=fsize)
			plt.text(xtxt, ytxt - (6 * mg), f"$E_c$: {self.Ec} {self.unitStress}", va='center',  fontsize=fsize)
			plt.text(xtxt, ytxt - (7 * mg), f"$E_s$: {self.Es} {self.unitStress}", va='center', fontsize=fsize)

		if type == 'detail':
			#Plot reinforcement
			beamSketch.rebarPlot(x0,y0,self.ds,self.bw,self.h,self.diaT,self.noT,self.HsT,self.diaC,self.noC,self.HsC,self.diaV)

			plt.text(xtxt, ytxt - (2 * mg), f"Cover: {self.ds}{self.unitLength}", va='center', fontsize=fsize)
			plt.text(xtxt, ytxt - (3 * mg), f"Tension Reinforcement: {self.noT}{self.barClass}-{self.diaT}", va='center', fontsize=fsize)
			plt.text(xtxt, ytxt - (4 * mg), f"Compression Reinforcement: {self.noC}{self.barClass}-{self.diaC}", va='center', fontsize=fsize)
			plt.text(xtxt, ytxt - (5 * mg), f"Shear Reinforcement: {self.barClass}{self.diaV}-{self.Vs}", va='center', fontsize=fsize)

		self.ax = plt.gca()
		self.ax.set_aspect('equal', adjustable='box')
		self.ax.axis("off")
		plt.show()



