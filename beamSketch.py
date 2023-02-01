import matplotlib.pyplot as plt
import matplotlib.patches as pat
import numpy as np


bw = 400
h = 800
ds = 50
diaT = 20
noT = 6
space = (bw - 2*ds - diaT) / (noT-1)
noC = 2
diaC = 16
spacec = (bw - 2*ds - diaC) / (noC-1)
diaV = 10
unitl = 'mm'
w = 10
L = 5



#General Boundaries


def beamPlot(x0,y0,bw,h):
    rect = pat.Rectangle((x0,y0), bw, h, fill=True, alpha=0.3, capstyle='round', edgecolor='black', linestyle = '-', facecolor = 'gray', linewidth =1)
    plt.gca().add_patch(rect)


def rebarPlot(x0,y0,ds,bw,h,diaT,noT,Hs,diaC,noC,Hcs,diaV,):
    fac = pat.FancyBboxPatch((x0+ds,y0+ds), bw-2*ds, h-2*ds, fill = False, boxstyle="round,pad=20",  facecolor = 'yellow', linewidth=2)
    plt.gca().add_patch(fac)
    barsT = []
    for i in range(noT):
        barSpace = i * Hs
        cir = pat.Circle((x0+ds+diaT/2 + barSpace, y0+ds+diaT/2), radius=diaT, color='red', alpha=1)
        barsT.append(cir)
    barsC = []
    for i in range(noC):
        barSpace = i * Hcs
        cir = pat.Circle((x0+ds+diaC/2 + barSpace, y0+h-ds-diaC/2), radius=diaC, color='red', alpha=1)
        barsC.append(cir)
    bars = barsT + barsC
    for i in (bars):
        plt.gca().add_patch(i)

def dimPlot(x0,y0,bw, h, pad):

    #Horizontal
    startH = (x0,y0 - pad)
    endH = (x0 + bw, y0 - pad)
    txtH = (x0 + bw/2, y0 - pad)
    alvH = 'bottom'
    alhH = 'center'
    plt.annotate(text='', xy=startH, xytext=endH, arrowprops=dict(arrowstyle='|-|', linewidth=1))
    plt.annotate(text='', xy=startH, xytext=endH, arrowprops=dict(arrowstyle='<|-|>', linewidth=1))
    plt.annotate(text=f'{bw}{unitl}', xy=txtH, rotation=0, va=alvH, ha=alhH, fontsize=10)

    #Vertical
    startV = (x0 - pad,y0)
    endV = (x0 - pad, y0 + h)
    txtV = (x0 -pad, y0 + h/2)
    alvV = 'center'
    alhV = 'right'
    plt.annotate(text='', xy=startV, xytext=endV, arrowprops=dict(arrowstyle='|-|', linewidth=1))
    plt.annotate(text='', xy=startV, xytext=endV, arrowprops=dict(arrowstyle='<|-|>', linewidth=1))
    plt.annotate(text=f'{h}{unitl}', xy=txtV, rotation=0, va=alvV, ha=alhV, fontsize=10)


def analysisPlot(interval, L, ULS, Vx, Mx, figsize=(7,9)):
    figsize = figsize
    pad = 70
    txtpad = 20
    fig, axs = plt.subplots(3,1, figsize=figsize, sharex=True)


    xlist = np.arange(0, L + interval, interval).tolist()
    # plt.xlim(- interval, L + interval)
    xtickslocation = xlist[::5]
    plt.xticks(xtickslocation)
    plt.xlim(- interval * 5, L + interval * 5)

    #Loading
    axs[0].vlines(xlist[::2], ULS, 0, linewidth=1,colors='red', linestyles='solid')
    endArrow = [10 for x in xlist[::2]]
    axs[0].scatter(xlist[::2], endArrow, marker=11, color='red')
    axs[0].text(L / 2, ULS + txtpad, f"w = {ULS} kN/m", fontsize=12, ha='center', va='center')
    axs[0].hlines(ULS, 0, L, colors='red', linewidth=1)
    axs[0].hlines(0, 0, L, colors='black', linewidth=1)
    axs[0].plot(0, -10, marker='^', ms=10, mfc='blue', mec='darkblue')
    axs[0].plot(L, -10, marker='o', ms=10, mfc='blue', mec='darkblue')
    axs[0].set_ylim([- ULS - pad, ULS + pad])
    axs[0].grid()

    #Shear
    axs[1].plot(xlist, Vx, color='red')

    axs[1].vlines(0, 0, Vx[0], colors='red', linewidth=1)
    axs[1].text(0, Vx[0] + txtpad, f"{round(Vx[-1],2)} kN", fontsize=12, ha='center', va='center')

    axs[1].vlines(L, 0, Vx[-1], colors='red', linewidth=1)
    axs[1].text(L , Vx[-1]- txtpad, f"{round(Vx[-1],2)}kN", fontsize=12, ha='center', va='center')

    axs[1].hlines(0,0,L,colors='black', linewidth=1)

    Vxmax = max(Vx,key=abs)
    Vxmax = abs(Vxmax)
    axs[1].set_ylim([- Vxmax - pad , Vxmax + pad])
    axs[1].grid()

    #Moment
    axs[2].plot(xlist, Mx, color='red', linewidth=1)
    axs[2].text( xlist[Mx.index(min(Mx))], min(Mx) - txtpad, f"{round(max(Mx, key=abs),2)} kN.m", fontsize=12, ha='center', va='center')
    axs[2].hlines(0, 0, L, linewidth=1, colors='black')
    Mxmax = max(Mx, key=abs)
    Mxmax = abs(Mxmax)
    print(Mxmax)
    axs[2].set_ylim([-Mxmax - pad, Mxmax + pad])
    axs[2].grid()
    # ax = plt.gca()
    # ax.set_aspect('equal', adjustable='box')
    plt.show()






