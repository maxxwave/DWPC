import numpy as np, matplotlib, json
matplotlib.use('Agg'); import matplotlib.pyplot as plt
C=['#257de0','#f47230','#00b894','#f2a600']; ink,ink2='#191919','#666666'
plt.rcParams.update({'font.size':9,'axes.edgecolor':'#bbbbbb','axes.labelcolor':ink,'xtick.color':ink2,'ytick.color':ink2,'font.family':'serif','mathtext.fontset':'cm'})
def style(a):
    a.grid(color='#eeeeee',lw=0.7); a.set_axisbelow(True)
    for s in ['top','right']: a.spines[s].set_visible(False)
Nv=10
def load(D):
    d=np.loadtxt(D+'sine_square_init'); x,y=d[:,0],d[:,1]; n=len(x)//2
    S=np.fromfile(D+'Signal.out',sep=' ').reshape(n,Nv); V=np.fromfile(D+'Valid_Signal.out',sep=' ').reshape(n,Nv)
    yt,yv=y[:n],y[n:]
    al=1e-3; W=np.linalg.solve(S.T@S+al*al*np.eye(Nv),S.T@yt)
    pt,pv=S@W,V@W
    acc=lambda p,t:int(((p>0.5).astype(int)==t.astype(int)).sum())
    return x,yt,S,pt,{'train':acc(pt,yt),'valid':acc(pv,yv),'n':n}
x,yt,S0,p0,r0=load('/tmp/rc/sine0/'); _,_,S1,p1,r1=load('/tmp/rc/sine/')
res={'0K':r0,'100K':r1}; print(res); json.dump(res,open('/tmp/rc/sine_res.json','w'))
S,pt=S0,p0
sl=slice(0,160); k=np.arange(160)
fig,ax=plt.subplots(4,1,figsize=(6.4,6.4),sharex=True)
a=ax[0]; style(a); a.step(k,x[sl],where='mid',color=ink,lw=1.0)
lab=yt[sl]
for i in k:
    if lab[i]>0.5: a.axvspan(i-0.5,i+0.5,color=C[1],alpha=0.12,lw=0)
a.set_ylabel('input $u_k$'); a.text(0.995,1.02,'shaded: target = 1 (square)',transform=a.transAxes,ha='right',va='bottom',fontsize=7.5,color=ink2)
a=ax[1]; style(a)
for j in range(Nv): a.plot(k,S[sl,j],color=C[0],lw=0.6,alpha=0.55)
a.set_ylabel('node state (nm)'); a.text(0.995,1.02,'10 virtual nodes, T = 0 K',transform=a.transAxes,ha='right',va='bottom',fontsize=7.5,color=ink2)
for a,p,T in [(ax[2],p0,'0 K'),(ax[3],p1,'100 K')]:
    style(a); a.step(k,lab,where='mid',color=ink2,lw=0.9,label='target')
    a.plot(k,p[sl],color=C[1],lw=1.2,label='readout')
    a.axhline(0.5,color=ink2,lw=0.6,ls=(0,(3,3))); a.set_ylabel('output'); a.set_ylim(-0.4,1.25)
    a.text(0.01,0.88,'T = '+T,transform=a.transAxes,fontsize=8,color=ink)
ax[2].legend(frameon=False,ncol=2,loc='upper right',fontsize=7.5); ax[3].set_xlabel('sample index $k$')
fig.tight_layout(); fig.savefig('/tmp/manual/figures/rc_sine.pdf'); fig.savefig('/tmp/rc/rc_sine.png',dpi=130); plt.close()
