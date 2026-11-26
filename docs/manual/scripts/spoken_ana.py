import numpy as np, json, sys, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
C=['#257de0','#f47230','#00b894','#f2a600']; ink,ink2='#191919','#666666'
plt.rcParams.update({'font.size':9,'axes.edgecolor':'#bbbbbb','axes.labelcolor':ink,'xtick.color':ink2,'ytick.color':ink2,'font.family':'serif','mathtext.fontset':'cm'})
D=sys.argv[1]; tag=sys.argv[2]; lam=2.1
raw=np.loadtxt('/tmp/rc/spoken/Signal_in.data'); U=raw[:,:1024]; lab=raw[:,1024].astype(int)
X=np.loadtxt(D+'/Processed_spoken_signal.txt'); assert X.shape==(2000,1024)
spk=np.repeat(np.arange(4),500); Y=np.eye(10)[lab]
def fit(A,Yt,l): return np.linalg.lstsq(A.T@A+l*l*np.eye(A.shape[1]),A.T@Yt,rcond=1e-10)[0]
def acc(A,W,l_): return int((np.argmax(A@W,1)==l_).sum())
R={'tag':tag}
# 1) exactly what the code reports: train on all 2000, score on same 2000
W=fit(X,Y,lam); R['in_sample']=acc(X,W,lab)
# 2) split.sh-style held-out: within each block of 50 (one digit of one speaker) first 25 train, last 25 test
idx=np.arange(2000); tr=(idx%50)<25; te=~tr
Z=lambda A:np.hstack([(A-A[tr].mean(0))/(A[tr].std(0)+1e-12),np.ones((len(A),1))])
for name,A in [('reservoir',X),('raw_input',U),('reservoir_std',Z(X)),('raw_input_std',Z(U))]:
    out={}
    for l in [1e-3,1e-1,1,2.1,10,100,1e3,1e4]:
        W=fit(A[tr],Y[tr],l); out[str(l)]=acc(A[te],W,lab[te])
    R['split_'+name]=out
# 3) leave-one-speaker-out at lambda giving best split accuracy
for name,A in [('reservoir',X),('raw_input',U)]:
    lbest=float(max(R['split_'+name],key=lambda k:R['split_'+name][k])); res=[]
    for s in range(4):
        W=fit(A[spk!=s],Y[spk!=s],lbest); res.append(acc(A[spk==s],W,lab[spk==s]))
    R['loso_'+name]={'lambda':lbest,'per_speaker':res}
# confusion matrix, split test set, reservoir, best lambda
lb=float(max(R['split_reservoir'],key=lambda k:R['split_reservoir'][k]))
W=fit(X[tr],Y[tr],lb); p=np.argmax(X[te]@W,1); CM=np.zeros((10,10),int)
for t,q in zip(lab[te],p): CM[t,q]+=1
R['cm']=CM.tolist(); print(json.dumps(R)); json.dump(R,open('/tmp/rc/spoken_%s.json'%tag,'w'))
fig=plt.figure(figsize=(6.8,3.2)); gs=fig.add_gridspec(2,2,width_ratios=[1.3,1])
i0=np.where(lab==3)[0][0]; k=np.arange(96)
a=fig.add_subplot(gs[0,0]); a.plot(k,U[i0,:96],color=ink2,lw=0.8,marker='o',ms=2); a.set_ylabel('input'); a.tick_params(labelbottom=False)
a.set_title('first 96 features of one "3" (Jackson)',fontsize=8.5,color=ink)
b_=fig.add_subplot(gs[1,0],sharex=a); b_.plot(k,X[i0,:96],color=C[0],lw=0.8,marker='o',ms=2); b_.set_ylabel('reservoir (nm)'); b_.set_xlabel('feature index $k$')
for q in (a,b_):
    q.grid(color='#eeeeee',lw=0.7); q.set_axisbelow(True)
    for sp in ['top','right']: q.spines[sp].set_visible(False)
ax=[None,fig.add_subplot(gs[:,1])]
a=ax[1]; im=a.imshow(CM,cmap='Blues'); a.set_xticks(range(10)); a.set_yticks(range(10))
a.set_xlabel('predicted digit'); a.set_ylabel('true digit')
for i in range(10):
    for j in range(10):
        if CM[i,j]: a.text(j,i,CM[i,j],ha='center',va='center',fontsize=6.5,color='white' if CM[i,j]>60 else ink)
a.set_title('held-out test set (%d / %d)'%(np.trace(CM),te.sum()),fontsize=8.5,color=ink)
fig.tight_layout(); fig.savefig('/tmp/manual/figures/rc_spoken_%s.pdf'%tag); fig.savefig('/tmp/rc/rc_spoken_%s.png'%tag,dpi=130)
