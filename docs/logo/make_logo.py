import numpy as np, math
NAVY='#14233F'; ORANGE='#F47230'; BLUE='#4FA3F7'; LIGHT='#E8EEF7'
def arrow(cx,cy,ang,L,col,w=5,head=9):
    # ang in degrees, 0 = +x; arrow centred on (cx,cy)
    a=math.radians(ang); dx,dy=math.cos(a)*L/2,-math.sin(a)*L/2
    x0,y0,x1,y1=cx-dx,cy-dy,cx+dx,cy+dy
    ux,uy=(x1-x0)/L,(y1-y0)/L; px,py=-uy,ux
    bx,by=x1-ux*head,y1-uy*head
    shaft=f'<line x1="{x0:.1f}" y1="{y0:.1f}" x2="{bx:.1f}" y2="{by:.1f}" stroke="{col}" stroke-width="{w}" stroke-linecap="round"/>'
    tri=f'<polygon points="{x1:.1f},{y1:.1f} {bx+px*head*0.62:.1f},{by+py*head*0.62:.1f} {bx-px*head*0.62:.1f},{by-py*head*0.62:.1f}" fill="{col}"/>'
    return shaft+tri
def icon(ox=0,oy=0,s=1.0):
    g=[f'<g transform="translate({ox},{oy}) scale({s})">']
    g.append(f'<rect x="0" y="0" width="200" height="200" rx="38" fill="{NAVY}"/>')
    # double-well potential E = x^4 - 2x^2 on x in [-1.65,1.65]
    xs=np.linspace(-1.62,1.62,200); E=xs**4-2*xs**2
    X=100+xs*50; Y=100-(E+1)*27
    d='M'+' L'.join(f'{x:.1f},{y:.1f}' for x,y in zip(X,Y))
    g.append(f'<path d="{d}" fill="none" stroke="{BLUE}" stroke-width="6" stroke-linecap="round" stroke-linejoin="round"/>')
    # DW in left well
    xm=-1.0; g.append(f'<circle cx="{100+xm*50:.1f}" cy="{100-12:.1f}" r="12" fill="{ORANGE}"/>')
    # nanowire strip with head-to-head transverse wall
    g.append(f'<rect x="22" y="128" width="156" height="44" rx="10" fill="{LIGHT}"/>')
    xs_=[40,64,88,112,136,160]; angs=[0,0,60,120,180,180]
    for x,a in zip(xs_,angs):
        col=ORANGE if a not in (0,180) else NAVY
        g.append(arrow(x,150,a,22,col,w=4.6,head=9))
    g.append('</g>'); return ''.join(g)
FONT="font-family=\"Inter Display, Inter, Fira Sans, DejaVu Sans, sans-serif\""
open('/tmp/logo/dwpc_icon.svg','w').write(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 200 200" width="400" height="400">{icon()}</svg>')
lock=f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 620 200" width="1240" height="400">
{icon()}
<text x="232" y="118" {FONT} font-size="104" font-weight="800" fill="{NAVY}" letter-spacing="2">DWPC</text>
<text x="236" y="160" {FONT} font-size="26" font-weight="500" fill="#4A5670" letter-spacing="0.5">Domain Wall Propagation Code</text>
</svg>'''
open('/tmp/logo/dwpc_logo.svg','w').write(lock)
