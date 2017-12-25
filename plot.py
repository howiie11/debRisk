from visual import *
import pandas as pd,numpy as np
from sys import exit

###############################################################
#CONSTANTS
###############################################################
BASE_DIR="./"

#SCENE PROPERTIES
scale=2.0 #PC
srange=2*scale #PC
axradius=scale/200.0
ambient=0.1

#VIEW PARAMETERS
phi0=radians(45.0)
dphi=-radians(360)

theta0=radians(45.0)
thetamax=radians(60.0)
dtheta=radians(60)

#SHOW NO SHOW
qaxis=1

#ANIMATION
trate=500
retain=50

###############################################################
#ROUTINES AND CLASSES
###############################################################
"""
Transform the universe coordinates into graphical coordinates.

Thus for instance in the graphical space z is going out from screen.

Input vector is the universe coordinates.
Output vector is the graphical coordinates.

Example:
  To change z-axis for x, yaxis for z and xaxis for y
  Use: (v[1],v[2],v[0])
  
"""
ta=lambda(v):vector(v[1],v[2],v[0])

def setView(i,n,phi0=0.0,theta0=np.pi/2,dphi=2*np.pi,dtheta=0.0,thetamax=0.0,
            anim=True):

    if not anim:i=0
    #Equivalen unit time
    t=(1.*i)/n

    #Angles
    theta=theta0+thetamax*np.sin(dtheta*t)
    phi=phi0+dphi*t

    #Coordinates on unit sphere
    z=-np.cos(theta)
    x=-np.sin(theta)*np.cos(phi)
    y=-np.sin(theta)*np.sin(phi)
    return dict(phi=phi,theta=theta,forward=ta((x,y,z)))

def setLight(phi=0.0,theta=0.0):
    z=np.cos(theta)
    x=np.sin(theta)*np.cos(phi)
    y=np.sin(theta)*np.sin(phi)
    return ta((x,y,z))

def circle(r=1,color=color.white,lw=0.1):
    qs=np.linspace(0,2*np.pi,100)
    path=[(r*np.cos(q),0,r*np.sin(q)) for q in qs]
    c=curve(pos=path,radius=lw)
    return c

###############################################################
#SCENE AND CAMMERA
###############################################################
view=setView(0,1,phi0=phi0,theta0=theta0)
light=setLight(phi=radians(45.0),theta=radians(45.0))
scene = display(title='Satellite Motion',
                x=0,y=0,width=1080, height=720,
                center=(0,0,0),
                background=(0,0,0),
                forward=view["forward"],
                range=srange,
                ambient=color.gray(ambient),
                lights=light
                )
#arrow(pos=vector(0,0,0),axis=scale*light,shaftwidth=axradius,color=color.red)

###############################################################
#AXIS AND SCALES
###############################################################
if qaxis:
    xaxis1=cylinder(pos=(0,0,0),axis=ta((scale,0,0)),radius=axradius)
    xaxis2=cylinder(pos=(0,0,0),axis=ta((-scale,0,0)),radius=axradius)
    xl=label(text='x',align='center',pos=ta((+scale,0,0)),box=False,opacity=0)
    yaxis1=cylinder(pos=(0,0,0),axis=ta((0,scale,0)),radius=axradius)
    yaxis2=cylinder(pos=(0,0,0),axis=ta((0,-scale,0)),radius=axradius)
    yl=label(text='y',align='center',pos=ta((0,+scale,0)),box=False,opacity=0)
    zaxis1=cylinder(pos=(0,0,0),axis=ta((0,0,scale)),radius=axradius)
    zaxis2=cylinder(pos=(0,0,0),axis=ta((0,0,-scale)),radius=axradius)
    zl=label(text='z',align='center',pos=ta((0,0,+scale)),box=False,opacity=0)
    
###############################################################
#EARTH
###############################################################
sphere(pos=ta((0,0,0)),radius=1.0,material=materials.earth)

###############################################################
#TRAJECTORY
###############################################################
orbit=np.loadtxt(BASE_DIR+"solution-elements.dat")
ts=orbit[:,0]
"""
path=[ta((x,y,z)) for x,y,z in zip(orbit[:,1],orbit[:,2],orbit[:,3])]
c=curve(pos=path,radius=axradius,color=color.yellow)
"""

sat=sphere(pos=ta((orbit[0,1],orbit[0,2],orbit[0,3])),
           radius=2*axradius,color=color.yellow,
           make_trail=True,retain=retain)

ev=scene.waitfor('keydown')
ntimes=len(ts)
for it,t in enumerate(ts):
    rate(trate)
    pos=ta((orbit[it,1],orbit[it,2],orbit[it,3]))

    view=setView(it,ntimes,
                 phi0=phi0,dphi=dphi,
                 theta0=theta0,dtheta=dtheta,thetamax=thetamax,anim=True)
    scene.forward=view["forward"]
    #scene.lights=sat.pos
       
    sat.pos=pos
