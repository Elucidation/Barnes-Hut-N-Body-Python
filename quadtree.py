from numpy import *
from random import random,seed
import Tkinter
from PIL import Image, ImageTk, ImageDraw
from time import time

seed(323)

G = 1.0 # Gravitational Constant
THRESHOLD = 1 # Number of bodies after which to subdivide Quad
MAXDEPTH = 10
THETA = 0.5 # Barnes-Hut ratio of accuracy
ETA = 0.5 # Softening factor


NUM_CHECKS = 0 # Counter

class QuadTree:
    """ Class container for for N points stored as a 2D Quadtree """
    root = None
    def __init__(self, bbox, N, theta = THETA):
        self.bbox = bbox
        self.N = N
        self.theta = theta
        self.reset()
    def reset(self):
        self.root = Quad(self.bbox)
    def generate(self):
        # Build up nodes of tree fresh for bodies
        self.reset()
        for x in xrange(self.N):
            # For each body, add to root
            self.root.addBody(x,0)

    def updateSys(self, dt):
        self.calculateBodyAccels()
        global VEL,POS
        VEL += ACC * dt
        POS += VEL * dt
        # for k in xrange(self.N):
        #     VEL[k] += ACC[k] * dt
        #     POS[k] += VEL[k] * dt
    def calculateBodyAccels(self):
        # Update ACC table based on current POS
        for k in xrange(self.N):
            ACC[k] = self.calculateBodyAccel(k)
    def calculateBodyAccel(self, bodI):
        return self.calculateBodyAccelR(bodI, self.root)

    def calculateBodyAccelR(self, bodI, node):
        # Calculate acceleration on body I
        # key difference is that body is ignored in calculations
        # print "calcAcc P:%s,M:%g -> %s" % (p,mass,node.center)
        acc = zeros(2,dtype=float)
        if (node.leaf):
            # print "Leaf"
            # Leaf node, no children
            for k in node.bods:
                if k != bodI: # Skip same body
                    acc += getForce( POS[bodI] ,1.0,POS[k],MASS[k])
        else:
            s = max( node.bbox.sideLength )
            d = node.center - POS[bodI]
            r = sqrt(d.dot(d))
            # print "s/r = %g/%g = %g" % (s,r,s/r)
            if (r > 0 and s/r < self.theta):
                # Far enough to do approximation
                acc += getForce( POS[bodI] ,1.0, node.com, node.mass)
            else:
                # Too close to approximate, recurse down tree
                for k in xrange(4):
                    if node.children[k] != None:
                        acc += self.calculateBodyAccelR(bodI, node.children[k])
        # print "ACC : %s" % acc
        return acc

    def calculateAccel(self, p):
        return self.calculateAccelR(p, self.root)
    def calculateAccelR(self, p, node):
        # Calculate acceleration on point p = [x,y]
        # print "calcAcc P:%s,M:%g -> %s" % (p,mass,node.center)
        acc = zeros(2,dtype=float)
        if (node.leaf):
            # print "Leaf"
            # Leaf node, no children
            for k in node.bods:
                acc += getForce(p,1.0,POS[k],MASS[k])
        else:
            s = max( node.bbox.sideLength )
            d = node.center - p
            r = sqrt(d.dot(d))
            # print "s/r = %g/%g = %g" % (s,r,s/r)
            if (r > 0 and s/r < self.theta):
                # Far enough to do approximation
                acc += getForce(p,1.0, node.com, node.mass)
            else:
                # Too close to approximate, recurse down tree
                for k in xrange(4):
                    if node.children[k] != None:
                        acc += self.calculateAccelR(p, node.children[k])
        # print "ACC : %s" % acc
        return acc

    def __repr__(self):
        return self.__str__()
    def __str__(self):
        return self.root.__str__()
        

def getForce(p1,m1,p2,m2):
    # need to d
    global NUM_CHECKS
    d = p2-p1
    r = sqrt(d.dot(d)) + ETA
    f = array( d * G*m1*m2 / r**3 )
    NUM_CHECKS += 1
    return f


class Quad:
    """ A rectangle of space, contains point bodies """
    def __init__(self,bbox,bod = None,depth=0):
        self.bbox = bbox
        self.center = bbox.center
        self.leaf = True # Whether is a parent or not
        self.depth = depth
        if bod != None: # want to capture 0 int also
            self.setToBody(bod)
            self.N = 1
        else:
            self.bods = []
            self.mass = 0.
            self.com = array([0,0], dtype=float)
            self.N = 0
            
        self.children = [None]*4 # top-left,top-right,bot-left,bot-right

    def addBody(self, idx,depth):
        # print "addBody(B%g,D%g)" % (idx,depth)
        # Recurse if you have a body or you have children
        if len(self.bods) > 0 or not self.leaf:
            # Not empty
            if (depth >= MAXDEPTH):
                self.bods.append(idx)
            else:
                # Subdivide tree
                subBods = [idx] # bodies to add to children
                if len(self.bods) > 0:
                    # if node has no children yet, move own body down to child
                    subBods.append(self.bods[0])
                    self.bods = []

                for bod in subBods:
                    quadIdx = self.getQuadIndex(bod)
                    if self.children[quadIdx]:
                        # child exists, recursively call 
                        self.children[quadIdx].addBody(bod,depth+1)
                    else:
                        # create child
                        subBBox = self.bbox.getSubQuad(quadIdx)
                        self.children[quadIdx] = Quad(subBBox, bod,depth+1)

                self.leaf = False

            # Update CoM
            bodyMass   = MASS[idx]
            # print idx
            # print "(%s * %g + %s * %g) / (%g + %g) = %s " % (self.com, self.mass, POS[idx], bodyMass, self.mass, bodyMass, ((self.com * self.mass + POS[idx] * bodyMass) / (self.mass + bodyMass))) 
            self.com   = (self.com * self.mass + POS[idx] * bodyMass) / (self.mass + bodyMass)
            self.mass += bodyMass
            # self.updateCOM()
        else:
            # Empty Quad, add body directly
            self.setToBody(idx)

        self.N += 1 # Number of bodies incremented
    
    def updateCOM(self):
        if self.leaf:
            self.mass = array(map(lambda x: MASS[x], self.bods)).sum()
            self.com = array(map(lambda x: POS[x]*MASS[x], self.bods)).sum(0) / self.mass
        else:
            self.mass = array(map(lambda child: child.mass if child else 0, self.children)).sum()
            # print array(map(lambda child: child.mass*child.com if child else zeros(2), self.children)).sum(0) / self.mass
            # raise "PAUSE"
            self.com   = array(map(lambda child: child.mass*child.com if child else zeros(2), self.children)).sum(0) / self.mass
        
    def setToBody(self,idx):
        self.bods = [idx]
        self.mass = float( MASS[idx].copy() )
        self.com  = POS[idx].copy()

    def getQuadIndex(self,idx):
        return self.bbox.getQuadIdx(POS[idx])
        
    def __repr__(self):
        return self.__str__()
    def __str__(self):        
        if len(self.bods) > 0:
            bodstring = str(self.bods)
        else:
            bodstring = "PARENT"
        if any(self.children):
            childCount = "C:%g," % sum(map(lambda x: 1 if x else 0, self.children))
            childStr = "\n"
            for x in xrange(4):
                childStr += ("-"*(self.depth+1))+str(x+1)+" "
                if self.children[x]:
                    childStr += str(self.children[x])
                if x < 3: 
                    childStr += "\n"

        else:
            childCount = ""
            childStr = ""
        return "D%g{N:%g,M:%g,%sCOM:%s,B:%s}%s" % (self.depth,self.N,round(self.mass,2),childCount,self.com.round(2),bodstring,childStr)


class BoundingBox:
    def __init__(self,box,dim=2):
        assert(dim*2 == len(box))
        self.box = array(box,dtype=float)
        self.center = array( [(self.box[2]+self.box[0])/2, (self.box[3]+self.box[1])/2] , dtype=float)
        self.dim = dim
        self.sideLength = self.max() - self.min()

    def max(self):
        return self.box[self.dim:]
    def min(self):
        return self.box[:self.dim]
    def inside(self,p):
        # p = [x,y]
        if any(p < self.min()) or any(p > self.max()):
            return False
        else:
            return True
    def getQuadIdx(self,p):
        # y goes up
        # 0 1
        # 2 3
        if p[0] > self.center[0]: # x > mid
            if p[1] > self.center[1]: # y > mid
                return 1
            else:
                return 3
        else:
            if p[1] > self.center[1]: # y > mid
                return 0
            else:
                return 2
    def getSubQuad(self,idx):
        # 0 1
        # 2 3
        # [x  y x2 y2]
        #  0  1  2  3
        b = array([None,None,None,None])
        if idx % 2 == 0:
            # Even #, left half
            b[::2] = [self.box[0], self.center[0]] # x - midx
        else:
            b[::2] = [self.center[0], self.box[2]] # midx - x2
        if idx < 2:
            # Upper half (0 1)
            b[1::2] = [self.center[1], self.box[3]] # midy - y2
        else:
            b[1::2] = [self.box[1], self.center[1]] # y - midy
        return BoundingBox(b,self.dim)

    def __repr__(self):
        return self.__str__()
    def __str__(self):
        return "<<%g,%g,%g,%g>>" % (self.box[0],self.box[1],self.box[2],self.box[3])

# class Body:
#     def __init__(self,info):
#         self.mass = info[0]
#         self.pos = array([info[1],info[2]])
#         self.vel = array([info[3],info[4]])
#     def __repr__(self):
#         return self.__str__()
#     def __str__(self):
#         return "BOD[M:%.2g,P:%.2g,V:%s]" % (self.mass,self.pos,self.vel)

##############################################
##############################################

def drawAccelGrid(drawPtr, gridSize):
    # Draw accel grid
    MAX_ACC = 1
    boundmax = max(BOUNDS.max())
    for x in xrange(gridSize):
        for y in xrange(gridSize):
            p = array([x+0.5,y+0.5],dtype=float) * BOUNDS.sideLength / gridSize
            acc = sys.calculateAccel(p)
            mag = sqrt((acc*acc).sum())
            acc *= (boundmax/gridSize)/mag # Unit value
            # if mag > MAX_ACC: 
            #     mag = MAX_ACC
            # print mag
            # print "%g : %s " % (mag, acc)
            # p = convertPos(p)
            drawLine(drawPtr, p, p+acc, (255,int(255*mag/MAX_ACC),0) )
            # drawPtr.line( (p[0],p[1],p[0]+acc[0],p[1]+acc[1]) , fill = )


def drawBodies(drawPtr):
    for x in xrange(N):
        if BOUNDS.inside(POS[x]):
            # print "Drawing point at %s" % convertPos(pos[x]) 
            # draw.point( convertPos(POS[x]), fill = (255,255,255) )
            p = convertPos(POS[x])
            r = MASS[x]*5
            drawPtr.ellipse( join(p-r,p+r) , fill = (0,0,0) )

def drawVels(drawPtr):
    for x in xrange(N):
        if BOUNDS.inside(POS[x]):
            drawLine(drawPtr, POS[x], POS[x]+VEL[x], (0,0,255) )

def drawAccs(drawPtr):
    for x in xrange(N):
        if BOUNDS.inside(POS[x]):
            drawLine(drawPtr, POS[x], POS[x]+ACC[x], (255,0,0) )

def join(p1,p2):
    return ([p1[0],p1[1],p2[0],p2[1]])

def drawBBOX(drawPtr, node):
    if node.depth == 0:
        drawCross(drawPtr,node.com)
    # if node.leaf and node.depth < 6:
    #     words = str(node.bods)
    #     wordpos = convertPos(  [node.bbox.box[0],node.bbox.box[3]]   ) 
    #     drawPtr.text(wordpos,words,fill=(0,150,0))
    p = join(convertPos(node.bbox.min()),
             convertPos(node.bbox.max()))
    # print "Draw box %s " % p
    drawPtr.rectangle( p , outline=(0,0,255) )

def drawLine(drawPtr, p1, p2, color):
    drawPtr.line( join( convertPos(p1) , convertPos(p2) ), fill=color )

CROSS_SIZE = 20 # pixels
def drawCross(drawPtr, p):
    p = convertPos(p)
    # print "Cross %s" % p
    drawPtr.line( (p[0]-CROSS_SIZE,p[1],p[0]+CROSS_SIZE,p[1]), fill=(0,255,0) )
    drawPtr.line( (p[0],p[1]-CROSS_SIZE,p[0],p[1]+CROSS_SIZE), fill=(0,255,0) )

def drawQuadTree(drawPtr, qtree):
    drawQuadTreeR(drawPtr, qtree.root)
def drawQuadTreeR(drawPtr,node):
    drawBBOX(drawPtr,node)
    for child in node.children:
        if child != None:
            drawQuadTreeR(drawPtr,child)

def convertPos(p):
    # From BOUNDS to IMAGE_SIZE and list format [x,y]
    # p = (x,y) -> [ix, iy]
    c = trunc( ( p - BOUNDS.min()) / (BOUNDS.max()-BOUNDS.min()) * array(IMAGE_SIZE) ).tolist()
    c[1] = IMAGE_SIZE[1] - c[1] # flip y axis
    return c

##############################################
##############################################

N = 100
BOUNDS = BoundingBox([0,0,10,10])

IMAGE_SIZE = (600,400)
GRID_SIZE = 10
CIRCLE_RADIUS = 10 # radius in pixels

# Global variables
# 2D Position
MASS = zeros(N,dtype=float)
POS = zeros((N,2),dtype=float)
VEL = zeros((N,2),dtype=float)
ACC = zeros((N,2),dtype=float)
for i in xrange(N):
    MASS[i] = 1 #random()*3
    POS[i] = BOUNDS.min() + array([random(),random()])*BOUNDS.sideLength
    # VEL[i] = array([random(),random()])*2-1

# MASS = array([1, 2])
# POS = array([ [0.5, 0.5], [9.5, 9.5] ])
# VEL = array([ [-1, -1], [1, 1] ])


sys = QuadTree(BOUNDS, N)
# sys.generate()
# sys.calculateBodyAccels()

# print sys
# for (m,p,v) in zip(MASS,POS,VEL):
#     print "M:%6s   P:%15s   V:%15s" % (m.round(2),p.round(2),v.round(2))
# print "----- SUM/AVG ----------------------------------"
# massTot = sum(MASS)
# posAvg = sum(array([MASS, MASS]).transpose() * POS, 0) / sum(MASS)
# velAvg = sum(array([MASS, MASS]).transpose() * VEL, 0) / sum(MASS)
# print "M:%6s   P:%15s   V:%15s" % (massTot.round(2),posAvg.round(2),velAvg.round(2))


# print sys.calculateAccel([0,0],1)
# print "Number of force calcs: %g" % NUM_CHECKS
DT = 0.1



# im.save("qtree.png", "PNG")
T = 0


def task_update():
    global T,boo
    print "Update t=%g-%g" % (T,T+DT)
    T += DT
    # print "Updating system by dt=%g" % DT
    ts = time()
    sys.generate()
    print "QUADTREE: %gms" % ((time()-ts)*1000)
    # print "Finished generating quad tree"
    # print "Done updating System."

    im = Image.new('RGB', IMAGE_SIZE, (255,255,255))
    draw = ImageDraw.Draw(im)

    ts = time()
    drawQuadTree(draw,sys)
    # drawAccelGrid(draw,GRID_SIZE)
    drawBodies(draw)
    print "DRAW: %gms" % ((time()-ts)*1000)
    # drawVels(draw)
    # drawAccs(draw)
    ts = time()
    sys.updateSys(DT)
    print "SYSUPD: %gms" % ((time()-ts)*1000)
    # print "Number of force calcs: %g" % NUM_CHECKS
    
    ts = time()
    boo = ImageTk.PhotoImage(im)
    # label_image = Tkinter.Label(root, image=boo)
    label_image.config(image=boo)
    label_image.pack()
    print "SHOWING: %gms" % ((time()-ts)*1000)
    root.after(10,task_update)

    print





def button_click_exit_mainloop (event):
    event.widget.quit() # this will cause mainloop to unblock.

root = Tkinter.Tk()
root.bind("<Button>", button_click_exit_mainloop)
label_image = Tkinter.Label(root)

print "Update t=%d-%d" % (T,T+DT)
T += DT
print "Updating system by dt=%g" % DT
sys.generate()
print "Finished generating quad tree"
sys.updateSys(DT)
print "Done updating System."

im = Image.new('RGB', IMAGE_SIZE, (255,255,255))
draw = ImageDraw.Draw(im)

drawQuadTree(draw,sys)
# drawAccelGrid(draw,GRID_SIZE)
drawBodies(draw)
drawVels(draw)
drawAccs(draw)

print "Number of force calcs: %g" % NUM_CHECKS
boo = None
# def update_image():
#     global boo
#     boo = ImageTk.PhotoImage(im)
#     # label_image = Tkinter.Label(root, image=boo)
#     label_image.config(image=boo)
#     label_image.pack()
# root.after(2000,task_update)

# update_image()
task_update()
print "Here"
# root.after(2000,task_update)
root.mainloop()

print "Done."