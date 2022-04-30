#!/usr/bin/env python
# coding: utf-8

# In[1]:



# with general position assumption
# assumes no negative x-axis and no points lie on it
# starting point < ending point assumed for intervals

import numpy as np
# from math import ceil
import math as math
from sympy import Symbol , Eq , solve , sqrt , N
from graphics import *
from operator import itemgetter
from itertools import chain


# In[2]:



def candidateR(P,p,q,a) :
    
    d = P.shape[1]
    n = P.shape[0]
    R = list()
    
    for i in range(n) :
        for k in range(n) :
            for t in range(p+q) :
                for sign_1 in (-1,1) :
                    for sign_2 in (-1,1) :
                        r = Symbol( "r" , positive=True )
                        rhs = P[k][0] - P[i][0] + t*a
                        lhs = sign_1*sqrt( r**2 - P[i][1]**2 ) + sign_2*sqrt( r**2 - P[k][1]**2 ) 
                        R_temp = solve( Eq( lhs , rhs ) ) 
                        R += R_temp
    R = set(R)
    
    intervals_r = dict()
    for r in R :
        for pt in P :
            if(r>pt[1]) :
                interval = [ pt[0] - math.sqrt(r*r-pt[1]*pt[1]) , pt[0] + math.sqrt(r*r-pt[1]*pt[1]) ]
                if r in intervals_r :
                    intervals_r.get(r).append(interval)
                else :
                    intervals_r[r] = [interval]
#         intervals_r[r].sort(key=itemgetter(0,1))
                
#     print(intervals_r)
    return intervals_r
   


# In[3]:



# I = np.array([[1.0,3.5,4.001,5,6,7,8.2,9.3,10.4,11,12.5]])
# I = [[1.0,5],[3.5,6],[4.001,9.3],[7,8.2],[10.4,12.5]]
# I = [[3,10],[2.1,12.0],[6,7.4],[5.4,8.2],[4.001,8.999],[1,14],[11,15],[17.5,100],[13,17],[18,20],[21,23],[24,26],[27,29],[28,30],[29,31],[30,32],[31,33]]
# F = np.sort(I,0)
# I[np.lexsort((I[:,1],I[:,0]))]

def getFaces(R) :
    S = dict()
    for r in R :
        I = sorted(R[r],key=itemgetter(0,1))
        J = sorted(R[r],key=itemgetter(1,0))

        # print( "I and J are :")
        # print(I)
        # print(J)
        # print("---")

        faces = []
        open_intervals = 1
        ci = 1
        cj = 0
        cf = 0
        previous = I[0][0]

        while( ci<len(I) ) :
            if( I[ci][0]<=J[cj][1] ) :
                if( open_intervals>0 ) :
                    faces.append([previous,I[ci][0]])
                    cf += 1
                previous = I[ci][0]
                ci += 1
                open_intervals += 1
            elif( I[ci][0]>J[cj][1] ) :
                faces.append([previous,J[cj][1]])
                previous = J[cj][1]
                cf += 1
                cj += 1
                open_intervals -= 1
        while( cj<len(J) ) :
            faces.append([previous,J[cj][1]])
            previous = J[cj][1]
            cf += 1
            cj += 1
            open_intervals -= 1
        S[r] = faces
        
    return S
    # print(S)


# In[4]:



def candidateP(faces,a) :
    
    M = []
    for i in range(len(faces)) :
        M.append([])
        for j in range(i,len(faces)) :
            steps = math.ceil((faces[j][0]-faces[i][0])/a)
            for k in range(3) :
                if((faces[i][0]+((steps+k)*a))<=faces[j][1]) :
                    M[i].append(faces[i][0]+((steps+k)*a))
    return sorted(set(chain(*M)))


# In[6]:



def colorcode(c) :
    if( c==0 or ( isinstance(c,str) and ( c.lower()=="red" or c.lower()=="r" ) ) ) :
        return 0 
    else :
        return 1

def othercolor(c) :
    if( c==0 or ( isinstance(c,str) and ( c.lower()=="red" or c.lower()=="r" ) ) ) :
        return 1
    else :
        return 0


# In[17]:



def runDP(DP,I,Ir,Ib,u,v,a,C,k,color) :
    
    # base cases
    if( Ir>=len(I) and Ib>=len(I) ) :
        # store value
        return True
    elif( k>=len(C) or ( Ir<len(I) and u<=0 ) or ( Ib<len(I) and v<=0 ) ) :
        # store value
        return False
    elif( DP[Ir][Ib][u][v][k][colorcode(color)] != -1 ) :
        return DP[Ir][Ib][u][v][k][colorcode(color)]
    else :
        i = min(Ir,Ib)
        while( i<len(I) and I[i][0]<=C[k] ) :
            if( I[i][1]<C[k] ) :
                DP[Ir][Ib][u][v][k][colorcode(color)] = False
                return False 
            i = i+1
            
        # find successor 'a' distance away from kth candidate point
        s = k+1
        while( s<len(C) and C[s]-C[k]<a ) :
            s = s+1
            
        if( colorcode(color)==0 ) :
            Ir2 = i
            Ib2 = Ib
            u2 = u-1 if (u>0) else u
            v2 = v
            next_r = k+1
            next_b = s
        else :
            Ir2 = Ir
            Ib2 = i
            u2 = u
            v2 = v-1 if (v>0) else v
            next_r = s
            next_b = k+1
          
        for i in range(next_r,len(C)) :
            if( runDP(DP,I,Ir2,Ib2,u2,v2,a,C,i,0) ) :
                DP[Ir][Ib][u][v][k][colorcode(color)] = True
                return True
        for i in range(next_b,len(C)) :
            if( runDP(DP,I,Ir2,Ib2,u2,v2,a,C,i,1) ) :
                DP[Ir][Ib][u][v][k][colorcode(color)] = True
                return True
        
    DP[Ir][Ib][u][v][k][colorcode(color)] = False
    return False
    


# In[19]:



def initializeDP(R,F,p,q,a) :

    for r in sorted(R.keys()) :
#         if (len(F[r])>=3) : # to see if candidateP works for cases where >=3 intervals exist i.e. meaningful cases
        print(r)
#             print(F[r])
        C = candidateP(F[r],a)
#             DP = np.zeros(len(R[r]),len(R[r]),p+1,q+1,len(C),2) 
        DP = np.full( (len(R[r])+1,len(R[r])+1,p+1,q+1,len(C),2), -1 )
#             print(C)
#             print("\n")
#             break 
        
        sortedI = sorted(R[r],key=itemgetter(0,1)) # do I need to sort again ?
    
        for k in range(len(C)) :
            if( runDP(DP,sortedI,0,0,p,q,a,C,k,"red") or runDP(DP,sortedI,0,0,p,q,a,C,k,"blue") ) :
                print(True)
                return True
        print(False)
    return False


# In[9]:



def drawRadii(P,R) :

    win = GraphWin('CandidateR', 1000, 600) 

    line = Line(Point(100, 300), Point(900, 300))
    line.setWidth(3)
    line.draw(win)
    
    message = Text(Point(300,700), "R")
    for i in range(17) :
        I = Circle(Point(100+(i*50),300),2)
        I.setFill("Black")
        I.draw(win)
    message.draw(win)
    for pt in P :
        point = Circle(Point(100+50*pt[0],300-50*pt[1]), 2)
        point.setFill('black')
        point.draw(win)
    for r in sorted(R.keys()) :
        circles = []
        for pt in P :
            circle = Circle(Point(100+50*pt[0],300-50*pt[1]), 50*N(r))
#             circle.setFill('black')
            circle.draw(win)
            circles.append(circle)
        win.getMouse()
        for c in circles :
            c.undraw()
        win.update()
        
    win.close()
    


# In[10]:



def getPoints() :
    
    # p1 = np.array([20,1])
    # p2 = np.array([21,1])
    p1 = np.array([9,2])
    p2 = np.array([6,1])
    p3 = np.array([1,1])
    P = np.array([p1,p2])
    np.sort(P,0)
    
    return P


# In[11]:


def runBinarySearch(R) :
    
    return 


# In[21]:



def main() :
    
    # add :
    # binary search
    # input
    # graphics

    p = 3
    q = 2
    a = 2
    
    P = getPoints()
    R = candidateR(P,p,q,a)
    drawRadii(P,R)
#     runBinarySearch(R)
    F = getFaces(R)
    print (initializeDP(R,F,p,q,a))
    
if __name__=="__main__":
    main()
    
    


# In[ ]:





# In[ ]:





# In[ ]:




