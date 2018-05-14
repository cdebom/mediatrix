"""
The following functions execute elementary geometrical calculations.
"""

from math import sin, cos ,sqrt, fabs, atan, tan 

def get_line_point_nearest2point(p,theta,c):
    """
    This Function calculates the distance from a line to a point.

    Input:
      - p       [float,float]  : a point p=[x,y]. 	
      - theta          float : angle defined by line angular coefficient.
      - c              float : linear coefficient.
    Output:
     - [float,flot] : The coordinates of the nearest point in the line to the point p.

    """
    if(theta!=2*atan(1)) and (theta!=0.):
        a=tan(theta)
        m=-1./a
        b=p[1]-m*p[0]
        x=(c-b)/(m-a)
        y=a*x+c
        dx=p[0]-x
        dy=p[1]-y
        return x,y
                
    elif (theta==2*atan(1)): # case of the perpendicular bisector is a vertical line a-->inf.
        dx=fabs(p[0]-c)
        return c,p[1]

    elif theta==0.: # case of the perpendicular bisector is a horizontal line a-->0
        dy=fabs(p[1]-c)
        return p[0],c

def two_points_to_line(p1,p2): 
    """
    Function that calculates the line defined by two points  p_i=[x_i,y_i].
    
    Input:
     - p1      [float,float] : first point.
     - p2      [float,float] : second point.
        
    Output:
     - float : angle defined by angular coefficient.
     - float : linear coefficient.
     
    """

    x1=p1[0]
    y1=p1[1]	
    x2=p2[0]
    y2=p2[1]	
    Dy=float(y2-y1)
    Dx=float(x2-x1)
    b1=0
    My=0
    Mx=0
    if (fabs(Dx)>0): 
        a=(Dy)/(Dx)
        theta=atan(a)
        b1=y1-a*x1			
    elif Dx==0:
        a=0
        theta=2*atan(1)
        b1=x1		
            
    return theta,b1	


def get_length_c(points): 
    """
    Function to calculate the length in a path determined by several points.

    Input:
     - points      <list> : list  of p_i points where p_i=x_i+J*y_i.
     
     
    Output: 
     - <float> : the length.  
     
    """
    
    L=0.
    err=0
    for j in range(0,len(points)-1):
        dl=abs(points[j]-points[j+1])
        L=L+dl

    return L	




def get_extrema_2loops( ximg, yimg, ref_position ):
    """
    Function to determine the two extrema of a set of points through a reference point. 

    The first extreme point (P1) is found by searching for the point furthest from the reference 
    point (usually is some definition of center). The second extreme point is the furthest one from 
    point P1.


    Input:
     - ximg        <list> : x coordinates of the image
     - yimg        <list> : y coordinates of the image
     - ref_position <int> : the (Python) position of the reference point in the ximg and yimg lists

    Output: 
     - <int> : position in the lists ximg and yimg of one extreme point
     - <int> : position in the lists ximg and yimg of the other extreme point

    """


    # define the reference point coordinates
    x_ref = ximg[ref_position] 
    y_ref = yimg[ref_position] 

    # find the furthest point from the reference point (1st extreme)
    furthest_pt1 = 0
    distmax = -1
    for i in range(0,len(ximg)): # loops over all points
        if (ximg[i] - x_ref)**2 + (yimg[i] - y_ref)**2 > distmax:
            distmax = (ximg[i] - x_ref)**2 + (yimg[i] - y_ref)**2
            furthest_pt1 = i

    # find the furthest point from the first extreme (2nd extreme)
    furthest_pt2 = 0
    distmax = -1
    for j in range(0,len(ximg)): # loops over all points
        if (ximg[j] - ximg[furthest_pt1])**2 + (yimg[j] - yimg[furthest_pt1])**2 > distmax:
            distmax = (ximg[j] - ximg[furthest_pt1])**2 + (yimg[j] - yimg[furthest_pt1])**2
            furthest_pt2 = j


    return furthest_pt1,furthest_pt2;





def three_points_to_circle(p1, p2, p3):

    """
    Function that calculates circle parameters defined by three points.
    
    Input:
	 - p1 [float,float] : pair (x,y) of coordinates of the 1st point
	 - p2 [float,float] : pair (x,y) of coordinates of the 2nd point
	 - p3 [float,float] : pair (x,y) of coordinates of the 3rd point
        
    Output:
     - float : center x coordinate
     - float : center y coordinate
     - float : circle radius
    """

    theta1,b1=define_perpendicular_bisector(p1,p2)
    theta2,b2=define_perpendicular_bisector(p2,p3)
    xc,yc,Flag=get_intercept_point(theta1,b1,theta2,b2)
    if Flag==1:
        return 0,0,0
    r=((p1[0]-xc)*(p1[0]-xc))+((p1[1]-yc)*(p1[1]-yc))
    r=sqrt(r)
    return xc,yc,r





def get_intercept_point(theta,b_t,phi,b_p):

    """
    Function that calculates the interception point from two straight lines.
    
    Input:
     - theta      float : angle defined by angular coefficient from first line.
     - b_t        float : linear coeficient from the first line.
     - phi        float : angle defined by angular coefficient from second line.
     - b_p        float : linear coeficient from the second line.
   
    Output:
     - float : intercept x coordinate.
     - float : intercept y coordinate.
     - int   : a Flag. if Flag= 1 the lines are parallels, otherwise Flag=0,
    """

    a_t=tan(theta)
    a_p=tan(phi)
    Flag=0
    if(fabs(a_t-a_p)<0.001):
        print "parallels found"
        Flag=1
        xi=0
        yi=0
    elif (theta==2*atan(1) and phi==0.):
        xi=b_t
        yi=b_p
    elif (theta==0. and phi==2*atan(1)):
        xi=b_p
        yi=b_t
    elif theta==2*atan(1):
        xi=b_t
        yi=a_p*xi+b_p
    elif phi==2*atan(1):
        xi=b_p
        yi=a_t*xi+b_t
    elif theta==0:
        yi=b_t
        xi=(yi-b_p)/a_p
    elif phi==0:   
        yi=b_p
        xi=(yi-b_t)/a_t
    else:
        xi=(b_p-b_t)/(a_t-a_p)
        yi=a_t*xi+b_t
    return xi,yi,Flag


def define_perpendicular_bisector(p1,p2): 
    """
    Function that calculates the perpendicular bisector line between two points 
    p_i=[x_i,y_i].
    
    Input:
     - p1      [float,float] : first point.
     - p2      [float,float] : second point.
        
    Output:
    - float : angle defined by angular coefficient.
    - float : linear coefficient.
     
    """

    x1=p1[0]
    y1=p1[1]	
    x2=p2[0]
    y2=p2[1]	
    Dy=float(y2-y1)
    Dx=float(x2-x1)
    b1=0
    My=0
    Mx=0
    if (fabs(Dy)>0): 
        m=-(Dx)/(Dy)
        theta=atan(m)
        My=float(y2+y1)/2.
        Mx=float(x2+x1)/2.
        
        b1=My-m*Mx			
    elif Dy==0:
        m=0
        theta=2*atan(1)
        b1=float(x2+x1)/2.		
            
    return theta,b1	

def two_points_to_line(p1,p2): 
    """
    Function that calculates the line defined by two points  p_i=[x_i,y_i].
    
    Input:
     - p1      [float,float] : first point.
     - p2      [float,float] : second point.
        
    Output:
     - float : angle defined by angular coefficient.
     - float : linear coefficient.
     
    """

    x1=p1[0]
    y1=p1[1]	
    x2=p2[0]
    y2=p2[1]	
    Dy=float(y2-y1)
    Dx=float(x2-x1)
    b1=0
    My=0
    Mx=0
    if (fabs(Dx)>0): 
        a=(Dy)/(Dx)
        theta=atan(a)
        b1=y1-a*x1			
    elif Dx==0:
        a=0
        theta=2*atan(1)
        b1=x1		
            
    return theta,b1	





def get_distance_from_line_to_point(p,theta,c):
    """
    This Function calculates the distance from a line to a point.

    Input:
      - p       [float,float]  : a point p=[x,y]. 	
      - theta          float : angle defined by line angular coefficient.
      - c              float : linear coefficient.
    Output:
     - float : The distance from the line to the point p.

    """
    if(theta!=2*atan(1)) and (theta!=0.):
        a=tan(theta)
        m=-1./a
        b=p[1]-m*p[0]
        x=(c-b)/(m-a)
        y=a*x+c
        dx=p[0]-x
        dy=p[1]-y
        Distance=sqrt(dx*dx + dy*dy)
                
    elif (theta==2*atan(1)): # case of the perpendicular bisector is a vertical line a-->inf.
        dx=fabs(p[0]-c)
        Distance=dx

    elif theta==0.: # case of the perpendicular bisector is a horizontal line a-->0
        dy=fabs(p[1]-c)
        Distance=dy
    return Distance 


def get_line_point_nearest2point(p,theta,c):
    """
    This Function calculates the distance from a line to a point.

    Input:
      - p       [float,float]  : a point p=[x,y]. 	
      - theta          float : angle defined by line angular coefficient.
      - c              float : linear coefficient.
    Output:
     - [float,flot] : The coordinates of the nearest point in the line to the point p.

    """
    if(theta!=2*atan(1)) and (theta!=0.):
        a=tan(theta)
        m=-1./a
        b=p[1]-m*p[0]
        x=(c-b)/(m-a)
        y=a*x+c
        dx=p[0]-x
        dy=p[1]-y
        return x,y
                
    elif (theta==2*atan(1)): # case of the perpendicular bisector is a vertical line a-->inf.
        dx=fabs(p[0]-c)
        return c,p[1]

    elif theta==0.: # case of the perpendicular bisector is a horizontal line a-->0
        dy=fabs(p[1]-c)
        return p[0],c
    


def width_ellipse(l,area):
    """
    Function to calculate wiidth of an ellipse of given area and length l.

    Input:
     - l          float : the ellipse length.
     - area       float : the ellipse area.
     
     
    Output: 
     - float : the width.  
     
    """

    w=area/(atan(1)*l)
    return w

def length_from_connected_dots(points): 
    """
    Function to calculate the sum of distance from a previous point to the next on 
    a list of points.

    Input:
     - points  [[float,float],...] : list  of p_i points where p_i=(x_i,y_i).
     
     
    Output: 
     - float : the length.  
     
    """
    x=[]
    y=[]
    for i in range(0,len(points)):
        x.append(points[i][0])
        y.append(points[i][1])
    l=0
    err=0
    for j in range(0,len(X)-1):
        dx=x[j]-x[j+1]
        dy=x[j]-y[j+1]
        dx=dx*dx
        dy=dy*dy
        dl=sqrt(dx+dy)
        l=l+dl

    return float(l)	



