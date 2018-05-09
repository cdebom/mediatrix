
import numpy as np
from el_geometry import get_extrema_2loops,get_length_c, three_points_to_circle, define_perpendicular_bisector, get_distance_from_line_to_point
from math import atan
import aplpy
from scipy.optimize import fmin
import itertools as it

def filamentation(image, method="medium",alpha=1,near_distance=(sqrt(2)/2), max_level=1000, use_extremes=True):
    """
    Function to perform the mediatrix decomposition method on a given object. 

    Input:
     - image_name   <str> : the image file name.
     - image_dir   <str> : the image directory. If it is on the same directory, directory=''.
     - method   <string> : possible values are 'medium'  or 'brightest'.
     - alpha      <float> : the factor alpha=l_i/w to stop the bisection.
     - near_distance      <float> : the distance to consider a point near to the perpendicular bisector.
     
    Output:
     - <dic> :  Dictionary structure. Each list item is a dictionary with information of corresponding to a mediatrix vector. The keys are 'theta' for the angle with x axis, 'linear_coefficient' for the linear coefficient from the line in the vector direction, 'origin' the point (x,y) of the vector origin, 'end' the point (x,y) of the vector, 'modulus' for vector modulus. The first item from the list has two extra keys 'id' wich contains the image_name and 'center' that keeps the objet center defined by the first mediatrix point in the first mediatrix level.
            
    """
    #image,hdr = getdata(image_dir+image_name, header = True )
    pixels=np.where(image>0)
    E1,E2=get_extrema_2loops(pixels[0], pixels[1], 0 )
    Area=len(pixels[1])
    p1=pixels[0][E1]+pixels[1][E1]*1j # the extreme points p_1 and p_2
    p2=pixels[0][E2]+pixels[1][E2]*1j
    keydots=[p1,p2]
    key_center=[p1,p2]
    key_center=find_keydots_c(p1,p2,pixels,image,key_center,Area, method=method,alpha=0,near_distance=near_distance,max_level=1,level=0)
    keydots=find_keydots_c(p1,p2,pixels,image,keydots,Area, method=method,alpha=alpha,near_distance=near_distance,max_level=max_level,level=0)
    #print "key center inside mediatrix"
    #print key_center
    #print "regular key center inside mediatrix"
    #print keydots
    #print keydots
    L=get_length_c(keydots)
    if use_extremes==False:
        if len(keydots)>3:
            keydots.remove(p1)
            keydots.remove(p2)
        else:
            print "Impossible to reject extremes, too few points" 
    mediatrix_vectors=find_mediatrix_vectors_c(keydots)
    #mediatrix_vectors['id']=image_name
    #medium=int(float(len(keydots))/2)
    if len(key_center)>2:
        mediatrix_vectors['center']=key_center[1]
    else:
        medium=int(float(len(keydots))/2)
        mediatrix_vectors['center']=keydots[medium]
    W=(len(pixels[0]))/(atan(1)*L)
    mediatrix_vectors['L/W']=L/W
    mediatrix_vectors['L']=L
    #x=[pixels[0][E1],mediatrix_vectors['center'].real,pixels[0][E2]]
    #y=[pixels[1][E1],mediatrix_vectors['center'].imag,pixels[1][E2]]
    p1_vec=[pixels[0][E1],pixels[1][E1]] # the extreme points p_1 and p_2
    p2_vec=[pixels[0][E2],pixels[1][E2]]
    p3_vec=[mediatrix_vectors['center'].real,mediatrix_vectors['center'].imag]
    x_c,y_c,r=three_points_to_circle(p1_vec,p3_vec,p2_vec)
    circle_center=x_c+y_c*1j
    mediatrix_vectors['circle_params']=[circle_center,p1,p2]

    return mediatrix_vectors



def find_keydots_c(p1,p2,image_pixels,image,keydots,area, method="medium",alpha=1,near_distance=(sqrt(2)/2),max_level=1000,level=0):
    """
    Function to calculate the keydots Points in Mediatrix Decomposition.
    
    Input:
     - p1      <array> : coordinates (x,y) of the first extreme point.
     - p2      <array> : coordinates (x,y) of the second extreme point.
     - image_pixels   <list> : list of points coordidates fitting the object.
     - image   <array> : the image matrix.
     - keydots  <array> : list with the two p_i extreme points and p_i=[p_i_x,p_i_y].
     - area  <array> : the object area.
     - method   <string> : possible values are 'medium' or 'brightest'.
     - alpha      <float> : the factor alpha=l_i/w.
     - near_distance      <float> : the distance to consider a point near to the perpendicular bisector.
     
    Output:
     - <array> : list  with the keydots points.  
     
    """
    level=level+1
    if (p1 in keydots) and (p2 in keydots):
        index1=keydots.index(p1)
        index2=keydots.index(p2)
        if  index1>index2:
            indexNext=index1
        else:
            indexNext=index2
    elif (p1 in keydots):
        indexNext=keydots.index(p1)
    elif (p1 in keydots):
        indexNext=keydots.index(p2)
    else:
        return keydots
	
    
    dl=abs(p1-p2)
    #print "complex"
    #print dl
    L=get_length_c(keydots)
    W=width_ellipse(L,area)

    if dl>(alpha*W) and len(keydots)<100 and level<=max_level:
        p1_r=[p1.real,p1.imag]
        p2_r=[p2.real,p2.imag]
        coefficients=define_perpendicular_bisector(p1_r,p2_r)
        p3,p3Flag=choose_near_point_c_vec(coefficients[0],coefficients[1],image_pixels,image,method,near_distance)
        #p3,p3Flag=choose_near_point_c(coefficients[0],coefficients[1],image_pixels,image,method,near_distance)
        
        #print "p3Flag"
        #print p3Flag
        if (p3Flag==0):
            
            if level==1:	
                #print "this is the chosen point"
                #print p3
                px=p3.real
                py=p3.imag
                test_points2=np.array([[px-1,py+1],[px,py+1],[px+1,py+1],[px-1,py],[px,py],[px+1,py],[px-1,py-1],[px,py-1],[px+1,py-1]])
                Ds=np.apply_along_axis(get_distance_from_line_to_point,1,test_points2,theta=coefficients[0],c=coefficients[1])
                Ds=np.array(Ds)

                """
                print "this is the neighbours and its intensities"
                print "this is the near distance "+str(near_distance)
                print str(p3.real-1)+' '+str(p3.imag+1)+' '+str(round(image[p3.real-1][p3.imag+1],4))+" "+str(Ds[0])+" up-left" 
                print str(p3.real)+' '+str(p3.imag+1)+' '+str(round(image[p3.real][p3.imag+1],4))+" "+str(Ds[1])+" up-center"
                print str(p3.real+1)+' '+str(p3.imag+1)+' '+str(round(image[p3.real+1][p3.imag+1],4))+" "+str(Ds[2])+" up-right"
                print str(p3.real-1)+' '+str(p3.imag)+' '+str(round(image[p3.real-1][p3.imag],4))+" "+str(Ds[3])+" mid-left"   
                print str(p3.real)+' '+str(p3.imag)+' '+str(round(image[p3.real][p3.imag],4))+" "+str(Ds[4])+" mid"
                print str(p3.real+1)+' '+str(p3.imag)+' '+str(round(image[p3.real+1][p3.imag],4))+" "+str(Ds[5])+" mid-right"
                print str(p3.real-1)+' '+str(p3.imag-1)+' '+str(round(image[p3.real-1][p3.imag-1],4))+" "+str(Ds[6])+" down-left"
                print str(p3.real)+' '+str(p3.imag-1)+' '+str(round(image[p3.real][p3.imag-1],4))+" "+str(Ds[7])+" down-center"
                print str(p3.real+1)+' '+str(p3.imag-1)+' '+str(round(image[p3.real+1][p3.imag-1],4))+" "+str(Ds[8])+" down-right"
                """

                                       
	
            if (not(p3 in keydots)):
                keydots.insert(indexNext,p3)
                keydots=find_keydots_c(p1,p3,image_pixels,image,keydots,area, method,alpha,near_distance,max_level,level)
                keydots=find_keydots_c(p3,p2,image_pixels,image,keydots,area, method,alpha,near_distance,max_level,level)
        else:
            pmed=(p1+p2)/2.
            if p1 in keydots: 
                keydots=find_keydots_c(p1,pmed,image_pixels,image,keydots,area, method,alpha,near_distance,max_level,level)
            if p2 in keydots:
                keydots=find_keydots_c(pmed,p2,image_pixels,image,keydots,area, method,alpha,near_distance, max_level,level)

    return keydots



def find_mediatrix_vectors_c(points): 
    """
    From a given set of points, this function returns the mediatrix decomposition vector between those points.
  
    Input:
     - points      <list> : list  of p_i points where p_i=(x_i,y_i).
     
    Output:
     - <list> : a list of dictionary structure. Each list item is a dictionary with information of corresponding to a mediatrix vector. The keys are 'theta' for the angle with x axis, 'linear_coefficient' for the linear coefficient from the line in the vector direction, 'origin' the point (x,y) of the vector origin, 'end' the point (x,y) of the vector, 'modulus' for vector modulus.  
     
    """
    mediatrix_vectors= {'origin': [] , 'end': [], }
    vectors=[]
    p1_r=[points[0].real,points[0].imag]
    p2_r=[points[len(points)-1].real,points[len(points)-1].imag]
    theta_ext,c_ext=two_points_to_line(p1_r,p2_r)
    for t in range(0,len(points)-1):
        p1_r=[points[t].real,points[t].imag]
        p2_r=[points[t+1].real,points[t+1].imag]
        #print p1_r, p2_r
        coefficients=define_perpendicular_bisector(p1_r,p2_r)
        
        origin=(points[t]+points[t+1])/2.
        modulus=abs(points[t]-points[t+1])
           

	if(coefficients[0]!=2*atan(1)):
            end1=origin + modulus*(cos(coefficients[0]))+modulus*(sin(coefficients[0]))*1j
            end2=origin - modulus*(cos(coefficients[0]))-modulus*(sin(coefficients[0]))*1j
            
            
        else:
            end1=origin+modulus*1j
            end2=origin-modulus*1j
        end1_r=[end1.real,end1.imag]
        end2_r=[end2.real,end2.imag]
        Dend1=get_distance_from_line_to_point(end1_r,theta_ext,c_ext)
        Dend2=get_distance_from_line_to_point(end2_r,theta_ext,c_ext)
        if Dend1<Dend2:
	    end=end1
        else:
            end=end2

        #mediatrix_vector = {'theta': coefficients[0], 'linear_coefficient': coefficients[1], 'origin': origin, 'end': end, 'modulus': modulus }
        mediatrix_vectors['origin'].append(origin) 
        mediatrix_vectors['end'].append(end)
        #vectors.append(mediatrix_vector)
                
        
    return mediatrix_vectors

def choose_near_point_c_vec(theta,c,object_pixels,image,method,near_distance): 
    FlagErr=0
    chosenX=0
    chosenY=0    
    test_points=np.array(np.transpose(object_pixels))
    
    #print "this are the Test points"
    #print test_points
    #D=get_distance_from_line_to_point(pixel,theta,c)
    Ds=np.apply_along_axis(get_distance_from_line_to_point,1,test_points,theta=theta,c=c)
    Ds=np.array(Ds)
    #print "this are D min e max"
    #print np.max(Ds)
    #print np.median(Ds)
    #print np.min(Ds)
    near_index=np.where(Ds<=near_distance)
    #print near_index[0]
    #print len(near_index[0])
    if len(near_index[0])==0:
        FlagErr=1
    else:    
        Ds=Ds[near_index]
        near_points=test_points[near_index]
        near_points=np.transpose(near_points)
        object_pixels_near=tuple((near_points[0],near_points[1]))
    
    

    if method=='brightest' and FlagErr==0:
        brights=image[object_pixels_near]
        index_chosen_arr=np.where(brights==np.max(brights))
        index_chosen=index_chosen_arr[0][0]
        chosenX=float(object_pixels_near[0][index_chosen])
        chosenY=float(object_pixels_near[1][index_chosen])    
    elif method=='brightestinline' and FlagErr==0:
        brights=image[object_pixels_near]
        index_chosen_arr=np.where(brights==np.max(brights))
        index_chosen=index_chosen_arr[0][0]
        X_aux=float(object_pixels_near[0][index_chosen])
        Y_aux=float(object_pixels_near[1][index_chosen]) 
        chosenX,chosenY=get_line_point_nearest2point([X_aux,Y_aux],theta=theta,c=c)
    elif method=='medium' and FlagErr==0:
        i,j=get_extrema_2loops( near_points[0], near_points[1], 0 )
        chosenX=float(near_points[0][i]+near_points[0][j])/2.
        chosenY=float(near_points[1][i]+near_points[1][j])/2.
    
    else:
        FlagErr=1
            
    return chosenX+chosenY*1j,FlagErr



def Evaluate_S_Statistic_on_matrix_c(mediatrix_data,obj_stamp, sigma_out=True,sigma=1,sigma_pre=0.5,Area_out=True):

    """
    Function to calculate the S estatistic measurements.
    
    Input:
     - mediatrix_data <list> : a list of dictionary structure. Each list item is a dictionary with information corresponding to a mediatrix vector. The keys are 'theta' for the angle with x axis, 'linear_coefficient' for the linear coefficient from the line in the vector direction, 'origin' the point (x,y) of the vector origin, 'end' the point (x,y) of the vector, 'modulus' for vector modulus. The first item from the list has an extra key 'id' wich contains the image file name. It is the output from Mediatrix_Decomposition.

     
    Output:
     - <dic> : dictionary with the S statistic measurements. The keys are  
     
    """
    guess=[mediatrix_data['center'].real,mediatrix_data['center'].imag]
    L=0
    theta=[]
    linear=[]
    L=mediatrix_data['L']
    #print " E o L eh"
    #print L
    for i in range(0,len(mediatrix_data['origin'])):
       origin=mediatrix_data['origin'][i]
       end=mediatrix_data['end'][i]
       delta=end-origin
       if delta.real!=0:
           a=float((delta.imag ))/delta.real
           theta.append(atan(a))
           b=end.imag-a*(end.real)
           linear.append(b)
       else:
           theta.append(2*atan(1))
           linear.append(end.imag-origin.imag)
    minM_r=fmin(M_function,guess,args=(theta,linear),maxiter=1000, disp=0)
    minM_val=M_function(minM_r,theta,linear)
    minM=minM_r[0]+minM_r[1]*1j
    R=abs(mediatrix_data['center']-minM)
    #print "esse e o ponto em que m eh minima"
    #print minM
    circle_center=mediatrix_data['circle_params'][0]
    #print "Este e o centro do circulo"
    #print circle_center
    center_comparison=abs(circle_center-minM)
    alpha=Find_angle_from_circle_section_c(mediatrix_data['circle_params'][0],mediatrix_data['circle_params'][1],mediatrix_data['circle_params'][2])
    S_output={'init': 0}
    try: 
        S_output['MinM_norm']=minM_val/(L*L)
    except:
        print " Impossible to define a lenght "#for "+str(mediatrix_data['id'])
        S_output['MinM_norm']=-1
    if R!=0:
        S_output['L/R']=L/R
    else:
        print "Impossible to define a radius and curvature"# for "+str(mediatrix_data['id'])
        S_output['L/R']=-1
    if R!=0 and alpha!=0:
        S_output['curvature_comparison']=(L/R)/alpha
    else:
        print "Impossible to compare curvature from circle section and mediatrix S"# for "+str(mediatrix_data['id'])
        S_output['curvature_comparison']=-1
    try:
        S_output['center_comparison']=center_comparison/L
    except:
        print "Impossible to compare center from circle method and mediatrix S"# for "+str(mediatrix_data['id'])
        S_output['curvature_comparison']=-1
    if sigma_out==True:
        
	search=np.array([])
        
        rc_center=minM
 
        
        sigma_points=[]
        step=sigma_pre
        
        lim=step
        sigma_points=find_rc_region_c(sigma_points,rc_center,step,lim=lim,theta=theta,linear=linear,minM=minM_val,sigma=sigma) 
        
        
        sigma_points=np.array(sigma_points)
        sigma_points=np.unique(np.round(sigma_points.real,0)+np.round(sigma_points.imag,0)*1j)


        E1,E2=get_extrema_2loops(sigma_points.real, sigma_points.imag, 0 )
        Extreme=[sigma_points[E1],sigma_points[E2]]
        sigma_lenght=get_length_c(Extreme)
        try:
            S_output['sigma_lenght']=(sigma_lenght)/L
        except:
            print "Impossible to define sigma_lenght for "+str(mediatrix_data['id'])
        sigma_minor_axis=len(sigma_points)/(4*atan(1)*sigma_lenght)
        #verificar essa eq.
        #print "Area de sigma, complex"
        #print sigma_points
        try:
            S_output['sigma_exc']=(sigma_lenght)/sigma_minor_axis
        except:
            print "Impossible to define sigma_exc for "+str(mediatrix_data['id'])
        if Area_out==True:
            
            pixels=np.where(obj_stamp>0)
            pixels_c=np.array(pixels[0]+pixels[1]*1j)
            intersection_points=np.intersect1d(pixels_c,sigma_points)
            S_output['intersection']=len(intersection_points)
                
    return S_output


def M_function(p,theta,b):
    M=0
    for i in range(0,len(theta)):
        #print "este e o theta real"+str(theta[i])
        aux=get_distance_from_line_to_point(p,theta[i],b[i])
        M=M+(aux*aux)
        #print aux
    M=M/float(len(theta))
    return M


def M_function_c(p,theta,b):
    M=0
       
    for i in range(0,len(theta)):
        D=0
        #print "este e o theta imag"+str(theta[i])
        if(theta[i]!=2*atan(1)) and (theta[i]!=0.):
            
            a=tan(theta[i])
            m=-1./a
            c=p.imag-m*p.real
            x=(b[i]-c)/(m-a)
            y=a*x+b[i]
            p_line=x+y*1j
            D=abs(p-p_line)
        elif (theta[i]==2*atan(1)): # case of the perpendicular bisector is a vertical line a-->inf.
            D=abs(p.real-b[i])
        elif theta[i]==0.: # case of the perpendicular bisector is a horizontal line a-->0
           D=abs(p.imag-b[i])
        #print D
        M+=(D*D)
    M=M/float(len(theta))
    return M


def Find_angle_from_circle_section_c(circle_center,circle_point1,circle_point2):

    m_point=(circle_point1+circle_point2)/2
    delta_adjacent=abs(m_point-circle_center)
    delta_opposite=abs(circle_point1-circle_center)
    alpha=atan(delta_opposite/delta_adjacent)
    alpha=2*alpha 
    return alpha



def find_rc_region_c(region,center,step,lim,theta,linear,minM,sigma):
    
    nex_level=0
    ini_region=len(region)
    X_search=np.arange(round(center.real,2)-lim, round(center.real,2)+lim,step)
    Y_search=np.arange(round(center.imag,2)-lim, round(center.imag,2)+lim,step)
    
    search=np.array([])
    #search_len=len(X_search)*len(Y_search)
   
    #for i in range(0,len(X_search)):
    #    piece=np.repeat(X_search[i],len(Y_search))+Y_search*1j
    #    search=np.concatenate((piece,search),axis=1)
    len_search=0
    piece=np.repeat(X_search[0],len(Y_search))+Y_search*1j
    search=np.concatenate((piece,search),axis=0)
    len_search+=len(piece)

    piece=np.repeat(X_search[-1],len(Y_search))+Y_search*1j
    search=np.concatenate((piece,search),axis=0)
    len_search+=len(piece)
    
    piece=X_search+np.repeat(Y_search[0],len(X_search))*1j
    search=np.concatenate((piece,search),axis=0)
    len_search+=len(piece)    

    piece=X_search+np.repeat(Y_search[-1],len(X_search))*1j
    search=np.concatenate((piece,search),axis=0)
    len_search+=len(piece)

    #search_len=len(X_search)*4
    search=np.reshape(search,len_search)
    #search=np.unique(search)     
    m_values_c=np.apply_along_axis(M_function_c,0,search,theta,linear)
    m_values_c=m_values_c- np.repeat(minM,len(search))
    for i in range(0,len(search)):
            if m_values_c[i]<=sigma:
               region.append(search[i])
    #print "esse e o tamanho do search"
    #print len(search)
    if len(region)>ini_region and len(region) < 2000:
        #print "entrando na recursividade"
        region=find_rc_region_c(region,center,step,lim+step,theta,linear,minM,sigma)
    return region






def plot_S_Statistic_apl(image_name,_id='',keydots=False,circle=False,rc=True, save=True,out_image='', args={},use_extremes=True):
    # FIXME to work in world coord==true
    """
    Function to make plot of S estatistic interesting measurements.
    
    Input:
     - mediatrix_data <dic> : the output from mediatrix decomposition.

     
    Output:
     - <bool> :   
     
    """
    opt={'increase': 2, 'relative_increase': True,'connected': False,'object_centered':True, 'out_type':'cutout', 'pmin':0.25 , 'pmax':99.75 , 'invert':True ,'out_title': 'Mediatrix Method', 'keys_color': "r" ,'alpha': 1 ,'max_level': 1000, 'near_distance': sqrt(2)/2, 'max_level': 1000, 'method':"brightest",'sigma':1,'sigma_pre':0.5, 'rc_color': 'm', 'world_coord': 'False'}
    opt.update(args)
     
    if out_image=='':
        out_image=image_name.replace(".fits","")+"_mediatrixS_plot.png"
 

    

    mediatrix_plot,mediatrix_data,image_ps=plot_mediatrixapl(image_name,_id=_id, keydots=keydots,circle=circle, save=False, args=opt,use_extremes=use_extremes)  
    
    

    guess=[mediatrix_data['center'].real,mediatrix_data['center'].imag]
    Length=0
    theta=[]
    linear=[]
    for i in range(0,len(mediatrix_data['origin'])):
       origin_x=mediatrix_data['origin'][i].real
       origin_y=mediatrix_data['origin'][i].imag
       end_x=mediatrix_data['end'][i].real
       end_y=mediatrix_data['end'][i].imag
       Length_aux=(origin_x - end_x)**2 + (origin_y - end_y)**2
       Length=Length+ sqrt(Length_aux)
       delta_x=float((end_x-origin_x ))
       if delta_x!=0:
           a=float((end_y-origin_y ))/delta_x
           theta.append(atan(a))
           b=end_y-a*(end_x)
           linear.append(b)
       else:
           theta.append(2*atan(1))
           linear.append(end_y-origin_y)
        
    MinM=fmin(M_function,guess,args=(theta,linear),maxiter=1000, disp=0)
    MinM_val=M_function(MinM,theta,linear)
    MinM_plot_X=[MinM[0],MinM[0]]
    MinM_plot_Y=[MinM[1],MinM[1]]
    



    mediatrix_plot.show_markers(MinM_plot_Y,MinM_plot_X,c='red',marker='.',zorder=2000)
    #print 'min m'
    #print MinM_plot_X
    #print MinM_plot_Y
    lim=round((float(Length)/2.)+5,2)
    sigma_X=[]
    sigma_Y=[]
    ps=[]
    Ms=[]
    pixels=where(image_ps>0)
    if rc==True:
        X_search=arange(round(MinM[0],2)-lim, round(MinM[0],2)+lim,opt['sigma_pre']*1)
        Y_search=arange(round(MinM[1],2)-lim, round(MinM[1],2)+lim,opt['sigma_pre']*1)


        search_aux=[X_search,Y_search]
        test_points=np.array(list(it.product(*search_aux)))

        #print "this are the Test points"
        #print test_points
    
        sigma=np.apply_along_axis(M_function,1,test_points, theta=theta,b=linear)
        sigma=np.array(sigma)
        sigma=sigma-MinM_val
        #print "this are the M_func values"
        #print sigma
        cr_index=np.where(sigma<1)
        sigma=sigma[cr_index]
        cr_points=test_points[cr_index]
    
        #print "this are the CR points"
        #print cr_points
        for k in range(0,len(cr_points)):
            p=cr_points[k]
            mediatrix_plot.show_markers(p[1],p[0],c=opt['rc_color'],marker=".",linewidths='0', zorder=999)
              

    if save==True:
        #xmin, xmax = xlim()
        #ymin, ymax = ylim()
        #mediatrix_plot.set_xlim(xmin-1*Length,xmax+1*Length)
        #mediatrix_plot.set_ylim(ymin-1*Length,ymax+1*Length)
        savefig(out_image)
        #mediatrix_plot.clear()
        return True
    else:
        #print " os pixels"
        #print pixels[0]
        #print pixels[1] 
        size0=abs(max(pixels[0])-min(pixels[0]))
        size1=abs(max(pixels[1])-min(pixels[1])) 
        return mediatrix_plot, [mediatrix_data['center'].real,mediatrix_data['center'].imag], [MinM[0],MinM[1]],[size0,size1]

def plot_mediatrixapl(image_name,_id='', keydots=False,circle=True, save=True, out_image='', args={},use_extremes=True):
    """
    Make a plot presenting the object, keydots and mediatrix vectors. If the input image name is
    not a postage stamp the code will read from the sextractor segmentation image the object
    position with given id and pixels intensity from sextractor objects image. The   function
    assumes that the segmentation and objects images names are *original_image_name*_seg.fits and
    *original_image_name*_obj.fits respectively.

    Input:
    - mediatrix_data <list> : the output from mediatrix_decomposition_on_matrix.
    - image_dir   <str> : the image directory. If it is on the same directory, directory=''.
    - keydots   <bool> : 'True' if you want to display the keydots and 'False' if you do not. 
    - colors   <dic> : set the plot colors. The possible keys are 'object', 'vector' and 'keydots'.
    - type <str> : cutout or object       
    Output:
     <bool>
         
    """
    opt={'increase': 2, 'relative_increase': True,'connected': False,'object_centered':True, 'out_type':'cutout', 'vmin':0 , 'invert':True ,'out_title': 'Mediatrix Decomposition', 'keys_color': "r" ,'alpha': 1 ,'max_level': 1000, 'near_distance': sqrt(2)/2, 'max_level': 1000, 'method':"brightest", 'world_coord': 'False', 'circle_radius': False, 'bkg_img': ''}
    opt.update(args)
    
    
   

    #image_seg_hdu=fits.open(image_segname)
    #image_obj_hdu=fits.open(image_objname)
    


    if opt['out_type']=='cutout':
        opt['object_centered']=False
    #else:
    #    opt['object_centered']=True
    type_arg=type(image_name) is str
    if type_arg:
        if out_image=='':
            out_image=image_name.replace(".fits","")+"_mediatrix_plot.png"
        if _id=='':
            image_ps,hdr=fits.getdata(image_name, header = True )
        else:
            image_segname=image_name.replace(".fits","")+"_seg.fits"
            image_objname=image_name.replace(".fits","")+"_obj.fits"
            image_seg,hdr = fits.getdata(image_segname, header = True )
            image_obj = fits.getdata(image_objname, header = False )
            image_ps,hdr=imagetools.segstamp(segimg=image_seg, objID=_id, objimg=image_obj, hdr=hdr, increase=opt['increase'], relative_increase=opt['relative_increase'], connected=opt['connected'], obj_centered=opt['object_centered'])
            del image_seg
            del image_obj
            gc.collect() 
            image_name_ps=image_name.replace(".fits","")+"_ps.fits"
            if opt['world_coord']=='True':
                fits.writeto(image_name_ps,image_ps.astype(float),header=hdr)
            
    else:
        image_ps=image_name.copy()
        if out_image=='':
            time_id=time.time()
            out_image=str(time_id)+"_mediatrix_plot.png"
 
    #image_ps,hdr=imcp.segstamp(segimg=image_seg, objID=_ids[i], objimg=image_obj, hdr=hdr, increase=2, relative_increase=True, connected=False, obj_centered=True)

    mediatrix_data=filamentation(image_ps, method=opt['method'],alpha=opt['alpha'],near_distance=opt['near_distance'],max_level=opt['max_level'], use_extremes=use_extremes) 
    
    if opt['bkg_img']!='':
        img,hdr=fits.getdata(opt['bkg_img'], header = True ) 
    else:
        img=image_ps.copy()

    #print "depois"
    #for j in range(0,len(img[0])):
    #    print "\n"
    #    for i in range(0,len(img[1])):
    #        print img[j][i]
    IDtime=str(time.time())
    #fits.writeto(ID+".test.fits",img.astype(float),header=None)
    pixels=np.where(image_ps>0)
    if opt['world_coord']=='True':
        if opt['out_type']=='cutout': 
            FitsPlot = aplpy.FITSFigure(image_name)
        else:
            FitsPlot = aplpy.FITSFigure(image_name_ps) #FIXME
            os.system(" rm "+image_name_ps)
    else:
        FitsPlot = aplpy.FITSFigure(img)
        #print " entrando no modo np array"
    smallest = np.amin(img)
    biggest = np.amax(img)


    
    if 'vmax' in opt.keys():
        if opt['vmax']=='Max':
           opt['vmax']=biggest
        if opt['vmin']=='Min':
           opt['vmin']=smallest
        FitsPlot.show_grayscale(vmin=opt['vmin'], vmax=opt['vmax'],invert=opt['invert'])
    else:
        if ('pmax' in opt.keys()) and ('pmin' in opt.keys()):
            FitsPlot.show_grayscale(pmin=opt['pmin'], pmax=opt['pmax'],invert=opt['invert'])
        else:
            FitsPlot.show_grayscale(invert=opt['invert'])
    #for i in range(0,100):
    #    print opt['invert']
    #    FitsPlot.show_grayscale(pmin=i*0.01, pmax=1,invert=False)
    #    FitsPlot.save("mediatrix_aplpy_withcuts/"+ID+"scaleMin"+str(i*0.01)+"Max"+str(1)+".png")
    #FitsPlot.show_grayscale(vmin=opt['vmin'], vmax=opt['vmax'],invert=opt['invert'])
    #print biggest
    #FitsPlot.save("mediatrix_aplpy_withcuts/"+IDtime+"scaleMin"+str(0)+"Max"+str(biggest)+".png")
    Length=0
    
      
    if keydots==True:
        #print pixels[0]
        pixels=np.array(pixels)
        #print pixels[0]
        #print type(pixels[0])
        #pixels[0]=pixels[0]+1
        #pixels[1]=pixels[1]+1
        E1,E2=get_extrema_2loops(pixels[0], pixels[1], 0 )
        Area=len(pixels[1])
        p1=pixels[0][E1]+ pixels[1][E1]*1j # the extreme points p_1 and p_2
        p2=pixels[0][E2]+ pixels[1][E2]*1j
        
        keydots=[p1,p2]
        keydots=find_keydots_c(p1,p2,pixels,image_ps,keydots,Area, method=opt['method'],alpha=opt['alpha'],near_distance=opt['near_distance'],max_level=opt['max_level'])
        
  
        if use_extremes==False:
            if len(keydots)>3:
                keydots.remove(p1)
                keydots.remove(p2)
            else:
                print "Impossible to reject extremes, too few points" 
        keyX=[]
        keyY=[]
        if opt['world_coord']=='True':
            for j in range(0,len(keydots)):
                
                keyX_aux,keyY_aux=FitsPlot.pixel2world(keydots[j].real,keydots[j].imag)
                keyX.append(keyX_aux)
                keyY.append(keyY_aux)
        else:
            for j in range(0,len(keydots)):
                
                keyX.append(keydots[j].real+1) #due to the image origin is at 1,1 not 0,0 like the matrix coordinates
                keyY.append(keydots[j].imag+1)
        
        FitsPlot.show_markers(keyY,keyX,c=opt['keys_color'],marker='.',zorder=500)
    if circle==True:
        if keydots==False:
            E1,E2=get_extrema_2loops(pixels[0], pixels[1], 0 )
        x_circ=[pixels[0][E1]+1,mediatrix_data['center'].real+1,pixels[0][E2]+1] #due to the image origin is at 1,1 not 0,0 like the matrix coordinates
        y_circ=[pixels[1][E1]+1,mediatrix_data['center'].imag+1,pixels[1][E2]+1] #due to the image origin is at 1,1 not 0,0 like the matrix coordinates
        FitsPlot.show_markers(y_circ,x_circ,c='g',marker='D',zorder=500)
        #print "as coordenadas sao y, x"
        #print mediatrix_data['center'].imag
        #print mediatrix_data['center'].real
        p1_vec=[pixels[0][E1],pixels[1][E1]] # the extreme points p_1 and p_2
        p2_vec=[pixels[0][E2],pixels[1][E2]]
        p3_vec=[mediatrix_data['center'].real+1,mediatrix_data['center'].imag+1] #due to the image origin is at 1,1 not 0,0 like the matrix coordinates
        x_c,y_c,r=three_points_to_circle(p1_vec,p3_vec,p2_vec)


        if r>0:
            FitsPlot.show_markers(np.array([p3_vec[1], y_c]), np.array([p3_vec[0], x_c]), layer=False, zorder=499,c='g',marker='.')
            #try: # circle works only in world coord
            #    xw,yw=FitsPlot.pixel2world(x_c,y_c)
                #print "this is the world coord"
                #print [xw,yw]
            #    xE1=pixels[0][E1]
            #    yE1=pixels[1][E1]
            #    xwE1,ywE1=FitsPlot.pixel2world(xE1,yE1)
            #    rw=sqrt((xw-xwE1)**2+(yw-ywE1)**2) #that relation is due to small angles
            #    FitsPlot.show_circles(yw, xw, rw, layer=False, zorder=499)
            #except:
            #    print "No WCS: Impossible to plot circle"
        else:
            print "impossible to define a circle "
    


        #A.scatter(keyY, keyX, s=20, c='b', marker='s')

    
    for i in range(0,len(mediatrix_data['origin'])):
        origin_x=mediatrix_data['origin'][i].real+1 #due to the image origin is at 1,1 not 0,0 like the matrix coordinates
        origin_y=mediatrix_data['origin'][i].imag+1 #due to the image origin is at 1,1 not 0,0 like the matrix coordinates
        end_x=mediatrix_data['end'][i].real+1 #due to the image origin is at 1,1 not 0,0 like the matrix coordinates
        end_y=mediatrix_data['end'][i].imag+1 #due to the image origin is at 1,1 not 0,0 like the matrix coordinates

        if opt['world_coord']==True:
            origin_x,origin_y=FitsPlot.pixel2world(origin_x,origin_y)
            end_x,end_y=FitsPlot.pixel2world(end_x,end_y)
        Length_aux=(origin_x - end_x)**2 + (origin_y - end_y)**2
        Length=Length+sqrt(Length_aux)
        d_x= end_x - origin_x
        d_y= end_y - origin_y
        #arr = Arrow(origin_y, origin_x, d_y, d_x, width=0.05*Length, fc=colors['vector'], ec='none',zorder=1000)
    #    print "vectors"
    #    print origin_x
    #    print origin_y
        FitsPlot.show_arrows(origin_y, origin_x, d_y, d_x,zorder=502 )
   
    #xmin, xmax = xlim()
    #ymin, ymax = ylim()
    #min_inc_axis=40
    #x_axis_length=(xmax+1*Length)-(xmin-1*Length)
    #y_axis_length=(ymax+1*Length)-(ymin-1*Length)
    #if  x_axis_length<min_inc_axis
    #A.axis("equal")
    #A.set_xlim(xmin-1*Length,xmax+1*Length)
    #A.set_ylim(ymin-1*Length,ymax+1*Length)    
    #ylabel("Y")
    #xlabel("X")
    #A.axis("equal")
    #title(out_title) 
    
    if save==True:
        FitsPlot.save(out_image)
        del FitsPlot
        del mediatrix_data
        del image_ps
        del img
        gc.collect
        return True
    else:
        return FitsPlot, mediatrix_data, image_ps.copy() 



