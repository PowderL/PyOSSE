""" this module is used to read in climatology data set 
    Authors: Dr. Liang Feng, University of Edinburgh
    History: v0.20, Nov.16, 2007
    Notes: 
      It is used to read the 'fake' average kernal 
""" 

from numpy import *
import gp_field as gpf
import time_module as tm
import gp_axis as gax
import geos_chem_def as gcd
def_ak_path=gcd.oco_ak_path
def_ak_prefix="aknorm"
def_axis_keys=['sza', 'od','Logp']  
def_axis_keys_no_vert=['sza', 'od']
def_axis_keys_no_vert.reverse()
def_axis_keys.reverse()
def_ak_name='oco-kernal'

def_ak_surf_type=['soil', 'ocean', 'snow', 'desert', 'conifer']
def_ak_pres=[  1., 70.0, 100.0, 200.0, 350.0, 450.0, \
               550.0, 650.0, 750.0, 850.0, 950.0, 1000.]
# def_ak_pres.reverse()
# print def_ak_pres

class average_kernal(gpf.gp_field):
    """ storage for average kernal data
    """
    
    def __init__(self, surface='snow', viewmode='nadir', \
                 gpname=None, do_debug=False):
        gpf.gp_field.__init__(self, surface, None)
        self.read_ak(surface, viewmode, do_debug)
        
    def read_ak(self, surface, viewmode, do_debug=False):
        """ read into the average kernal data
        Arguments:
        surface the surface type
        """
        if (viewmode=='nadir'):
            fname=def_ak_prefix+"_"+surface+".dat_1"
        else:
            fname=def_ak_prefix+"_"+surface+".dat_lg_glint"
        
        fpath=def_ak_path
        fc=open(fpath+fname, "r")
        sza=list()
        od=list()
        data=list()
        line1=fc.readline()
        lines=fc.readlines()
        fc.close()
        
        for line in  lines:
            terms=line.split()
            val=float(terms[0])
            if (not val in sza):
                sza.append(val)
            sval=terms[1].strip()
            sval=sval[2:]
            val=0.1*float(sval)
            if (val>1.0):
                val=0.1*val
            
            
            if (not val in od):
                od.append(val)
            for term in terms[2:]:
                data.append(float(term))
        # I need a order key list 
        
        sza=array(sza)
        nsza=size(sza)
        od=array(od)
        nod=size(od)
        pres=array(def_ak_pres)
        npres=size(pres)
        data=array(data)
        
        data=reshape(data, [nsza, nod, npres])
        if (do_debug):
            db_isels=[[0,1], [1,0]]
            for isel in db_isel:
                db_pos=isels[ii]
                # print 'sza, od', sza[db_pos[0]],od[db_pos[1]]
                # print 'ak:', data[db_pos[0],db_pos[1],:]
                
        
        ax_sza=gax.gp_axis('sza', sza)
        # print 'ax1'
        ax_od=gax.gp_axis('od', od)
        # print 'ax2'
        
        ax_logp=gax.gp_axis('Logp', log10(pres))
        axis_set=[ax_sza, ax_od, ax_logp]
        gpname=def_ak_name
        gpunit='NONE'
        self.add_gp(gpname, data)
        self.set_grid(axis_set)
        self.set_gp_attr(gpname, 'UNIT', gpunit)
        del data

    def get_ak_prof(self, rpos):
        """ this function is used to get a profile vs pressure or vertical grid
        Arguments:
        the name of the gp to be read
        """
        gpname=def_ak_name
        ref_x, profs=self.gp_get_prof(gpname, rpos, xaxis=2) # choose pressure

        return ref_x, profs


if __name__=='__main__':
    """ that is a test
    """
    from pylab import *
    from numpy import *
    ak=average_kernal('ocean')
    rpos=list()
    rpos=[10.0, 0.25]
    logp, prof=ak.get_ak_prof(rpos)
    # print shape(logp)
    #  print shape(prof)
    print 'interpolated ak at sza, od', rpos[:]
    print prof

    data=ak.get_gp(def_ak_name)
    p1=data[0,0,:]
    p2=data[0,1,:]
    p3=(p1+p2)/2.0
    print 'the exact ak at sza, od', rpos[:]
    print p3
    pres=10**(logp)
    # subplot(1,2,1)
    # semilogy(prof, pres)
    # semilogy(p1, pres, '--')
    # semilogy(p2, pres, '--')
    # ylim([1000.0, 1.0])
    
    # title('Intp ak vs nearby ak')
    ak2=average_kernal('conifer')
    logp2, prof2=ak2.get_ak_prof(rpos)

    ak3=average_kernal('desert')
    logp3, prof3=ak3.get_ak_prof(rpos)

    ak4=average_kernal('snow')
    logp4, prof4=ak4.get_ak_prof(rpos)


    ak5=average_kernal('soil')
    logp5, prof5=ak5.get_ak_prof(rpos)

    
    
    # subplot(2,2,1)
    semilogy(prof, pres)
    semilogy(prof2, pres)
    semilogy(prof3, pres)
    semilogy(prof4, pres)
    semilogy(prof5, pres)
    
    ylim([1000.0, 10.0])
    legend(['ocean', 'conifer', 'desert', 'snow', 'soild'])
    title('OCO average kernel')
    savefig('oco_ak.png')
    show()
    
    
    
    
    
    
    
    

   
    
    
    
    
