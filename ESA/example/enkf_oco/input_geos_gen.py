import numpy as npy 
import ESA.util.time_module as tm

def create_new_input_file(run_step,\
                              run_path,\
                              data_path,\
                              yyyy_lst,\
                              doy_lst,\
                              pbuse_lst,\
                              pb_start=1,\
                              member_start=1, \
                              member_end=12, \
                              pbflnm='CO2_EMISSION_EN',\
                              tmpfile="input.geos.temp", \
                              newfile="input.geos.new" , \
                              newext=None,\
                              enr_step_lst=None, \
                              enr_st_doy=None,\
                              enr_st_yyyy=None,\
                              enr_pbuse_lst=None,\
                              do_bk_run=4,\
                              **keywords \
                              ):
    
    """ 
    create the new input to drive an  ensemble or single tracer run
    
    Inputs:
    ===================================================
    1. run_step:<int>: step of ensemble run 
    2. run_path:<str>: directory for job run 
    3. data_path:<str>:
    
    4. yyyy_lst:<list, t:int> year for starting time at each step  
    5. doy_lst:<list, t:int>: day of year for starting doy at each step 
    
    5. pbuse_lst:<list, t:int>: '1' means flux perturbation file  will be used 
    6. pb_start:<int>: the starting layer (or region) of the flux perturbation file 
    7. member_start:<int>: the starting member of the ensemble run 
    8. member_end:<int>:  the ending member of the ensemble run 
    9. pbflnm:<str>: name of the pb file 
    10. tmpfile:<str>: name of input.geos template 
    11. newfile:<str>: name of the new input.geos 
    12. newext:<str>:  extension for output file names
    13. enr_step_lst:<list>: ensemble run steps 
    14. enr_st_doy:<list>: ensemble run starting doy.   
        If it is not given, doy_lst will be used 
    15. enr_st_yyyy: <list>: ensemble run starting doy.  
        If it is not  given, yyyy_lst will be used
    
    16. enr_pbuse_lst: <list>: ensemble run flux perturbation switch 
    17. do_bk_run:<int>: run type
       5--run forced by posterior flux 
       4--ensemble run 
       3--single tracer run with prior flux 
    
   
    18. keywords:<dict>: dictionary for extra inputs 
    
    reserved keywords:
    -->1. time_start:<str>: the starting time (yyyymmdd)
    -->2. time_end:<str>: the ending time (yyyymmdd)
    
    """
    
    # S1: find start and end time 
    
    ntime=npy.size(yyyy_lst)
    if ('time_start' in keywords): 
        ts=keywords['time_start']
        
        syyyy=ts[0:4]
        yst=int(syyyy)
        
        smm=ts[4:6]
        mmst=int(smm)
        
        sdd=ts[6:8]
        ddst=int(sdd)
        
        dst=tm.day_of_year(yst, mmst, ddst)
        
        
    else:
        
        yst=yyyy_lst[0]
        dst=doy_lst[0]

    tst=tm.doy_to_utc(dst, sec=0, yyyy=yst)
    tst=tst.replace('-', '')
    tst=tst.replace(':', '')
    print 'input_geos_gen->starting time', tst

    if (enr_step_lst==None):
        # #c: make a dummpy step list
        
        iday0=doy_lst[0]
        for iday in doy_lst[1:]:
            if (iday>iday0):
                enr_step_lst=iday-iday0
                break
            else:
                iday0=iday
        
    if ('time_end' in keywords): 
        
        te=keywords['time_end']
        syyyy=te[0:4]
        yed=int(syyyy)
        
        smm=te[4:6]
        mmed=int(smm)
        
        sdd=te[6:8]
        dded=int(sdd)
        
        ded=tm.day_of_year(yed, mmed, dded)
    else:
        
        yed=yyyy_lst[1] #-1
        ded=doy_lst[1]  #-1
    
    
    tend=tm.doy_to_utc(ded, sec=0, yyyy=yed)
    tend=tend.replace('-', '')
    tend=tend.replace(':', '')
    
    print 'input_geos_gen->end time', tend
    
    # S3: read in input.geos.template file 
    
    fin=open(tmpfile, "r")
    lines=fin.readlines()
    print  len(lines)
    fin.close()

    # S4: open new file 

    fout=open(newfile, "w")
    section_start=0
    colwidth=25
    line_count=0

    ntracer=member_end-member_start+1

    # S5: set extension for file outputs 
    if (newext==None):
        enafix=r'ST%3.3d.EN%4.4d-EN%4.4d' % (run_step, member_start, member_end)
    else:
        enafix=newext
    
    
    # S6: loop over input lines, and replace default values with current settings 

    for line in lines:
        # #c: check for section heads
        
        if ("SIMULATION MENU" in line):  # 1. simulation  
            section_start=1
            line_count=0
            fout.write(line)
        elif ("TRACER MENU" in line):   # 2. tracer 
            section_start=2
            line_count=0
            fout.write(line)
        elif ("ND51 MENU" in line):    # 3.  ND51
            section_start=3
            line_count=0
            fout.write(line)
        elif("OUTPUT MENU" in line):   #10. output  
            section_start=10
            line_count=0
            fout.write(line)
        elif ("ENSEMBLE MENU" in line):    #4 ensemble  
            section_start=4
            line_count=0
            fout.write(line)
        elif ("DIAGNOSTIC MENU" in line):  # 5 diagnostic
             section_start=5
             line_count=0
             fout.write(line)
        elif("GAMAP MENU" in line):        #6 GAMAP menu 
            section_start=6
            line_count=0
            fout.write(line)
        
        elif (section_start>0):

            # #c: if the line is in one of the sections. 
            # #c: replace the place-holder with inputs
            

            line_count=line_count+1
            line_head=line[0:colwidth]
            line_left=line[colwidth:]

            if (section_start==1):    #  simulation  
                if (line_count==1):
                    new_line=line_head+" "+tst+"   "
                    fout.write(new_line+"\n")
                elif (line_count==2):
                    line_head=line[0:colwidth]
                    new_line=line_head+" "+tend+"   "
                    fout.write(new_line+"\n")
                elif (line_count==3):
                    print line_left
                    line_left=line_left.replace('$RUNPATH', run_path)
                    new_line=line_head+line_left
                    fout.write(new_line)
                elif (line_count==4):
                    print line_left
                    line_left=line_left.replace('STYYY.ENXXXX-ENXXXX', enafix)
                    new_line=line_head+line_left
                    fout.write(new_line)
                elif (line_count==6):  
                    line_left=line_left.replace('$DATAPATH', data_path)
                    line_left=line_left.replace('STYYY.ENXXXX-ENXXXX', enafix)
                    new_line=line_head+line_left
                    fout.write(new_line)
                elif (line_count==13):
                    line_left=line_left.replace('$RUNPATH', run_path)
                    new_line=line_head+line_left
                    fout.write(new_line)
                elif (line[0:2]=='--'):
                    line_count=0
                    section_start=0
                    fout.write(line)
                else:
                    fout.write(line)
            elif (section_start==2): # tracer menu 
                if (line_count==1):
                    fout.write(line)
                elif(line_count==2):
                    line_head=line[0:colwidth]
                    snum=r'%4d' % (ntracer)
                    new_line=line_head+" "+snum+"   "+"\n"
                    fout.write(new_line)
                elif (line_count==3):
                    fout.write(line)
                elif (line_count==4):
                    line_head=line[0:colwidth]
                    tracer_no=line[colwidth:colwidth+5]
                    line_left=line[colwidth+5:]
                    tracer_id=1
                    tracers=list()
                    
                    for itracer in range(member_end-member_start+1):
                        st1=r'Tracer #%3.3d ENSEM' % (itracer+member_start)
                        lst1=len(st1)
                        new_head=st1[:]+' '*(colwidth-1-lst1)+":"
                        st1=r'%5d' % (tracer_id)
                        lst1=len(st1)
                        new_tracer_no=(5-lst1)*' '+st1
                        new_line=new_head+new_tracer_no+line_left
                        fout.write(new_line)
                        tracers.append(tracer_id)
                        tracer_id=tracer_id+1
                        
                elif (line[0:2]=='--'):
                    line_count=0
                    section_start=0
                    fout.write(line)
                else:
                    pass
            elif (section_start==3): # ND51 menu 
                if (line_count==1):
                    fout.write(line)
                elif (line_count==2):
                    line_left=line_left.replace('$DATAPATH', data_path)
                    line_left=line_left.replace('STYYY.ENXXXX-ENXXXX', enafix)
                    new_line=line_head+line_left
                    fout.write(new_line)
                elif (line_count==3):
                    line_left=""
                    for itracer in [1, member_end-member_start+1]:
                        stracer=r' %d' % (itracer)
                        line_left=line_left+stracer
                                            
                    line_left=line_left+' '+'196 198 199 200 201'  # output other information
                    new_line=line_head+line_left+' \n'
                    print new_line
                    # tx=raw_input()
                    
                    fout.write(new_line)
                    
                    
                elif (line[0:2]=='--'):          
                    line_count=0
                    section_start=0
                    fout.write(line)
                else:
                    fout.write(line)
            elif (section_start==4): #ensemble  menu
                if (line_count==1):
                    st1=r' %d  ' % member_start
                    new_line=line_head+st1+' \n'
                    fout.write(new_line)
                elif (line_count==2):
                    st1=r' %d  ' % member_end
                    new_line=line_head+st1+' \n'
                    fout.write(new_line)
                elif (line_count==3):
                    st1=r' %d  ' % do_bk_run
                    new_line=line_head+st1+' \n'
                    fout.write(new_line)

                elif (line_count==4):
                    st1=r' %d ' % ntime
                    new_line=line_head+st1+' \n'
                    fout.write(new_line)
                
                elif (line_count==5):
                    st1=""
                    if (enr_st_doy==None):
                        enr_st_doy=doy_lst
                    for doy in enr_st_doy:
                        sdoy=r'%d' % doy
                        st1=st1+' '+sdoy
                    st1=st1+' '
                    new_line=line_head+st1+' \n'
                    fout.write(new_line)
           
                elif (line_count==6):
                    st1=""
                    if (enr_st_yyyy==None):

                        for yyyy in yyyy_lst:
                            syyyy=r'%d' % yyyy
                            st1=st1+' '+syyyy
                        
                    else:

                        for yyyy in enr_st_yyyy:
                            syyyy=r'%d' % yyyy
                            st1=st1+' '+syyyy
                    
                    st1=st1+' '
                    new_line=line_head+st1+' \n'
                    fout.write(new_line)
                
                elif (line_count==7):   
                    st1=r' %d  ' % pb_start
                    new_line=line_head+st1+' \n'
                    fout.write(new_line)
                
                elif (line_count==8):
                    st1=""
                    if (enr_pbuse_lst==None):
                        for pbuse in pbuse_lst:
                            spbuse=r'%d' % pbuse
                            st1=st1+' '+spbuse
                        
                    else:
                        
                        for pbuse in enr_pbuse_lst:
                            spbuse=r'%d' % pbuse
                            st1=st1+' '+spbuse
                    
                    st1=st1+' '
                    new_line=line_head+st1+' \n'
                    fout.write(new_line)
                
                elif (line_count==9): # 'flnm'
                
                    new_line=line_head+' '+pbflnm.strip()+'\n'
                    fout.write(new_line)
                elif (line[0:2]=='--'):          
                    line_count=0
                    section_start=0
                    fout.write(line)
            
            elif(section_start==5): #diagnostic  menu
                
                if (line_count==1):
                    line_left=line_left.replace('$DATAPATH', data_path)
                    line_left=line_left.replace('STYYY.ENXXXX-ENXXXX', enafix)
                    new_line=line_head+line_left
                    fout.write(new_line)
                
                elif (line[0:2]=='--'):          
                    line_count=0
                    section_start=0
                    fout.write(line)
                else:
                    fout.write(line)
            
            elif(section_start==6): #gamap  menu
                
                if (line_count==1):
                    line_left=line_left.replace('$DATAPATH', data_path)
                    line_left=line_left.replace('STYYY.ENXXXX-ENXXXX', enafix)
                    new_line=line_head+line_left
                    fout.write(new_line)
                elif (line_count==2):
                    line_left=line_left.replace('$DATAPATH', data_path)
                    line_left=line_left.replace('STYYY.ENXXXX-ENXXXX', enafix)
                    new_line=line_head+line_left
                    fout.write(new_line)
                elif (line[0:2]=='--'):          
                    line_count=0
                    section_start=0
                    fout.write(line)
                else:
                    fout.write(line)

            elif (section_start==10): #output menu
                
                if (line_count==1):
                    id_cnt=0
                    yyyy_list=list()
                    mm_list=list()
                    dd_list=list()
                    
                    for idoy in doy_lst:
                        cyyyy, cmm, cdd=tm.doy_to_time_array(idoy, yyyy_lst[id_cnt])
                        yyyy_list.append(cyyyy)
                        mm_list.append(cmm)
                        dd_list.append(cdd)
                        id_cnt=id_cnt+1
                    
                    yyyy_list=npy.array(yyyy_list)
                    mm_list=npy.array(mm_list)
                    dd_list=npy.array(dd_list)
                    
                # #c: force diagnostic outputs at end of the rerun 
                
                if ((line_count>=1) and (line_count<13)):
                    line_left=line_left.replace('3', '0')
                    sel_mm=npy.where(mm_list==line_count)
                    
                    if (npy.size(sel_mm)>0):
                        if (npy.size(sel_mm)==1):
                            sel_mm=[npy.squeeze(sel_mm)]
                        else:
                            sel_mm=npy.squeeze(sel_mm)
                            
                        for isel in sel_mm:
                                            
                            sel_dd=dd_list[isel]
                            ispace=0
                            # skip the line
                            
                            for ichar in line_left:
                                if (ichar<>' '):
                                    break
                                else:
                                    ispace=ispace+1
                            
                            sel_dd=sel_dd+ispace
                            
                            print line_left
                            
                            if (sel_dd>=len(line_left)):
                                sel_dd=len(line_left)-1
                        
                            if (sel_dd==1):
                                tline_1='3'
                            
                            else:
                                tline_1=line_left[0:sel_dd-1]
                                tline_1=tline_1+'3'
                            
                                tline_1=tline_1+line_left[sel_dd:]
                            line_left=tline_1
                            
                    new_line=line_head+line_left
                    fout.write(new_line)        
                elif (line[0:2]=='--'):  # the section end line          
                    line_count=0
                    section_start=0
                    fout.write(line)

            else:
                fout.write(line)
                    
            
            

        else:
            fout.write(line)
    
    fout.close()
    
    print 'reach the end'
    #  tx2=raw_input()


# <<< TEST >>> 

if (__name__=="__main__"):
    run_path='./'
    data_path='./enkf_output/'
    yyyy_lst=[2009, 2009, 2009]
    doy_lst=[1, 32, 34]
    pbuse_lst=[1,1,1]
    pb_start=1
    
    
    create_new_input_file(1,\
                              './',\
                              './enkf_output/',\
                              yyyy_lst,\
                              doy_lst,\
                              pbuse_lst,\
                              pb_start=pb_start,\
                              member_start=1, \
                              member_end=12, \
                              pbflnm='CO2_EMISSION_EN',\
                              tmpfile="input.geos.temp", \
                              newfile="input.geos.new" , \
                              enr_step_lst=None, \
                              enr_st_doy=None,\
                              enr_st_yyyy=None,\
                              enr_pbuse_lst=None,\
                              do_bk_run=4)
                              

    

                  
                
