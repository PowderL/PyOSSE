import satellite_operator as sat_op
import surface_operator as sf_op
import gv_operator as gv_op
import geos_chem_def as gcdf

class obs_info:
    def __init__(self, obs_name, datapath, stype, mean_update, hm_update):
        self.obs_name=obs_name
        self.datapath=datapath
        self.stype=stype
        self.hm_update=hm_update
        self.mean_update=mean_update
        
        try:
            self.id=gcdf.instrument_id[self.obs_name]
        except KeyError:
            self.id=-999
        
def get_obs_list():

    obs_list=list()
     
    # aqua
    if (len(gcdf.view_type_list)==0):
        obs_name=gcdf.view_type+'_'+gcdf.view_mode
        datapath=gcdf.obs_path
        
        stype='satellite'
        aqua_obs=obs_info(obs_name, datapath, stype)
        obs_list.append(aqua_obs)
        
        
        if (gcdf.add_ersl):
            
            obs_name='ersl_flask'
            datapath='./esrl_obs/flask'
            stype='surface'
            ersl_obs=obs_info(obs_name, datapath, stype)
            
            obs_list.append(ersl_obs)
    else:
        nobs_type=len(gcdf.view_type_list)
        for iobs in range(nobs_type):
            view_type=gcdf.view_type_list[iobs]
            view_mode=gcdf.view_mode_list[iobs]
            # datapath=gcdf.run_path+'/'+view_type+'_obs/'
            datapath=gcdf.obs_path_list[iobs]
            
            stype=gcdf.view_station_list[iobs]
            if (view_mode.strip()<>""):
                obs_name=view_type+'_'+view_mode.strip()
            
            else:
                obs_name=view_type
            hm_update=gcdf.hm_update_list[iobs]
            mean_update=gcdf.mean_update_list[iobs]
            
            aqua_obs=obs_info(obs_name, datapath, stype,  mean_update, hm_update)
            obs_list.append(aqua_obs)
        
    return obs_list

def get_obs_op_list(obs_list):
    op_list=list()
    for oif in obs_list:
        if (oif.stype=='satellite'):
            op=sat_op.satellite_xgp(oif.obs_name, oif.datapath, \
                                    oif.id, oif.mean_update, oif.hm_update)
        elif (oif.stype=='surface'):
            op=sf_op.surface_gp(oif.obs_name, oif.datapath, oif.id, \
                                oif.mean_update, oif.hm_update)
        elif (oif.stype=='gv'):
            op=gv_op.gv_gp(oif.obs_name, oif.datapath, oif.id, \
                           oif.mean_update, oif.hm_update)
        else:
            print 'obs type error', oif
            ttt=raw_input()
        
        op_list.append(op)
    return op_list


    

    

    
        
        
