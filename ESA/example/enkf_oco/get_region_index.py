from numpy import *
def get_region_index(regname):
    

    t3_land_reg=['ice',\
                     'North American Boreal',\
                     'North American Temperate',\
                     'South American Tropical',\
                     'South American Temperate',\
                     'North Africa',\
                     'South Africa',\
                     'Eurasia Boreal', \
                     'Eurasia Temperate',\
                     'Tropical Asia', \
                     'Australia',  \
                     'Europe']
    
    t3_ocean_reg=['North Pacific Temperate', \
                      'West Pacific Tropics', \
                      'East Pacific Tropics', \
                      'South Pacific Temperate',\
                      'Northern Ocean',\
                      'Northern Atlantic Temperate',\
                      'Atlantic Tropics', \
                      'South Atlantic Temperate',\
                      'South Ocean', \
                      'Indian Tropical', \
                      'South Indian Temperate']
    
    
    # boreal land
    
    if (regname=='nam'):
        # NA
        sum_reg_idx=[0+1, 1+1]
        stitle='N. America)'

    elif (regname=='nat'):
        # NA
        sum_reg_idx=[1+1]
        stitle='N. America temperate'
    elif (regname=='nab'):
        # NA
        sum_reg_idx=[0+1]
        stitle='N. America boreal'
        
    elif (regname=='sat'):
        # S
        sum_reg_idx=[1+2]
        stitle='S. America tropical'
    elif (regname=='sam'):
        # S
        sum_reg_idx=[1+3]
        stitle='S. America '
        
    elif (regname=='naf'):
        sum_reg_idx=[1+4]
        stitle='N. Africa'
        
    elif (regname=='saf'):
        # TA 
        sum_reg_idx=[1+5]
        stitle='S. Africa'


    elif (regname=='eab'):
        # TA
        
        sum_reg_idx=[1+6]
        stitle='Eurasia boreal'

    
    elif (regname=='eat'):
        # TA 
        sum_reg_idx=[1+7]
        stitle='Eurasia temperate'

    elif (regname=='ta'):
        # TA 
        sum_reg_idx=[1+8]

        stitle='Tropical Asia'
        
    elif (regname=='aus'):
        sum_reg_idx=[1+9]
        stitle='Australia'

    elif (regname=='eu'):
        sum_reg_idx=[1+10]
        stitle='Europe'
        
    elif (regname=='land'):
        sum_reg_idx=range(0, 12)
        stitle='Land'
    elif (regname=='ocean'):
        sum_reg_idx=range(12, 23)
        stitle='Ocean'
    elif (regname=='global'):
        sum_reg_idx=arange(0, 23)
        stitle='Global'
        
    elif (regname=='north_land'):
        sum_reg_idx=[0+1,1+1, 6+1, 7+1, 10+1]
        stitle='North land'
    
    elif (regname=='tropical_land'):
        sum_reg_idx=[1+2,1+4, 1+8]
        stitle='Tropical land'
        
    elif (regname=='south_land'):
        sum_reg_idx=reg_ids=[1+3, 1+5, 1+9]
        stitle='South land'
    
    
    elif (regname=='npto'):
        sum_reg_idx=reg_ids=[1+11]
        stitle='North Pacific Temperate'

    
    elif (regname=='wpto'):
        sum_reg_idx=reg_ids=[1+12]
        stitle='West Pacific Tropics'

        
    elif (regname=='epto'):
        sum_reg_idx=reg_ids=[1+13]
        stitle='east Pacific Tropics'

    
    elif (regname=='spto'):
        sum_reg_idx=reg_ids=[1+14]
        stitle='South Pacific Temperate'

    
    elif (regname=='no'):
        sum_reg_idx=reg_ids=[1+15]
        stitle='Northern Ocean'
    
    elif (regname=='nato'):
        sum_reg_idx=reg_ids=[1+16]
        stitle='Northern Atlantic Temperate'
        

    elif (regname=='atto'):
        sum_reg_idx=reg_ids=[1+17]
        stitle='Atlantic tropical'

    

    elif (regname=='sato'):
        sum_reg_idx=reg_ids=[1+18]
        stitle='South Atlantic Temperate'


    
    elif (regname=='so'):
        sum_reg_idx=reg_ids=[1+19]
        stitle='South Ocean'


    
    
    elif (regname=='idto'):
        sum_reg_idx=reg_ids=[1+20]
        stitle='Indian Tropical'


    elif (regname=='sidto'):
        sum_reg_idx=reg_ids=[1+20]
        stitle='South Indian Temperate'
    

    elif (regname=='tropical_oceans'):
        sum_reg_idx=reg_ids=[1+12, 1+13, 1+17, 1+20]
        stitle='Tropical Oceans'
        
    
     
    elif (regname=='north_oceans'):
        sum_reg_idx=reg_ids=[1+11, 1+14, 1+15, 1+16]
        stitle='North Oceans'

    elif (regname=='south_oceans'):
        sum_reg_idx=reg_ids=[1+18, 1+19, 1+21]
        stitle='South Oceans'

        


    

    
        
    return sum_reg_idx, stitle



             
              



