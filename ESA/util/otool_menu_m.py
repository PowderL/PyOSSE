"""
class and functions for extracting options from  xml or text configuration  files

  Authors: L. Feng, Edinburgh University
  History: v0.9, 2012.12.08
  History: v0.95, 2013.03.24
  

  
  Class
  =================================================================
  1. menu_cl: class for menu list (ie., configure options)
  
  Functions:
  ==================================================================
  1. txt_to_menu: read text file to class menu_cl 
  2. xml_to_menu: read xml file to  class menu_cl

  3. update_menu_from_text: update menu class using text lines 
  4. update_menu_from_element: update menu class using element 
  
"""

import  numpy as npy
import ESA.util.time_module as  tm
import ESA.util.geo_constant as gc
import ESA.util.otool_obj as oob
import ESA.util.message_m as msm

#====<CLASS>====
 

class menu_cl:
   
    """
    Members:
    --------------------------------------------------------
    1. menu_item:<dict>: menu chain 
    2. attr_dict:<dict>: attributes
    

    Functions:
    --------------------------------------------------------
    1. __getitem__: overriding index function ('[]') to check dictionary entry of self.menu_item
    2. delete_menu:  remove item from self.menu_item
    3. update_menu:  Add or update menu item
    4. get_attr:  get attribute 
    5. set_attr:  set attribute 
    6. print_menu:  print out what is in the menu chain
    
    
    Notes:
    -------------------------------------------------------------
    1. menu_cl is designed to be a container for dictionary trees so it can have a nested structure like:
       Menu P---->Menu C1 
       ---->Menu C2
       etc
       Consquently item under sub-menu (for example menu C1) will be accessed using the 
       full name of menu chain separated by '.' (for example  P.C1.itemname)  
    
    2. special words:   __load; __dict; __list;
       ---> __load: If an item value is given as a string in form of  '__load:module_name:module_member' ,  
       self.update_menu will try to import this module, and then assign its corresponding member to the menu item. 
    
       ---> __dict: if an item value is given as a string in form of '__dict:a=x, b=y, ..., 
       self.update_menu will create a dictionary as {'a':x, 'b':y,,.} under the item's name
       
       ---> __list: if an item value is given as a string in form of '__list:a, b, ..., 
       self.update_menu will create a list as [a, b, c] under the item's name
         
    """
    
    
    def __init__(self, **menu_dict):
        """initialization
        
        Inputs:
        ----------------------------------
        1. menu_dict:<dict>: menu 

        
        """
        self.menu_item={}
        for keynm in menu_dict:
            keyval=menu_dict[keynm]
            self.update_menu(keynm, keyval)
        self.attr_dict={}
        self.attr_dict.update({'ot_type':oob.ot_menu})
        
    def get_attr(self, attr_name):
        """
        get attribute of the object 
        
        Inputs:
        ------------------------------------------
        1.attr_name:<str>: attribut name 

        Returns:
        ------------------------------------------
        1. val:<obj>: attribute
        
        """
        attr_name=attr_name.strip()
        
        if (attr_name in self.attr_dict):
            return self.attr_dict[attr_name]
        else:
            msm='No attribute found for '+attr_name
            
            return None 
    
    

    def set_attr(self, attr_name, attr_value):
        
        """
        set attribute  of the object 
        
        Inputs:
        ------------------------------------------
        1.attr_name:<str>: attribute name 
        2.attr_value:<obj>: attribute value 
        
        """
        attr_name=attr_name.strip()
        self.attr_dict.update({attr_name:attr_value})
        
        
    
    def __getitem__(self, itemname):
        """
        Overriding index function

        Inputs:
        ----------------------------------------------------
        1. itemname: <str>: name of the menu item to be checked 


        Outputs:
        --------------------------------------------------
        1. value: the menu value if found 


        Notes:
        --------------------------------------
        1. if itemname contains '.' (for example  M1.M2.M3), the value will be searched down through chain 
        of menu classes (such as M1-->M2--M3)
        
        
        
        """
        if ('.' in itemname):
            # # we need to look through the menu chain 
            
            menu_chain_lst=itemname.split('.')
            nsub=len(menu_chain_lst)
            # search down through top menu 
            
            value=self.menu_item
            for isub in range(nsub):
                mname=menu_chain_lst[isub]
                value=value[mname]
                if (value==None):
                    return value
                
            return value
        
        else:
            
            if (itemname in self.menu_item):
                value=self.menu_item[itemname]
                return value
            else:
                msg='No entry found: '+itemname
                msm.show_err_msg(msg)
                return None
    
    def update_menu(self, menu_name, menu_val):
        """
        Add or update menu item

        Inputs:
        --------------------------------------------
        1. menu_name:<str>: name of the menu 
        2. menu_val: <obj>: value 
        
        
        Notes:
        --------------------------------------
        1. if itemname contains '.' (for example  M1.M2.M3), the value will be added or updated down through 
        the chain of menu classes (such as M1-->M2--M3)
        
        2. special words:   __load; __dict; __list;
       ---> __load: If an item value is given as a string in form of  '__load:module_name:module_member' ,  
       self.update_menu will try to import this module, and then assign its corresponding member to the menu item. 
    
       ---> __dict: if an item value is given as a string in form of '__dict:a=x, b=y, ..., 
       self.update_menu will create a dictionary as {'a':x, 'b':y,,.} under the item's name
       
       ---> __list: if an item value is given as a string in form of '__list:a, b, ..., 
       self.update_menu will create a list as [a, b, c] under the item's name
       
       ---> __range: if an item value is given as a string in form of '__range:a, b, c, 
       self.update_menu will create a array as arange(a, b, c) under the item's name
       
       
       """
        
        if ('.' in menu_name):
            menu_chain_lst=menu_name.split('.')
            nsub=len(menu_chain_lst)
            cmenu=self.menu_item
            for isub in range(nsub-1):
                mname=menu_chain_lst[isub]
                cmenu=cmenu[mname]
                if (cmenu==None):
                    return cmenu
            
            last_name=menu_chain_lst[nsub-1]
            
            cmenu.update_menu(last_name, menu_val)
                        
        else:
            # #T: assign the value to menu item 
            if (oob.get_ot_type(menu_val)==oob.ot_str):
                if (menu_val=='None'):
                    # None is given
                    
                    self.menu_item.update({menu_name:None})
                
                elif (menu_val=='True'):
                    self.menu_item.update({menu_name:True})
                
                elif (menu_val=='False'):
                    self.menu_item.update({menu_name:False})
                
                elif ('__load:' in menu_val):
                    # #c: if value is from other module 
                    opname, modname, varname=menu_val.split(':')
                    # print opname, modname, varname
                    # print menu_name
                    
                    modname=modname.strip()
                    varname=varname.strip()
                    tmp_md=__import__(modname, globals(), locals(),-1)
                    val=getattr(tmp_md, varname)
                    self.menu_item.update({menu_name:val})
                
                elif ('__dict:' in menu_val):
                    # #c: if value is one dictionary     
                    
                    opname, svalue=menu_val.split(':')
                    svalue=svalue.strip()
                    
                    if (svalue==''):
                        dict_val={}
                    elif (svalue=='Nil'):
                        dict_val={}
                    else:
                        terms=svalue.split(',')
                        dict_val={}
                        for sterm in terms:
                            skn, skval=sterm.split('=')
                            skn=skn.strip()
                            skval=skval.strip()
                            
                            try:
                                tx=eval(skval)
                            except:
                                tx=skval
                    
                            dict_val.update({skn:tx})
                        
                    self.menu_item.update({menu_name:dict_val})
                    
                elif ('__list:' in menu_val):
                    # #c: if value is one dictionary     
                    
                    opname, svalue=menu_val.split(':')
                    svalue=svalue.strip()
                    
                    if (svalue==''):
                        lst_val=[]
                    elif (svalue=='Nil'):
                        lst_val=[]
                    else:
                        terms=svalue.split(',')
                        lst_val=[]
                        for skval in terms:
                            skval=skval.strip()
                            try:
                                tx=eval(skval)
                            except:
                                tx=skval
                            
                        
                            lst_val.append(tx)
                    
                    self.menu_item.update({menu_name:lst_val})

               
                elif ('__range:' in menu_val):
                    # #c: if value is given a array create by menu 
                    
                    
                    opname, svalue=menu_val.split(':')
                    svalue=svalue.strip()
                    
                    if (svalue==''):
                        lst_val=[]
                    elif (svalue=='Nil'):
                        lst_val=[]
                    else:
                        terms=svalue.split(',')
                        if (len(terms)==2):
                            val_st=eval(terms[0])
                            val_end=eval(terms[1])
                            lst_val=npy.arange(val_st, val_end)
                        elif (len(terms)==3):

                            val_st=eval(terms[0])
                            val_end=eval(terms[1])
                            dv=eval(terms[2])

                            
                            lst_val=npy.arange(val_st, val_end, dv)
                        
                    
                    self.menu_item.update({menu_name:lst_val})

                
               
                else:
                    # #c: simple obj
                    try:
                        tx=eval(menu_val)
                    except:
                        tx=menu_val
                    
                    self.menu_item.update({menu_name:tx})
            
            else:
                # #c: simple obj
                
                self.menu_item.update({menu_name:menu_val})

    def delete_menu(self, menu_name):
        
        """
        
        remove item from self.menu_item
        Inputs:
        ------------------------------------
        1. menu_name:<str>: name of the menu to be deleted
        
        """

        
        if ('.' in menu_name):
            
            menu_chain_lst=menu_name.split('.')
            
            nsub=len(menu_chain_lst)
            cmenu=self.menu_item
            for isub in range(nsub-1):
                mname=menu_chain_lst[isub]
                cmenu=cmenu[mname]
                if (cmenu==None):
                    return cmenu
            
            
            last_name=menu_chain_lst[nsub-1]
            cmenu.delete_menu(last_name)
            
        else:
            del self.menu_item[menu_name]
    
        
    def print_menu(self, parentname=''):
        
        """

        print out what is in the menu chain
        
        """
        for cname in self.menu_item:
            
            
            cmenu=self.menu_item[cname]
            
            
            if (type(cmenu)==type(menu_cl)):
                
                print parentname+cname+':', cmenu 
                
            else:
                if (oob.get_ot_type(cmenu)==oob.ot_menu):
                    cparentname=parentname+cname+'.'
                    cmenu.print_menu(cparentname)
                else:
                    print parentname+cname+':', cmenu 
                    
#<<< FUNCTIONS >>> 

def update_menu_from_element(parent_menu, parent_mname, element):
    """
    Updatemenu class using element 
    
    Inputs:
    -------------------------------------------
    1. parent_menu:<menu_cl>: menu class 
    2. parent_mname:<str>: name of the menu class
    3. element:<elment>: element class
    
    Returns:
    ------------------------------------
    1. parent_menu:<menu_cl>: updated menu class 

    
    """
    
    if (element.tag=='MENU'):
        # #c: element is a (sub) menu block 
        
        if ('name' in element.attrib):
            mname=element.attrib['name']
        else:
            mname='CM'
        
        if (len(parent_mname)>0):
            full_mname=parent_mname+'.'+mname
        else:
            full_mname=mname
            
        
        
        # ## add the attributes to menu 
        cmenu=menu_cl(**element.attrib)
        
        parent_menu.update_menu(full_mname, cmenu)
        # add other item to attributes
        for sc in element: 
            if (sc.tag=='MENU'):
                # call back the function 
                
                update_menu_from_element(parent_menu, full_mname, sc)
            else:
                full_itemname=full_mname+'.'+sc.tag # item name with parents. 
                parent_menu.update_menu(full_itemname, sc.text)
    
    else:
        
        full_itemname=element.tag
        parent_menu.update_menu(full_itemname, element.text)
        
        
        # other things

        
        
        
    return parent_menu



def xml_to_menu(flnm):
    """
    function to convert xml file into menu 

    Inputs:
    -------------------------------------
    1.flnm :<str>: name for xml
    
    Returns:
    ----------------------------------
    1. root_menu:<menu_cl>: menu read from file

    """
    
    import xml.etree.ElementTree as ET
    tree = ET.parse(flnm)
    root = tree.getroot()

    root_menu=menu_cl()
    root_mname=""
    for child in root:
        root_menu=update_menu_from_element(root_menu, root_mname, child)
    
    

    return root_menu

def update_menu_from_text(parent_menu, lines):
    
    """
    Update menu class using text lines 
    
    Inputs:
    -------------------------------------------
    1. parent_menu:<menu_cl>: parent menu 
    2. lines: lines will be go through
    
    Returns:
    -----------------------------------------------------
    1. parent_menu:<menu_cl>: update parent menu 
    """
    nline=len(lines) # line number 
    icur=0 # index for current line 

    # S1: loop over each line 
    
    for iline in range(nline):
        if (icur==nline):
            # #c:when reach the end
            break
        
        line=lines[icur]
        line=line.strip()
        icur=icur+1
        
        if (line==''):
            # empty line 
            pass
        elif (line[0]=='#'):
            # #T: comment line or control (sub menu) line 
            
            if (len(line)<5):
                # #c: short comment line 
                pass
            
            elif (line[0:5]=='#MENU'):
                # #c: encounter menu block 
                
                icnt=1
                blk_line=list()
                print 'processing menu begin:', line
                print 'line +1:', lines[icur+1]
                print 'line +2:', lines[icur+2]
                print  '---------------------------------------------'   
                print '  '
                
                
                # fetch lines within the block 
                for  cline in lines[icur:]:
                    cline=cline.strip()
                    icur=icur+1
                    
                    if (cline==''):
                        pass
                    elif (cline[0]=='#'):
                        if (len(cline)<5):
                            pass

                        elif (cline[0:5]=='#MENU'):
                            icnt=icnt+1
                            
                        elif (cline[0:5]=='#MEND'): # block end
                            icnt=icnt-1
                            if (icnt==0):
                                break # block search done
                        else:
                            pass
                    
                    blk_line.append(cline)
                
                
                # creat a sub menu 
                cmenu=menu_cl()
                update_menu_from_text(cmenu, blk_line)
                mname=cmenu['name']
                
                if (mname==None):
                    mname='SM'
                # add cmenu to cmenu
                
                parent_menu.update_menu(mname, cmenu)
            else:
                # ##c: comment line 
                
                pass
            # #c: endif (len(line)<5)

        else:
            # #T: menu item line 
            
            mname, mval=line.split('|')
            parent_menu.update_menu(mname, mval)
        
        # ##c: end if (line=='')
            
    
    return parent_menu

def preprocess_text(lines):
    """
    preprocess text line 
    Inputs:
    ----------------------------------------------------
    1. lines:<list, t:str>: lines to be prepressing 
    
    Returns:
    --------------------------------------------
    1. new_lines<list, t:str>: lines after being  prepressing
    2. para_dict:<dict>: parameters found. 
    

    Notes:
    1. the parameter block begin with #SETP and end with #ENDP
    2. defined parameters will be used to replace $PARANAME$ in text lines. 


    """
    nline=len(lines)
    # loop over all line 

    txt_line=list()
    para_dict={}
    # search parameter 
    
    is_in_blk=0
    icur=0
    for iline in range(nline):
        line=lines[icur]
        line=line.strip()
        
        if (line==''):
            pass
        elif (line[0]=='#'):
            
            if (len(line)<5):
                txt_line.append(line)
            elif (line[0:5]=='#SETP'):
                is_in_blk=is_in_blk+1
            elif (line[0:5]=='#ENDP'):
                is_in_blk=is_in_blk-1
            else:
                txt_line.append(line)
        else:
            if (is_in_blk>0):
                vname, val=line.split('=')
                vname=vname.strip()
                val=val.strip()
                para_dict.update({'$'+vname+'$':val})
            else:
                txt_line.append(line)
        
        
        icur=icur+1
    
    # end for
    # replace parameter 
    
    new_lines=list()
    for line in txt_line:
        for varname in para_dict:
            var=para_dict[varname]
            line=line.replace(varname, var)
        new_lines.append(line)
    
        
    return new_lines, para_dict

    

def txt_to_menu(flnm):

    """
    function to convert txt file into menu 

    Inputs:
    
    -------------------------------------
    1.flnm :<str>: name for txt file 
    
    Returns:
    ----------------------------------
    1. root_menu:<menu_cl>: menu read from file


    
    
    """
    
    fl=open(flnm, 'r')
    lines=fl.readlines()
    fl.close()
    lines, para_dict=preprocess_text(lines)
    root_menu=menu_cl()
    
    root_menu=update_menu_from_text(root_menu, lines)
    
    return root_menu


    
#<<< TEST >>> 
    

if (__name__=='__main__'):
    z=menu_cl(x=30, y=40)
    cl_menu=menu_cl(a=10, b=20, z=z)
    
    print cl_menu['a']
    print cl_menu['b']
    print cl_menu['z.x']

    cl_menu.update_menu('z.t', 30)
    
    
    cl_menu.print_menu()
    
    
    cl_menu.delete_menu('z.x')


    cl_menu.print_menu()

    cl_menu.delete_menu('z')

    cl_menu.print_menu()
    
    
    root_menu=xml_to_menu('test_menu.xml')
    root_menu.print_menu()

    print root_menu.attr_dict
    
    f=root_menu['P1.C1.fdoy']
    doy=f(2010, 5, 1)
    print 'doy', doy
    
    
    print '---text file test ---------'
    
    root_menu2=txt_to_menu('../instrument/vob_menu.txt')
    root_menu2.print_menu()
    
