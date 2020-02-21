 
# This file is auto-generated from a Python script that parses a PhysiCell configuration (.xml) file.
#
# Edit at your own risk.
#
import os
from ipywidgets import Label,Text,Checkbox,Button,HBox,VBox,FloatText,IntText,BoundedIntText,BoundedFloatText,Layout,Box
    
class UserTab(object):

    def __init__(self):
        
        micron_units = Label('micron')   # use "option m" (Mac, for micro symbol)

        constWidth = '180px'
        tab_height = '500px'
        stepsize = 10

        #style = {'description_width': '250px'}
        style = {'description_width': '25%'}
        layout = {'width': '400px'}

        name_button_layout={'width':'25%'}
        widget_layout = {'width': '15%'}
        units_button_layout ={'width':'15%'}
        desc_button_layout={'width':'45%'}

        param_name1 = Button(description='tumor_transition_rate', disabled=True, layout=name_button_layout)
        param_name1.style.button_color = 'lightgreen'

        self.tumor_transition_rate = FloatText(
          value=0.0022956841138659324,
          step=0.0001,
          style=style, layout=widget_layout)

        param_name2 = Button(description='tumor_max_necrosis_rate', disabled=True, layout=name_button_layout)
        param_name2.style.button_color = 'tan'

        self.tumor_max_necrosis_rate = FloatText(
          value=0.002777777777777778,
          step=0.0001,
          style=style, layout=widget_layout)

        param_name3 = Button(description='elastic_rate', disabled=True, layout=name_button_layout)
        param_name3.style.button_color = 'lightgreen'

        self.elastic_rate = FloatText(
          value=0.05 ,
          step=0.01,
          style=style, layout=widget_layout)

        param_name4 = Button(description='plastic_rate', disabled=True, layout=name_button_layout)
        param_name4.style.button_color = 'tan'

        self.plastic_rate = FloatText(
          value=0.0005,
          step=0.0001,
          style=style, layout=widget_layout)

        param_name5 = Button(description='max_ECM_displacement', disabled=True, layout=name_button_layout)
        param_name5.style.button_color = 'lightgreen'

        self.max_ECM_displacement = FloatText(
          value=0.75,
          step=0.1,
          style=style, layout=widget_layout)

        param_name6 = Button(description='tumor_max_pressue', disabled=True, layout=name_button_layout)
        param_name6.style.button_color = 'tan'

        self.tumor_max_pressue = FloatText(
          value=1.0,
          step=0.1,
          style=style, layout=widget_layout)

        units_button1 = Button(description='1/min', disabled=True, layout=units_button_layout) 
        units_button1.style.button_color = 'lightgreen'
        units_button2 = Button(description='1/min', disabled=True, layout=units_button_layout) 
        units_button2.style.button_color = 'tan'
        units_button3 = Button(description='1/min', disabled=True, layout=units_button_layout) 
        units_button3.style.button_color = 'lightgreen'
        units_button4 = Button(description='1/min', disabled=True, layout=units_button_layout) 
        units_button4.style.button_color = 'tan'
        units_button5 = Button(description='micron', disabled=True, layout=units_button_layout) 
        units_button5.style.button_color = 'lightgreen'
        units_button6 = Button(description='', disabled=True, layout=units_button_layout) 
        units_button6.style.button_color = 'tan'

        desc_button1 = Button(description='The proliferation rate of tumor cells', disabled=True, layout=desc_button_layout) 
        desc_button1.style.button_color = 'lightgreen'
        desc_button2 = Button(description='The maximum necrosis rate of tumor cells', disabled=True, layout=desc_button_layout) 
        desc_button2.style.button_color = 'tan'
        desc_button3 = Button(description='The elastic force parameter for parenchyma', disabled=True, layout=desc_button_layout) 
        desc_button3.style.button_color = 'lightgreen'
        desc_button4 = Button(description='The plastic reorganization force parameter for parenchyma', disabled=True, layout=desc_button_layout) 
        desc_button4.style.button_color = 'tan'
        desc_button5 = Button(description='The maximum mechanical displacement for parenchyma', disabled=True, layout=desc_button_layout) 
        desc_button5.style.button_color = 'lightgreen'
        desc_button6 = Button(description='The maximum pressure threshold for tumor cells proliferation', disabled=True, layout=desc_button_layout) 
        desc_button6.style.button_color = 'tan'

        row1 = [param_name1, self.tumor_transition_rate, units_button1, desc_button1] 
        row2 = [param_name2, self.tumor_max_necrosis_rate, units_button2, desc_button2] 
        row3 = [param_name3, self.elastic_rate, units_button3, desc_button3] 
        row4 = [param_name4, self.plastic_rate, units_button4, desc_button4] 
        row5 = [param_name5, self.max_ECM_displacement, units_button5, desc_button5] 
        row6 = [param_name6, self.tumor_max_pressue, units_button6, desc_button6] 

        box_layout = Layout(display='flex', flex_flow='row', align_items='stretch', width='100%')
        box1 = Box(children=row1, layout=box_layout)
        box2 = Box(children=row2, layout=box_layout)
        box3 = Box(children=row3, layout=box_layout)
        box4 = Box(children=row4, layout=box_layout)
        box5 = Box(children=row5, layout=box_layout)
        box6 = Box(children=row6, layout=box_layout)

        self.tab = VBox([
          box1,
          box2,
          box3,
          box4,
          box5,
          box6,
        ])

    # Populate the GUI widgets with values from the XML
    def fill_gui(self, xml_root):
        uep = xml_root.find('.//microenvironment_setup')  # find unique entry point
        vp = []   # pointers to <variable> nodes
        if uep:
            for var in uep.findall('variable'):
                vp.append(var)

        uep = xml_root.find('.//user_parameters')  # find unique entry point
        self.tumor_transition_rate.value = float(uep.find('.//tumor_transition_rate').text)
        self.tumor_max_necrosis_rate.value = float(uep.find('.//tumor_max_necrosis_rate').text)
        self.elastic_rate.value = float(uep.find('.//elastic_rate').text)
        self.plastic_rate.value = float(uep.find('.//plastic_rate').text)
        self.max_ECM_displacement.value = float(uep.find('.//max_ECM_displacement').text)
        self.tumor_max_pressue.value = float(uep.find('.//tumor_max_pressue').text)


    # Read values from the GUI widgets to enable editing XML
    def fill_xml(self, xml_root):
        uep = xml_root.find('.//microenvironment_setup')  # find unique entry point
        vp = []   # pointers to <variable> nodes
        if uep:
            for var in uep.findall('variable'):
                vp.append(var)

        uep = xml_root.find('.//user_parameters')  # find unique entry point
        uep.find('.//tumor_transition_rate').text = str(self.tumor_transition_rate.value)
        uep.find('.//tumor_max_necrosis_rate').text = str(self.tumor_max_necrosis_rate.value)
        uep.find('.//elastic_rate').text = str(self.elastic_rate.value)
        uep.find('.//plastic_rate').text = str(self.plastic_rate.value)
        uep.find('.//max_ECM_displacement').text = str(self.max_ECM_displacement.value)
        uep.find('.//tumor_max_pressue').text = str(self.tumor_max_pressue.value)
