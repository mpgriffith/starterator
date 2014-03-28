#!/usr/bin/env python
import MySQLdb
from gi.repository import Gtk
import ConfigParser
import os
import subprocess
import sys
import utils
from uiStarterate import StarteratorEnterInformation
from uiDialogs import DatabaseInfoDialog, PreferencesDialog
import time
import phamgene

"""
GUI Application for Starterator
Process:
    Open Starterator, StarteratorWindow appears
    Choose one of 5 choices, a new dialog pops up based on that choice
    Enter in relevant information, press Starterate
        -Closing this dialog when Starterator is running will cause starterator to stop some time after exit
        -Not instantly, but it checks pretty often-depends on what is happening
    When done, new dialog shows with a link to final report
    If error occurs, dialog shows up
    If unable to connect to the database, DatabaseInfoDialog appears until correct info added
"""

MENU_UI = """
<ui>
  <menubar name='MenuBar'>
    <menu action='File'>
        <menuitem action="ViewReports" />
      <menuitem action='FileQuit' />
    </menu>
    <menu action='Edit'>
      <menuitem action='Preferences' />
    </menu>
    <menu action='Help'>
      <menuitem action='Contents'/>
      <separator />
      <menuitem action='About'/>
    </menu>
  </menubar>
</ui>
"""


class StarteratorWindow(Gtk.Window):
    def show_choice_button(self):
        self.vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
        self.add(self.vbox)
        hbox = Gtk.Box(spacing= 6)
        self.choice_list = Gtk.ListStore(str)
        for item in self.choices:
            self.choice_list.append([item])
        self.choice_label = Gtk.Label('Choose what you want to starterate: ')

        self.choice_box = Gtk.ComboBox().new_with_model(self.choice_list)
        self.choice_button = Gtk.Button("OK")
        renderer_text = Gtk.CellRendererText()
        self.choice_box.pack_start(renderer_text, True)
        self.choice_box.add_attribute(renderer_text, "text", 0)
        self.choice_button.connect('clicked', self.on_choice_changed)
        self.vbox.pack_start(self.choice_label, False, False, 0)
        hbox.pack_start(self.choice_box, False, False, 0)
        hbox.pack_start(self.choice_button, False, False, 0)
        self.vbox.pack_start(hbox, False, False, 0)

  
    def __init__(self):
        Gtk.Window.__init__(self, title="Starterator")

        self.choices = ['Whole Phamerated Phage', 'Whole Unphamerated Phage', 
                        'One Phamerated Gene', 'One Unphamerated Gene', 'Pham']
  
        # self.set_border_width(10)
        self.config_info = utils.get_config()
        self.set_icon_from_file(utils.icon_file)
        self.vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        self.add(self.vbox)
        action_group = Gtk.ActionGroup("my_actions")

        self.add_file_menu_actions(action_group)
        self.add_edit_menu_actions(action_group)
        self.add_help_menu_actions(action_group)

        uimanager = self.create_ui_manager()
        uimanager.insert_action_group(action_group)

        menubar = uimanager.get_widget("/MenuBar")
        box = Gtk.Box(orientation=Gtk.Orientation.VERTICAL)
        box.pack_start(menubar, False, False, 0)
        self.vbox.add(box)
        vbox = Gtk.Box(orientation= Gtk.Orientation.VERTICAL, spacing=6)
        hbox = Gtk.Box(spacing= 6)
        self.choice_list = Gtk.ListStore(str)
        for item in self.choices:
            self.choice_list.append([item])
        self.choice_label = Gtk.Label('Choose what you want to starterate: ')
        self.choice_box = Gtk.ComboBox().new_with_model(self.choice_list)
        renderer_text = Gtk.CellRendererText()
        self.choice_box.pack_start(renderer_text, True)
        self.choice_box.add_attribute(renderer_text, "text", 0)
        self.choice_button = Gtk.Button("OK")
        self.choice_button.connect('clicked', self.on_choice_changed, self.choice_box)
        self.vbox.pack_start(self.choice_label, False, False, 0)
        hbox.pack_start(self.choice_box, True, True, 0)
        hbox.pack_start(self.choice_button, True, True, 0)
        vbox.add(hbox)
        self.vbox.pack_start(vbox, True, True, 0)
        self.vbox.show()
        # self.config_info = {}
        # print utils.get_configuration
        self.check_blast_type()
        db = self.attempt_db_connect()


    def create_ui_manager(self):
        uimanager = Gtk.UIManager()

        # Throws exception if something went wrong
        uimanager.add_ui_from_string(MENU_UI)

        # Add the accelerator group to the toplevel window
        accelgroup = uimanager.get_accel_group()
        self.add_accel_group(accelgroup)
        return uimanager

    def add_file_menu_actions(self, action_group):
        action_filemenu = Gtk.Action("File", "File", None, None)
        action_group.add_action(action_filemenu)

        action_fileview = Gtk.Action("ViewReports", "View Completed Reports", None, None)
        action_fileview.connect("activate", self.on_view_clicked)
        action_group.add_action(action_fileview)

        action_filequit = Gtk.Action("FileQuit", None, None, Gtk.STOCK_QUIT)
        action_filequit.connect("activate", Gtk.main_quit)
        action_group.add_action(action_filequit)

    def add_help_menu_actions(self, action_group):
        action_menu = Gtk.Action("Help", "Help", None, None)
        action_group.add_action(action_menu)

        action_contents = Gtk.Action("Contents", "Contents", None, None)
        action_contents.connect("activate", self.on_contents_clicked)
        action_group.add_action(action_contents)

        action_contents = Gtk.Action("About", "About", None, None)
        action_contents.connect("activate", self.on_about_clicked)
        action_group.add_action(action_contents)

    def add_edit_menu_actions(self, action_group):
        action_menu = Gtk.Action("Edit", "Edit", None, None)
        action_group.add_action(action_menu)
        action_contents = Gtk.Action("Preferences", "Preferences", None, None)
        action_contents.connect("activate", self.on_pref_clicked)
        action_group.add_action(action_contents)
    
    def on_view_clicked(self, button):
        print os.path.abspath(self.config_info["final_file_dir"])
        os.system("xdg-open \""+ os.path.abspath(self.config_info["final_file_dir"])+"\"")

    def on_contents_clicked(self, button):
        root_dir = utils.help_files 
        print root_dir
        uri = 'ghelp:%s/index.page' % root_dir
        print uri
        time_now = int(time.time())

        Gtk.show_uri(None, uri, time_now)

    def on_about_clicked(self, button):
        pass

    def on_pref_clicked(self, data):
        dialog = PreferencesDialog(self, self.config_info)
        response = dialog.run()
        if response == Gtk.ResponseType.APPLY:
            utils.write_to_config_file(self.config_info)
        dialog.destroy()

    def check_blast_type(self):
        blast_dir = self.config_info['blast_dir']
        blast_cmd = self.config_info['blast_dir'] + 'makeblastdb'
        print os.path.join(blast_dir,'blastp')
        try:
            print blast_cmd
            subprocess.check_call([blast_cmd, '-help'])
            print blast_cmd, "worked"
            self.config_info['legacy_blast'] = False
        except:
            try:
                print os.path.join(blast_dir,'formatdb')
                self.config_info['legacy_blast'] = True
            except:
                dialog = Gtk.MessageDialog(self, 0, Gtk.MessageType.ERROR,
                    Gtk.ButtonsType.CANCEL, "BLAST is not installed")
                dialog.format_secondary_text(
                "Please install BLAST.")
                dialog.run()
                dialog.destroy()
                sys.exit(1)
                
    
    def get_configuration(self):
        config = ConfigParser.RawConfigParser()
        print utils.STARTERATOR_PATH + "/extras/starterator.config"
        config.read(utils.STARTERATOR_PATH + "/extras/starterator.config")
        print config
        self.config_info = dict(config.items('Starterator'))
        # self.config_info['intermediate_file_dir'] = os.path.abspath(self.config_info['intermediate_file_dir'])+ '/'
        # self.config_info['final_file_dir'] = os.path.abspath(self.config_info['final_file_dir']) + '/'
        # self.config_info['protein_db'] = os.path.abspath(self.config_info['protein_db']) + '/'
      

    def show_folders(self, name, readable_name):
        dialog = Gtk.FileChooserDialog("Please choose the folder for storing %s" %readable_name, self,
            Gtk.FileChooserAction.CREATE_FOLDER,
            (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
             Gtk.STOCK_OPEN, Gtk.ResponseType.OK))
        response = dialog.run()
        if response == Gtk.ResponseType.OK:
            print "Open clicked"
            print "File selected: " + dialog.get_filename()
            self.config_info[name] = dialog.get_filename()
            print name, self.config_info[name]
        elif response == Gtk.ResponseType.CANCEL:
            print "Cancel clicked"
        dialog.destroy()

    def db_connect(self):
        db = MySQLdb.connect(self.config_info['database_server'], 
                self.config_info['database_user'],
                self.config_info['database_password'],
                self.config_info['database_name'])
        return db
    
    def attempt_db_connect(self):
        try:
            print 'attempting to connect', self.config_info
            db = MySQLdb.connect(self.config_info['database_server'], 
                    self.config_info['database_user'],
                    self.config_info['database_password'],
                    self.config_info['database_name'])
            db.close()
        except:
            db_dialog = DatabaseInfoDialog(self, self.config_info['database_server'], 
                    self.config_info['database_name'],
                    self.config_info['database_user'])
            response = db_dialog.run()
            if response == Gtk.ResponseType.OK:
                self.config_info['database_server'] = db_dialog.info['db_server'] 
                self.config_info['database_user'] = db_dialog.info['db_user'] 
                self.config_info['database_password'] = db_dialog.info['db_password']
                self.config_info['database_name'] = db_dialog.info['db_name'] 
            db_dialog.destroy()
            self.attempt_db_connect()
        else:
            utils.write_to_config_file(self.config_info)

    
    def clean_up_files(self):
        for f in os.listdir(self.config_info['intermediate_file_dir']):
            os.remove(os.path.join(self.config_info['intermediate_file_dir'], f))

    def check_files(self):
        if self.config_info["intermediate_file_dir"] == "?":
            self.show_exception("Intermediate files")
            return False
        if self.config_info["final_file_dir"] == "?":
            self.show_exception("Final Reports")
            return False
        # if self.config_info["protein_db"] == "?":
        #     self.show_exception("Protein Database")
        #     return False
        return True
    
    def show_exception(self, name):
        dialog = Gtk.MessageDialog(self, 0, Gtk.MessageType.ERROR,
            Gtk.ButtonsType.OK, "%s folder is not indicated!" % name)
        dialog.format_secondary_text("Please choose folder in the Preferences.")
        dialog.run()
        dialog.destroy()


    def on_choice_changed(self, button, combo):
        # if self.info_entry_box:
        #     del self.info_entry_box
        if not self.check_files():
            return
        self.config_info['intermediate_file_dir'] = os.path.abspath(self.config_info['intermediate_file_dir'])+ '/'
        self.config_info['final_file_dir'] = os.path.abspath(self.config_info['final_file_dir']) + '/'
        self.config_info['protein_db'] = os.path.abspath(self.config_info['protein_db']) + '/'
        db = self.db_connect()
        phamgene.check_protein_db(self.config_info["count"])
        tree_iter = combo.get_active_iter()
        if tree_iter != None:
            model = combo.get_model()
            choice = model[tree_iter][0]
            dialog = StarteratorEnterInformation(self, choice, self.config_info)