#include <gtk/gtk.h>
#include <glib.h>
#include "run_anmmpi.h"
#include "tell_user.h"
#include "pdbutil_hiv.h"
#include "miscutil_hiv.h"
#include "readcommandline.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* Inputs to gui:
   Hessian calculation:
   -c      Specify cutoff distance
   -k      Specify force constant
   -p      Print projection matrix to file
   -h      Print block Hessian matrix to file

   Decomposition:
   -n      Specify number of modes
   --lapack Use LAPACK instead of BLZPACK for decomposition
   -?      Specify number of processors for BLZPACK

   Perturbation scanning:
   -s      Specify perturbation scaling factor
   -ko     Specify knockouts

   General:
   Specify output directory
   Specify file format
*/

/* Is this GUI or command line? */
int g_is_gui = 1;
GtkTextView *textview_com;

/* gtk objects */
typedef struct {
  GtkWidget *window_main;
  GtkWidget *dialog_filechooser;
  GtkWidget *checkbutton_is_blockfile;
  GtkWidget *checkbutton_print_projection_matrix;
  GtkWidget *checkbutton_print_block_hessian;
  GtkWidget *checkbutton_calculate_modes;
  GtkEntry *entry_filename;
  GtkEntry *entry_cutoff;
  GtkEntry *entry_gamma;
  GtkEntry *entry_num_modes;
  GtkEntry *entry_num_processors;
  GtkLabel *status_text;

  Local_CL_Opt params;
} gui_vars;


char *InputToCommand(const Local_CL_Opt *CL, int num_processors)
{
  static char *garp;
  char *bunk;
  
  garp = calloc(4096, sizeof(char));
  bunk = calloc(99, sizeof(char));
  strcpy(garp, "mpirun -n ");
  sprintf(bunk, "%d anmmpi_cl ", num_processors);
  strcat(garp, bunk);
  strcat(garp, CL->infile);
  if(g_anm_cutoff != g_anm_cutoff_default){
    strcat(garp, " -c ");
    sprintf(bunk, "%.3lf", g_anm_cutoff);
    strcat(garp, bunk);
  }
  if(CL->num_modes != g_num_modes_default){
    strcat(garp, " -n ");
    sprintf(bunk, "%d", CL->num_modes);
    strcat(garp, bunk);
  }
  if(g_anm_eta != g_anm_eta_default){
    strcat(garp, " -s ");
    sprintf(bunk, "%.3lf", g_anm_eta);
    strcat(garp, bunk);
  }
  if(g_anm_gamma != g_anm_gamma_default){
    strcat(garp, " -k ");
    sprintf(bunk, "%.3lf", g_anm_gamma);
    strcat(garp, bunk);
  }
  if(CL->print_prj_mtx != 0) strcat(garp, " -p");
  if(CL->print_block_hessian != 0) strcat(garp, " -h");
  if(CL->print_vectors != 1) strcat(garp, " -q");
  if(CL->is_blzpack != 1) strcat(garp, " --lapack");
  garp[strlen(garp)]='\0';
  free(bunk);

  return garp;
}

/* 'run_job' executes the given input command on the system and redirects stdout to the GUI textview.
   The input string should include redirection of stdout and stderr. See 'InputToCommand'. Standard
   output will be redirected to the textview object indicated by 'textwindow'. */
void run_job(gchar *command, GtkTextView *textwindow){  
  gchar *output_txt, *error_txt;
  GError *error = NULL;
  GtkTextBuffer *buffer;
  GtkTextMark *mark;
  GtkTextIter iter;

  /* TODO: Run this asynchronously, so that the textbox updates as the program runs. */
  g_spawn_command_line_sync(command, &output_txt, &error_txt, NULL, &error);
  if (error != NULL){
    tell_user("\n**Error in 'run_job':\n");
    tell_user("%s\n",error->message);
    return;
  }

  buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(textwindow));
  mark = gtk_text_buffer_get_insert(buffer);
  gtk_text_buffer_get_iter_at_mark(buffer, &iter, mark);
  gtk_text_buffer_insert(buffer, &iter, output_txt, -1);
  mark = gtk_text_buffer_get_insert(buffer);
  gtk_text_buffer_get_iter_at_mark(buffer, &iter, mark);
  gtk_text_buffer_insert(buffer, &iter, error_txt, -1);
  gtk_text_view_scroll_to_mark(GTK_TEXT_VIEW(textwindow), mark, 0.0, TRUE, 0.0, 0.0);
  return;
}

const G_MODULE_EXPORT set_input_filename(GtkButton *button, gui_vars *input)
{
  input->params.infile = gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(input->dialog_filechooser));
  gtk_entry_set_text(GTK_ENTRY(input->entry_filename), input->params.infile);
}

const G_MODULE_EXPORT on_checkbutton_calculate_modes_toggled(GtkButton *button, gui_vars *input)
{
  if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(input->checkbutton_calculate_modes)))
    gtk_entry_set_text(GTK_ENTRY(input->entry_num_modes), "20");
  else
    gtk_entry_set_text(GTK_ENTRY(input->entry_num_modes), "0");
}




/* This function sends input from the GUI to run an ANM job. */
const G_MODULE_EXPORT on_bt_run_clicked(GtkButton *button, gui_vars *input)
{
  const char *entry_string;
  char *command_line;
  int num_processors;

  /* Set variables to default values */
  input->params.nko = 0;
  input->params.print_prj_mtx = 0;
  input->params.print_block_hessian = 0;
  input->params.add_cypa = 0;
  input->params.is_blzpack = 1;

  /* Get the name of the file to process. */
  input->params.infile = gtk_entry_get_text(input->entry_filename);

  /*
  if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(input->checkbutton_is_blockfile))) 
    input->params.is_blocks = 1;
  else input->params.is_blocks = 0;
  */
  if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(input->checkbutton_print_projection_matrix))) 
    input->params.print_prj_mtx = 1;
  else input->params.print_prj_mtx = 0;
  if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(input->checkbutton_print_block_hessian))) 
    input->params.print_block_hessian = 1;
  else input->params.print_block_hessian = 0;
  if(gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(input->checkbutton_calculate_modes))) 
    input->params.print_vectors = 1;
  else input->params.print_vectors = 0;

  /* Cutoff */
  entry_string = gtk_entry_get_text(input->entry_cutoff);
  sscanf(entry_string, "%lf", &g_anm_cutoff);

  /* Gamma */
  entry_string = gtk_entry_get_text(input->entry_gamma);
  sscanf(entry_string, "%lf", &g_anm_gamma);

  /* Number of modes */
  entry_string = gtk_entry_get_text(input->entry_num_modes);
  sscanf(entry_string, "%d", &input->params.num_modes);

  /* Number of processors */
  entry_string = gtk_entry_get_text(input->entry_num_processors);
  sscanf(entry_string, "%d", &num_processors);

  command_line = InputToCommand(&input->params, num_processors);
  tell_user("Command is '%s'\n",command_line);
  run_job(command_line, textview_com);
}

void print_to_textview(const char *text, GtkTextView *textwindow)
{
  GtkTextBuffer *buffer;
  GtkTextMark *mark;
  GtkTextIter iter;

  buffer = gtk_text_view_get_buffer(GTK_TEXT_VIEW(textwindow));
  mark = gtk_text_buffer_get_insert(buffer);
  gtk_text_buffer_get_iter_at_mark(buffer, &iter, mark);
  gtk_text_buffer_insert(buffer, &iter, text, -1);

  return;
}

int tell_user(const char *fmt,...)
{
  int n;
  int size = 1024;     /* Guess we need no more than 1024 bytes */
  char *p, *np;
  va_list ap;

  if(rank==0){

    if((p = malloc(size)) == NULL) return 1;

    while (1) {

      /* Try to print in the allocated space */
      va_start(ap, fmt);
      n = vsnprintf(p, size, fmt, ap);
      va_end(ap);

      /* Check error code */
      if (n < 0) return 1;

      /* If that worked, print the string to the appropriate output. */
      if (n < size){
	if(g_is_gui == 1) print_to_textview(p, textview_com);
	else{
	  printf("%s",p);
	  fflush(stdout);
	}
	free(p);
	return 0;
      }
      /* Else try again with more space */
      size = n + 1;       /* Precisely what is needed */
      if ((np = realloc (p, size)) == NULL) {
	free(p);
	return 1;} 
      else
	p = np;
    }
  }
  return 0;
}


int main(int argc, char *argv[])
{
  GtkBuilder *builder;
  GError *error = NULL;
  gui_vars gui;
  gchar *output_txt, *error_txt;

  g_anm_cutoff = g_anm_cutoff_default;
  g_anm_gamma = g_anm_gamma_default;
  g_anm_eta = g_anm_eta_default;
  g_bound_resid = g_bound_resid_default;

  /* TODO: Add a switch that will determine whether to call the CLI or GUI */
  if(argc>1){
    g_is_gui = 0;
    ReadCommandLine(argc, argv, &(gui.params));
    run_anmmpi(&gui.params);
    if(gui.params.nko>0) free_ivector(gui.params.KO, 1, gui.params.nko);
  }
  else{
    gtk_init(&argc, &argv);

    /* Construct a GtkBuilder instance and load UI description */
    builder = gtk_builder_new_from_file("anmmpi.glade");
    /*
    builder = gtk_builder_new ();
    if(!gtk_builder_add_from_file(builder, "anmmpi.glade", &error)){
      g_warning("%s\n",error->message);
      g_free(error);
      return 1;
    }
    */

    /* Introduce the widgets to the program. This only needs to be done for widgets that 
       pass information to/from the main program. Look and feel of GUI is handled in glade. */
    gui.window_main = GTK_WIDGET(gtk_builder_get_object(builder, "window_main"));
    gui.dialog_filechooser = GTK_WIDGET(gtk_builder_get_object(builder, "dialog_filechooser"));
    gui.checkbutton_is_blockfile = GTK_WIDGET(gtk_builder_get_object(builder, "checkbutton_is_blockfile"));
    gui.entry_filename = GTK_ENTRY(gtk_builder_get_object(builder, "entry_filename"));
    gui.entry_cutoff = GTK_ENTRY(gtk_builder_get_object(builder, "entry_cutoff"));
    gui.entry_gamma = GTK_ENTRY(gtk_builder_get_object(builder, "entry_gamma"));
    gui.entry_num_processors = GTK_ENTRY(gtk_builder_get_object(builder, "entry_num_processors"));
    gui.checkbutton_print_projection_matrix = GTK_WIDGET(gtk_builder_get_object(builder, "checkbutton_print_projection_matrix"));
    gui.checkbutton_print_block_hessian = GTK_WIDGET(gtk_builder_get_object(builder, "checkbutton_print_block_hessian"));
    gui.checkbutton_calculate_modes = GTK_WIDGET(gtk_builder_get_object(builder, "checkbutton_calculate_modes"));
    gui.entry_num_modes = GTK_ENTRY(gtk_builder_get_object(builder, "entry_num_modes"));
    textview_com = GTK_TEXT_VIEW(gtk_builder_get_object(builder, "textview_com"));
    gui.status_text = GTK_LABEL(gtk_builder_get_object(builder, "status_text"));

    /* Connect signal handlers to the constructed widgets. */
    gtk_builder_connect_signals(builder, &gui);
    g_object_unref(G_OBJECT(builder));
    g_spawn_command_line_sync("nproc --all", &output_txt, &error_txt, NULL, &error);
    output_txt[strlen(output_txt)-1]='\0';
    gtk_entry_set_text(GTK_ENTRY(gui.entry_num_processors), output_txt);
    gtk_widget_show(gui.window_main);
    gtk_main();


  }
  return 0;
}
