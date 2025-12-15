# Import custom scripts
from . import helpers
import logomaker

# Import python modules (visualization)
import matplotlib.pyplot as plt

# Import python modules (data manipulation)
import pandas as pd 
import numpy as np
import os

# Import config information
config = helpers.read_json()




###########################################################################
# Class that generate motif-logo and perform data analysis on DNA sequences
class mlogo():
    """
    Class that generate motif-logo and perform data analysis on DNA sequences.
    The data can be grouped by the "group_by" argument (column) and generate 
    seperate plot for each subject 
    """

    ######################
    # Initiating the class
    def __init__(self,
                 seq_file,
                 group_by:str = None,
                 divide_subject:bool = True,
                 subject_col:str = "subject_id",
                 sequence_col:str = "sequence",
                 germline_col:str = "germline"):
        """
        seq_file : str -> file path / df object of sequencing data.
        group_by : str -> group the dataset by this column. will create number of subplots according
                          to the number of unique values in the grouped by column.
        divide_subject : bool -> Generate different plot for each subject under `subject_id` column.
        subject_col : str -> name of the subject column, defualt is `subject_id`. 
        sequence_col : str -> name of the DNA sequence column, defualt it `sequence`.
        germline_col : str -> name of the germline DNA sequence column, defualt it `germline`.
        """
        arguments_input = [group_by,
                           divide_subject,
                           subject_col,
                           sequence_col,
                           germline_col]
        
        self.init_dict = {i:j for i,j in zip(["group_by", "divide_subject", "subject_col", "sequence_col", "germline_col"], arguments_input)}
        
        
        # Loading sequences dataset
        try:
            if isinstance(seq_file, str):        
                self.seq_df = pd.read_csv(seq_file, index_col=0)
                print(f"> Dataset loaded (seq_file = '{seq_file}')")

            elif isinstance(seq_file, pd.DataFrame):
                 self.seq_df = seq_file
                 print("> Dataset loaded")

        except:
                print(f"> Invalid input, please make sure seq_file argument entred correctly.\n  (Invalid: seq_file = '{seq_file}')")

        # Translating germline and sequence DNA into amino acid sequence
        self.seq_df["sequence_aa"] = self.seq_df["sequence"].apply(helpers.nt_transalte_104) # Somatic sequence
        self.seq_df["germline_aa"] = self.seq_df["germline"].apply(helpers.nt_transalte_104) # Germline sequence
        
    ###############################################################
    # Plotting the motif logo according to the amino acids position
    def motif_logo(self,
                    by : str, # "sequence" for somatic data or "germline" for germline data
                    aa_start: int,
                    aa_end: int,
                    remove_pos : list = None,
                    fig_title : str = "",
                    save_fig : bool = True,
                    save_name : str = "",
                    yaxis_metric : str ="information"):
        """
        by : str -> "sequence" or "germline" motif.
        aa_start : int -> Start of the amino acid motif, relative to the amino acid sequence start.
        aa_end : end -> End of the amino acid motif, relative to the amino acid sequence start.
        fig_title : str -> Title of the output figure.
        remove_pos : list -> list of integers, this amino acid positions will be removed from the plot.
        save_fig : bool -> if True will save the figure into the output folder, as defined in congif.json.
        save_name : str -> Name of the output plot (file name).
        yaxis_metric -> Metric to be used in the plot y-axis: must be `counts`, `probability`, `weight`, or `information`.
                            (https://logomaker.readthedocs.io/en/latest/implementation.html#logomaker.alignment_to_matrix)
        """
            
        # Translating the amino acid sequences
        ##seq_list = self.seq_df[f"{by}_aa"].apply(lambda X : X[aa_start:aa_end])
            
        # 1st step -> Setting initial shape of 1
        # Final shape of the sub-plot will be placed in this variable
        self.subplot_axidx = [1, 1]

        # 2nd step -> getting the shape of the subplot by the values of grouped by / n-subjects
        # Getting unique group by column values (n-cols)
        group_by_value = self.init_dict["group_by"]
        if isinstance(group_by_value, str): 
             gp_unique = self.seq_df[group_by_value].unique() # Unique values in the group by sequence.
             self.subplot_axidx[0] = len(gp_unique) # Number of unique values, will be number of columns.
             cond_gp = [self.seq_df[group_by_value] == i for i in gp_unique] # Conditions of unique values.
            
        # If not group by any column, number of column = 1.
        else:
             cond_gp = [np.full((self.seq_df.shape[0],1), True).flatten()]
            
        # Getting subjects values (n-rows)
        subject_bool, subject_col = self.init_dict["divide_subject"], self.init_dict["subject_col"]
        subj_unique = self.seq_df[subject_col].unique()
        if subject_bool:
             subj_len = len(subj_unique) # Number of subjects
             self.subplot_axidx[1] = subj_len # Number of rows in the subplot
             cond_subj = [self.seq_df[subject_col] == i for i in subj_unique]

        # if subject not to be divided, rows in sub-plot will be equal to 1.
        else:
             cond_subj = [np.full((self.seq_df.shape[0],1), True).flatten()]

        # Condition matrix
        self.cond_matrix = [cond_gp, cond_subj]
        self.cond_index = [[i,j] for i in range(0, len(cond_gp)) for j in range(0,len(cond_subj))]


        # 3rd step -> getting the condition indexes (1d or 2d, depend on shape)
        # Getting the shape of the sub-plots based on the condition shape
        if ([1, 1] == self.subplot_axidx):
            self.cond_matrix = [np.full((self.seq_df.shape[0],1), True).flatten()]
               
        else:
            try:
                reverse_i = {1:0, 0:1}
                i = self.subplot_axidx.index(1)
                self.cond_matrix = self.cond_matrix[reverse_i[i]]

            except:
                self.cond_matrix = [(i & j) for i in cond_gp for j in cond_subj]

                    
        # Figure initiation
        n_cols, n_rows = self.subplot_axidx[0], self.subplot_axidx[1]
        fig, axs = plt.subplots(ncols=n_cols, 
                                nrows=n_rows,
                                squeeze=False, # Making the axes index 2d even if they are one 1d
                                figsize = (n_cols * 7, n_rows * 3))
        
        motif_data = []
        unique_aa = np.delete(np.unique(list(helpers.codon_dic_updated.values())), 0)
        for i, ic in zip(self.cond_index, self.cond_matrix):
            temp_metric = yaxis_metric
            index_col = i[0]
            index_row = i[1] 

            ax = axs[index_row, index_col]

            ### Temp dataframe itiration -> creation of array of sequences according to the condition
            temp_df = self.seq_df[ic]
            temp_seqs = temp_df[f"{by}_aa"].apply(lambda X : X[aa_start : aa_end]).values

            # Visualization code running 
            if __name__ == "__main__": 
                print(temp_df.ab_target.unique(), temp_df.subject_id.unique())

            ### Logomaker sequence aligment and matrix creation
            # Creation of the sequences the Matrix, to_type='probability' are other options.
            aa_counts = logomaker.alignment_to_matrix(sequences=temp_seqs, 
                                                      to_type='counts', 
                                                      characters_to_ignore=".-"
                                                     )
            
            # Chaning the index name to the amino acid positions
            aa_counts.index = list(range(aa_start+1, aa_end+1))
            if isinstance(remove_pos, list):
                try:
                    rem_pos = [int(i) for i in remove_pos]
                    aa_counts = aa_counts[aa_counts.index.isin(rem_pos) == False]
                except:
                    print(f"> Failed to remove amino acid positions {remove_pos}, invalid argument input.")
            
            # Convert the count matrixto different matrix is needed (yaxis_metric argument)
            if (temp_metric != "counts") & (temp_metric in ["probability", "weight", "information"]):
                # There is a bug that prevent the transformation of counts matrix to information in some cases
                try:

                    og_rows = aa_counts[aa_counts.sum(axis=1) == 0]
                    aa_counts = aa_counts[aa_counts.sum(axis=1) != 0]

                    aa_metric = logomaker.transform_matrix(aa_counts, 
                                                           from_type='counts', 
                                                           to_type=yaxis_metric
                                                            )
                    
                    aa_metric = pd.concat([aa_metric, og_rows]).sort_index()

                except:
                        print(f"Error converting the counts df to {temp_metric} - index = [{index_col},{index_row}] type. \nproblomatic dataframe saved to output folder")
                        path_error = f"{config["output_folder"]}\\error_df_index"
                        temp_metric = "counts"

                        if os.path.exists(path_error) is False:
                            os.mkdir(path_error)

                        aa_metric.to_csv(path_error + f"\\motif_counts_[{index_col},{index_row}]_.csv")

            # Preparing motif dataset to be saved  
            temp_data = aa_metric.reindex(columns=unique_aa).fillna(0)
            temp_data.index = temp_data.index.rename(f"[{str(index_row)}, {str(index_col)}")
            temp_data["fig_index"] = f"[{index_row},{index_col}]"
            
            # If grouped by column or subject, add this information column to the dataframe
            if isinstance(group_by_value, str):
                temp_data[group_by_value] = gp_unique[index_col]
            if subject_bool:
                temp_data["subject"] = subj_unique[index_row]

            motif_data.append(temp_data) # temp data
            
            # Revising x-axis values for plotting
            aa_metric.index = list(range(0, len(aa_metric)))

            ### Logomaker sub-plot creation
            logo = logomaker.Logo(aa_metric,
                                  color_scheme='skylign_protein', # 'chemistry' colors by property (polar, acidic, etc.)
                                  vpad=.1,
                                  width=.8,
                                  ax=ax,
                                  )

            ### Sub-plot properties - style the Logo
            # Enabling only the left and bottom subplot lines
            logo.style_spines(visible=False)
            logo.style_spines(spines=['left', 'bottom'], visible=True)

            # Setting lables
            if_info = ""
            if temp_metric == "information":
                 if_info = " (Bits)"
                 logo.ax.set_ylim(0,4.5) # Defining y-limit if metrix = bits

            logo.ax.set_ylabel(f"{temp_metric}{if_info}".capitalize(), fontsize=12) #y-lavbel
            logo.ax.set_xlabel("Amino Acid Position", fontsize=12) #x-label
            logo.ax.set_title(f"", fontsize=14) #subplot label 

            # x-axis labels -> returning to original values of actual aa positions
            aa_index = aa_counts.index.values
            ticks_pos, tick_names = range(0, len(aa_index)), [str(i) for i in aa_index]
            logo.ax.set_xticks(ticks_pos)
            logo.ax.set_xticklabels(tick_names)
          
        ### Global figure properties
        # Set column titles (top row)
        if isinstance(group_by_value, str):
            for ax, col in zip(axs[0,:], gp_unique):
                ax.set_title(col, pad=20, fontsize=15)

        # Set row titles (right side)
        if subject_bool:
            pad = 5 # Padding from the plot edge
            for ax, row in zip(axs[:, -1], subj_unique):
                ax.annotate(row, xy=(1, 0.5), xytext=(pad, 0), fontsize=15,
                            xycoords='axes fraction', textcoords='offset points',
                            ha='left', va='center', rotation=270)
        
        # fig title and layout
        fig.suptitle(fig_title, fontsize=15)
        plt.tight_layout()
        
        # > save figure
        name = f"{str(aa_start+1)+str(aa_end+1)}logo"
        if save_fig:
            output_folder = config["output_folder"]

            if os.path.exists(output_folder) is False:
                os.mkdir(output_folder)
                print(f"{output_folder} folder was created.")
            
            name = f"{save_name}_{str(aa_start+1)+str(aa_end+1)}logo"
            plt.savefig(f"{output_folder}\\motif_figure\\{name}.png",  bbox_inches='tight')
            print(f"Plot saved as {name} at {output_folder} folder")
        
        pd.concat(motif_data).to_csv(f"{output_folder}\\motif_data\\{name}.csv")
        print("Raw motif dataframe information saveed into motif_data folder")

 
        return self.cond_matrix
