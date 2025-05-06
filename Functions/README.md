Functions in this folder work with data in the following format.

Data from a single mouse is represented in a N_samples x (N_neurons+4) data matrix. N_samples is the total number of samples over all trials and all sessions for a given mouse. N_neurons is the number of neurons.

Each of the first N_neurons columns contain the calcium fluorescence data (dF/F). The last four columns respectively contain the following data.
1. The external variable. For PPC data (Driscoll et al.), this is position or heading in the maze. For V1 data, this is the grating orientation and the time bin of the movie.
2. The trial index.
3. The trial type. For PPC, values of 2 and 3 correspond to left and right turns, respectively. For V1, values of 2 and 3 correspond to grating and movie stimuli.
4. The session number.
