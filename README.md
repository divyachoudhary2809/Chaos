# Chaos in bacterial stress responses
The folder contains 4 folders, three of which contain code to run the simulations for the following:<br>
(1) Oxidative stress response : Oxidative_stress<br>
(2) General stress response : General_stress<br>
(3) Oxidative stress response with added noise : Noise_model <br>

and the last folder is for experimental data analysis:<br>
(4) Experimental data analysis: <br>

Simulation folders contain the following python files: <br>
(1)Oxidative_stress contains -<br>
	&emsp;&emsp;&emsp;&emsp;(a) imagesFromModel_function_1001.py<br>
	&emsp;&emsp;&emsp;&emsp;(b) plotting.py<br>
	&emsp;&emsp;&emsp;&emsp;(c) response_model.py<br>
	&emsp;&emsp;&emsp;&emsp;(d) main_file.py<br>
	&emsp;&emsp;&emsp;&emsp;(e) run_oxidative.py<br>
(2) General_stress contains - <br>
	&emsp;&emsp;&emsp;&emsp;(a) imagesFromModel_function_1001.py<br>
	&emsp;&emsp;&emsp;&emsp;(b) plotting.py<br>
	&emsp;&emsp;&emsp;&emsp;(c) abstract_response_model.py<br>
	&emsp;&emsp;&emsp;&emsp;(d) main_file.py<br>
	&emsp;&emsp;&emsp;&emsp;(e) run_general.py<br>
(3) Noise_model contains - <br>
	&emsp;&emsp;&emsp;&emsp;(a) discreteresponsemodel.py<br>
	&emsp;&emsp;&emsp;&emsp;(b) plotting.py<br>
	&emsp;&emsp;&emsp;&emsp;(c) main_file_noisy.py<br>
	&emsp;&emsp;&emsp;&emsp;(d) run_noise_oxidative_stress_model.py<br>

Here, <br>
(a) The 'run_oxidative.py', 'run_general.py' and 'run_noise_oxidative_stress_model.py' are codes to run the simulations for oxidative, general stress response or noisy oxidative stress response model respectively based on user desired inputs.<br>
(b) The 'imagesFromModel_function_1001.py' file contains codes for coding images to collate as movies in Fiji for response outputs from the model. The images are saved in an output folder that the code makes 'Output_images' in the folder where the code is being run.<br>
(c) The 'plotting.py' contains functions to plot different variables as user desires. By default it plots the GrxA (for oxidative stress) / enzyme (for general stress) dynamics over time and the poincare plots for data points 400 minutes after treatment time.<br>
(d) The 'response_model.py', 'abstract_response_model.py' or 'discreteresponsemodel.py' files contain the Ordinary differential equations (ODEs) for the oxidative and general stress response respectively.<br>
(e) Finally, the 'main_file.py' or 'main_file_noisy.py' contains the code to run all 3 models together -- the growth rate model + the cell-cell interaction model + the stress response model.<br>

(4) Experimental data analysis contains - <br>
	&emsp;&emsp;&emsp;&emsp;(a) data_analysis.py<br>
Here, <br>
The data_analysis.py file contains code to encode for the following to produce the experimental data figures. The RAW experimental data files input to the code are uploaded on Oxford Research Archive (https://doi.org/10.5287/ora-b7dw9pmqd):<br>
(a) Plot GrxA fluorescence expression traces for individual mother cells growing in microfluidic growth channels <br>
(b) Mean elongation rate and PgrxA expression for cells at different positions in the growth trench  <br>
(c) Mean PgrxA autocorrelation (DPgrxA ACF) for all mother cells and DPgrxA ACF for individual cells  <br>
(d) Poincare plots and PgrxA expression traces for alive and dead cells  <br>
(e) Interdivision time versus DPgrxA ACF peak time for individual cells <br>
 

Refer to 'Chaos in a bacterial stress response' Divya Choudhary, Kevin R Foster*, Stephan Uphoff* for more details about the models and data collection + analysis of experimental data.<br>
*=Corresponding author
