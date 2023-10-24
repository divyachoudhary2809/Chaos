# Chaos in bacterial stress responses
The folder contains 3 folders that contain code to run the simulations for:
(1) Oxidative stress response : Oxidative_stress
(2) General stress response : General_stress
(2) Oxidative stress response with added noise : Noise_model 

Both folders contain the following python files:
(1)Oxidative_stress contains -
	(a) imagesFromModel_function_1001.py
	(b) plotting.py
	(c) response_model.py
	(d) main_file.py
	(e) run_oxidative.py
(2) General_stress contains - 
	(a) imagesFromModel_function_1001.py
	(b) plotting.py
	(c) abstract_response_model.py
	(d) main_file.py
	(e) run_general.py
(3) Noise_model contains - 
	(a) discreteresponsemodel.py
	(b) plotting.py
	(c) main_file_noisy.py
	(d) run_noise_oxidative_stress_model.py

Here, 
(a) The 'run_oxidative.py', 'run_general.py' and 'run_noise_oxidative_stress_model.py' are codes to run the simulations for oxidative, general stress response or noisy oxidative stress response model respectively based on user desired inputs.
(b) The 'imagesFromModel_function_1001.py' file contains codes for coding images to collate as movies in Fiji for response outputs from the model. The images are saved in an output folder that the code makes 'Output_images' in the folder where the code is being run.
(c) The 'plotting.py' contains functions to plot different variables as user desires. By default it plots the GrxA (for oxidative stress) / enzyme (for general stress) dynamics over time and the poincare plots for data points 400 minutes after treatment time.
(d) The 'response_model.py', 'abstract_response_model.py' or 'discreteresponsemodel.py' files contain the Ordinary differential equations (ODEs) for the oxidative and general stress response respectively.
(e) Finally, the 'main_file.py' or 'main_file_noisy.py' contains the code to run all 3 models together -- the growth rate model + the cell-cell interaction model + the stress response model.
