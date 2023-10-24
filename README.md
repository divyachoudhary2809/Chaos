# Chaos in bacterial stress responses
The folder contains 3 folders that contain code to run the simulations for:<br>
(1) Oxidative stress response : Oxidative_stress<br>
(2) General stress response : General_stress<br>
(2) Oxidative stress response with added noise : Noise_model <br>

Both folders contain the following python files: <br>
(1)Oxidative_stress contains -<br>
	(a) imagesFromModel_function_1001.py<br>
	(b) plotting.py<br>
	(c) response_model.py<br>
	(d) main_file.py<br>
	(e) run_oxidative.py<br>
(2) General_stress contains - <br>
	(a) imagesFromModel_function_1001.py<br>
	(b) plotting.py<br>
	(c) abstract_response_model.py<br>
	(d) main_file.py<br>
	(e) run_general.py<br>
(3) Noise_model contains - <br>
	(a) discreteresponsemodel.py<br>
	(b) plotting.py<br>
	(c) main_file_noisy.py<br>
	(d) run_noise_oxidative_stress_model.py<br>

Here, <br>
(a) The 'run_oxidative.py', 'run_general.py' and 'run_noise_oxidative_stress_model.py' are codes to run the simulations for oxidative, general stress response or noisy oxidative stress response model respectively based on user desired inputs.<br>
(b) The 'imagesFromModel_function_1001.py' file contains codes for coding images to collate as movies in Fiji for response outputs from the model. The images are saved in an output folder that the code makes 'Output_images' in the folder where the code is being run.<br>
(c) The 'plotting.py' contains functions to plot different variables as user desires. By default it plots the GrxA (for oxidative stress) / enzyme (for general stress) dynamics over time and the poincare plots for data points 400 minutes after treatment time.<br>
(d) The 'response_model.py', 'abstract_response_model.py' or 'discreteresponsemodel.py' files contain the Ordinary differential equations (ODEs) for the oxidative and general stress response respectively.<br>
(e) Finally, the 'main_file.py' or 'main_file_noisy.py' contains the code to run all 3 models together -- the growth rate model + the cell-cell interaction model + the stress response model.<br>


Refer to 'Chaos in a bacterial stress response' Divya Choudhary, Kevin R Foster, Stephan Uphoff for more details about the models.
