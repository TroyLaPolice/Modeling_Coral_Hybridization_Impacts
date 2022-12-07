# Agent-Based Evolutionary Simulations to Infer Coral Hybridization Impacts on Population Genetic Diversity and Adaptability




#### File descriptions:

###### coral_model.slim
This is the current working model used to run the simulations

###### automated_jobsub.py
This python script is run to output jub submission scripts that can be submitted to the cluster. It requires a formatted text file with the parameter inputs. It parses this file and builds job submission scripts based on the parameter combinations provided. It then generates a final master file that submits all jobs at once. 

###### coral_model_spatial.slim
This is a test version of the model trialed to simulate local adaptation and habitat preference. 

###### laughing_bird_caye.png
This is the map file used to simulate realistic reef landscape on the caye

###### simple.png
This is the simplified map file
