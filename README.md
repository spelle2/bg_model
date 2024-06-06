# Neurocomputational model of Basal Ganglia to simulate Reversal Learning Taks

This repository contains a model of basal ganglia able to simulate both deterministic and probabilistic reversal learning tasks. Moreover, with this model we were able to 
explore the dopaminergic system's influence on learning, particularly in Parkinson's Disease (PD), where medication leads to impaired adaptability.

## Basal Ganglia model
<p align="center">
  <img src="BG_MODEL.png" alt="fig"/>
</p>

The model implemented is a representation of the human BG.
The whole network comprises several neural units, characterized by first-order low-pass dynamics to reproduce the integrative property of the membrane and a sigmoidal relation for output activity to represent the presence of lower and upper saturation values for neuronal activity.
The model comprises a sensory representation (S) and a motor representation (C) in the cortex. The sensory representation’s neurons correspond to the stimulus presented to the network, while the motor neurons in the cortex encode possible actions. In addition to the previous neural units downstream the cortex, the model includes the striatum, further subdivided into Go and NoGo pathways, the globus pallidus pars externa (Gpe), the globus pallidus pars interna (Gpi), the thalamus (T), all with a neuron count in a 1-to-1 relationship with cortex, the cholinergic interneuron (ChI) and the subthalamic nucleus (STN), modeled as single neurons representing the entire population activity.
The model integrates three primary pathways: direct (via Go neurons), indirect (via NoGo neurons), and hyperdirect (via STN).


(Baston and Ursino, 2015; Schirru et al., 2022).
## Deterministic Tasks
<p align="center">
  <img src="Deterministic task.png" alt="fig"/>
</p>
In the original experiment presented by Cools et al. (Cools et al., 2006), participants (PD-OFF, PD-ON, Control) predicted outcomes in a simulated card game (Fig. 4). Participants were instructed to imagine themselves as a casino boss observing a player during a card game (see Cools et al. 2006 for more details). 
Participants were presented with two images per trial; one of them is highlighted. By pressing a corresponding button, they had to predict the outcome of the player game: a win (green button) or loss (red button) associated with the highlighted image. This setup allows participants to learn these associations over time. During the task, reversal phases were introduced where the previously learned associations were switched. Now, the image previously associated with a win represents a loss and vice versa. This tests the subject’s adaptability and learning under changed conditions. 
Crucially, the participant predicts "a posteriori" whether the player has won or lost, meaning they could not influence the results, only predict what happened.
The task consists of separate blocks, each comprising 120 trials, during which a reversal occurs several times: in detail, a reversal takes place once an appropriate amount of knowledge is achieved. The criterion of knowledge is a predefined number of consecutive correct responses (ranging from 5 to 9, selected randomly) to prevent predictability of the reversal.
