# Neurocomputational model of Basal Ganglia to simulate Reversal Learning Taks

This repository contains a model of basal ganglia able to simulate both deterministic and probabilistic reversal learning tasks. Moreover, with this model we were able to 
explore the dopaminergic system's influence on learning, particularly in Parkinson's Disease (PD), where medication leads to impaired adaptability.
[The link to the official paper will be available soon]

## Basal Ganglia model
<p align="center">
  <img src="BG_MODEL.png" alt="fig"/>
</p>

The model implemented is a representation of the human BG.
The whole network comprises several neural units, characterized by first-order low-pass dynamics to reproduce the integrative property of the membrane and a sigmoidal relation for output activity to represent the presence of lower and upper saturation values for neuronal activity.
The model comprises a sensory representation (S) and a motor representation (C) in the cortex. The sensory representation’s neurons correspond to the stimulus presented to the network, while the motor neurons in the cortex encode possible actions. In addition to the previous neural units downstream the cortex, the model includes the striatum, further subdivided into Go and NoGo pathways, the globus pallidus pars externa (Gpe), the globus pallidus pars interna (Gpi), the thalamus (T), all with a neuron count in a 1-to-1 relationship with cortex, the cholinergic interneuron (ChI) and the subthalamic nucleus (STN), modeled as single neurons representing the entire population activity.
The model integrates three primary pathways: direct (via Go neurons), indirect (via NoGo neurons), and hyperdirect (via STN).


## Deterministic Task
<p align="center">
  <img src="Deterministic Task.png" alt="fig2"/>
</p>
Participants were instructed to imagine themselves as a casino boss observing a player during a card game (see Cools et al. 2006 for more details). In particular, they were presented with two images per trial; one of them is highlighted. By pressing a corresponding button, they had to predict the outcome of the player game: a win (green button) or loss (red button) associated with the highlighted image. This setup allows participants to learn these associations over time. During the task, reversal phases were introduced where the previously learned associations were switched. Now, the image previously associated with a win represents a loss and vice versa. This tests the subject’s adaptability and learning under changed conditions. 
Crucially, the participant predicts "a posteriori" whether the player has won or lost, meaning they could not influence the results, only predict what happened.
The task consists of separate blocks, each comprising 120 trials, during which a reversal occurs several times: in detail, a reversal takes place once an appropriate amount of knowledge is achieved. The criterion of knowledge is a predefined number of consecutive correct responses (ranging from 5 to 9, selected randomly) to prevent predictability of the reversal.

## Probabilistic Task
The probabilistic task is the one presented by(Cools et al., 2001), two different colors are directly associated with rewards or punishment in a probabilistic way. In detail, selecting a color resulted in a reward 80% of the time and a punishment 20% of the time. In contrast, choosing the other color is rewarded with a 20% probability and punished with an 80% probability. 


## Citing
```bibtex
@article{Baston2015,
  title = {A Biologically Inspired Computational Model of Basal Ganglia in Action Selection},
  volume = {2015},
  ISSN = {1687-5265, 1687-5273},
  url = {https://doi.org/10.1155/2015/187417},
  DOI = {10.1155/2015/187417},
  journal = {Computational Intelligence and Neuroscience},
  author = {Baston,  Chiara and Ursino,  Mauro},
  year = {2015},
  pages = {1-24}
}
```
```bibtex
@article{Schirru2022,
  title = {Phasic Dopamine Changes and Hebbian Mechanisms during Probabilistic Reversal Learning in Striatal Circuits: A Computational Study},
  volume = {23},
  ISSN = {1422-0067},
  url = {https://doi.org/10.3390/ijms23073452},
  DOI = {10.3390/ijms23073452},
  journal = {International Journal of Molecular Sciences},
  author = {Schirru, Miriam and Véronneau-Veilleux, Florence and Nekka, Fahima  and Ursino,  Mauro},
  year = {2015},
  month = {03}
  pages = {3452}
}
```

