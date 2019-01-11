STEPS


1) CALCULATE EDGE WEIGHTS

  Steps 1.1-1.3 are precomputation steps for 1.4

  1.1) Compute mutex, mutex_mod, mutex_avgmod
  $ python run compute_mutex_values.py
  It generates the files in ./hint/out/edge_weights/pre/

  1.2) Compute cov, cov_mod, cov_avgmod
  $ python run compute_cov_values.py
  It generates the files in ./hint/out/edge_weights/pre/

  1.3) Compute density1, density2, density3
  $ python run compute_density_values.py
  It generates the files in ./hint/out/edge_weights/pre/

  1.4) Compute the edge weights for a specific model
  $ python run compute_edge_weights.py
  It generates the files in ./hint/out/edge_weights/


2) RANDOM WALK

  2.1) Random walk for HotNet2
  $ python random_walk_hotnet.py
  It generates the .npz matrices in ./hint/out/hotnet2_randomwalk/

  2.2) Random walk for our model
  $ python random_walk.py
  It generates the .npz matrices in ./hint/out/our_randomwalk/


3) FIND THE CONNECTED COMPONENTS

  3.1) HotNet2 connected components
  Connected components for hotnet are in ./hint/out/connected_components/hotnet/ taken from the resources provided from the HotNet2 website/

  3.2) Find the connected components
  $ python connected_components.py
  Change the initial value of threshold_start accordingly.
  It generates the connected components in ./hint/out/connected_components/
  
  3.3) Split large modules
  $ python connected_components_isolarge.py
  Check related connected component to determine threshold values
  It generates the connected components in ./hint/out/connected_components_isolarge/


4) EVALUATE THE MODELS

  4.1) Evaluate the connected components separately
  $ python evaluate_connected_components.py
  Creates evaluation files in ./hint/out/evaluation

  4.2) Evaluating of the connected components needed to create the tables/graphs
  $ python evaluate_connected_components_tab.py
  Creates evaluation files in ./hint/out/evaluation_tab

  4.3) Create the score tables
  $ python calculate_scores.py
  Creates score tables in ./hint/out/evaluation_tab

  4.4) Create the weighted score tables with weighted inverse coverage
  $ python calculate_weighted_scores.py
  Creates weighted score tables in ./hint/out/evaluation_tab
  To calculate weighted coverage use calculate_weighted_scores_wcov.py
