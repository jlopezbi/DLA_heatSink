# Coral Simulation

This is a rhinoPython script that generates coraly shapes.

It is basically diffusion limited aggregation (DLA) on a mesh

To run this get rhino 5, make sure rhinoPython is set up, and then in rhino type "runpythonscript." Navigate to the folder where the script is saved and run meshDLA.py. It will ask to select a 'seed' mesh which will act as the base mesh. The seed mesh will be modified by the algorithm. Be patient, it tends to take a few minutes to run a simulation and the only indiciation it is running is the spinny beach ball of death (mac).  

There is an included 3dm file that you can use for the seed.

This is based on work by Jaap Kaandorp and Roeland Merks.


## Simulation Method
The simulation starts with an arbitrary ‘seed’ composed of points, called nodes, which represent polyps, connected by edges, which represent the physical connections between polyps. This structure is called a mesh. Spheres, representing nutrients, move in a pseudo-random walk, biased downwards. The random walk is an approximation of diffusion, which is one of the mechanism for the transport of coral-nutrients. When a sphere contacts a node, that ‘nourished’ node grows away from its neighbors perpendicular to the mesh surface. The immediate neighbors of that node also grow, but by an amount proportional to there distance from the ‘nourished’ node. This extremely mechanistic algorithm produces surprisingly complex forms which bear some resemblance to stony corals.  

## Heat-Sink Design Idea
The present model could be altered to simulate heat conduction and convective effects (fluid flow) around the mesh. This simulation would serve to calculate the ‘nutrient’ level at each node in the mesh. More complex rules might be developed, such as vertices releasing a chemical that impedes growth of other neighbors and growth based on velocity or temperature gradients. [Personal Communication with Roeland Merks]

## Future Directions
- Try out more distant-neighbor nutrient sharing
- Try out more complicated polyp behavior (i.e. not growing along vert normal)


More info and images [here](www.joshlopezbinder.com/works/coral-simulation)
