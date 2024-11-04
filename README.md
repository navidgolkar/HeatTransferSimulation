# HeatTransferSimulation
<p>A simulation of a chip that is connected to a base. The chip is producing heat and from 3 sides is insulated, but from the one side a wind is blowing.</p>
<p>The problem is shown in the picture below: </br><img src="images/problem.png" alt="Problem" height="400"></br>The initial values are considered accordingly and can be seen in the code.</p>
<p>The code creates 365 nodes that are numbered according to the problem description. The simulation runs and solves the Steady State and then based on that answer calculates the best time step for unsteady solution of the problem.</br>The end result of the simuation is something like this:</br><img src="images/answer1.png" alt="Answer" height="400"></br><img src="images/answer2.png" alt="Answer" height="400"></p>
