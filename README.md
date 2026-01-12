
<h3> Units </h3>

- Weight: N/mm^3
- Mass: N*s^2/mm = 1000kg = 1Mg = Tonnes
- Pressure: MPa = N/mm^2
- Length: mm
- Force: N
- Time: s
- Power: mW
- Current: mA 
- Energy: mJ

<h3> To Do </h3>
<h5> As of 11/5/24 </h5>

# TO DO:
- C rating constraint (10C)
- Fatigue Constraint (attempt)
- composite fatigue (brittle fracture) (treat as black aluminum?)
- cycle counting? Max principle stress
- Forcing function (two triangle)
- Optimal control theory ()
- python library cycle counting
- (9915) yale logo

- Read Liu paper  
- Cycle Counting library
- Office space (115) Get desk 
- Lab 151

- Fix forcing function (probably sign issue)
- Add prelift force (rev up)
- Therory of vibrations with applications thomson 3rd page 63

- Landing mission profile (forcing function goes to 0) (done)
- create logging file for optimization
- create functions for optimization problem
- FBD 
- Create plots for forcing function and stress response mid flight

- Create I beam cross-section
- (b) * 2/3 width I beam to represent hub
- Add Y-symm boundary condition
- Remove old BC if overwritten
- Fix woehler_curve_data
- structure project better

TO DO 8/14:
- Create plot of UAV section with labeled nodes and elements
- Fix woehler_curve_data
- structure project better
- scipy.optimize.minimize to run all constraint
- prepare for turning constraints on and off 
- create functions for optimization
- Start formulating thesis (write middle first)
- Grid on all plots
- n-hub elements for tapered cross section?

8/21:
- MDDO 
- Formulate objective and constraint
- NLP formulation 
- Create natural frequency constraint (will be active)
- First natural frequency > forcing function * FOS 
- forcing function * FOS  / First natural frequency < 1
- nallow = 1/FOS 
- 3/8in allowable
- FOS > FOS allowable
- eric sandrin gentic algorithms

- 8/28: 
- Make the test function work for the optimizer
- Take FOS into account

- 9/4:
- Work on energy constraint (done)
- Format constraints in terms of FOS
- Check energy constraint fomula (done)
-  put design variables in the numerator if possible 

- Fixed Energy constraint
- More detailed test print statements
- made hub height a parameter

9/11:
- Use harmonic solver instead of natural frequency (done)
- "allowable_disp_to_amplistude_ratio" (done)
- U_mag/N (done)
- Normalize the design variables (done)
- Test upper and lower bounds on design variables (done)
- Outline thesis
  - write constraints in english (write in terms of FOS (1.2))
- figures

- find abaqus to determine stress state in strut (done

9/25 
- Check orientation of the hub (done)
- Test with 8 force model
- Look into optimizer documentation 
- Test with differnt algorithms 
- Move abaqus files into folder

10/2
Questions:
How in depth should I get with my specific implmentation in the thesis. Should I talk about how my program works? no.

Change parameters to check results (done)
  experiments with optimizer (Make stress active) (focus on vehicle)
Make plots (turn off dark mode) (starting)
Read thesis (done)
Move abaqus files into folder (done)
Look at VT template for latex thesis (grad school site) (done)

10/9:

Todo make deesign variablees active on inactive from input
Seperate force functions outside
visualize design space
Consider using it practically 
Test the full model no sector 
Differeentiate symetric and asymmetric modes

To Do:

Plot Design Space 
  - Check frequency constraint (done)
  - Plot frequency as Z (done)

Test the full model (done)
  - Check natural frequencies 

Allow to set design variables to parameters. (done)

Start writing FEA section of thesis
  - Annotate plots in powerpoint

Find Exact Deadline
  - Create plan for thesis

Thesis Deadline
  10/23 (Finish FEA section)
  10/29 (Finish Drone Geometry section)
  11/6 (Finish Main content)
  11/14 (Results)
  11/21 (Finish thesis draft)
  11/28 (Revision)
  12/1/25 (Finish thesis final draft)
  12/3/25 (Minimum)

To Do:
  Save opt object to file
  Make CSV files for results (done)
  Dataframes for results (done)
  Change name of frequncy constraint on plots (done)
  Use different line styles
  
  Test points around solution to test sensitivity (done))

Refine outline 

- Talk about the specific implementation of the UAV FEA
  - Archetecture dont rederive 
  - 
  - 564 citation (article)
  - make sure it is reeproducable
  - Pngs tables
  - Color tables
  - Highlight Sector model
  - Addreess fact that I can not fit battery in reality 
  - mass + I1 + I2 + I3

- Date 10/30
  - Find justification for mission profile numbers
  - rand's document (dp 10 to red line and back)
  - Maximum: 12m/s Limit Cruise speed
  - 5600m ( distance across palisades fire)
  - Cite nyt and google earth
  - 100m treetop height
  - 1500ft (487ft) https://apps.nationalmap.gov/viewer/
  - https://www.arcgis.com/home/item.html?id=ba0631e5417f4c5f8c02a834d121bd81
  - cite the email 
  - cite nick's thesis for orginal UAV design
  - Sample points for displacement constraint to validate solution 
  - mass + I1 + I2 + I3 (working on it)
  - Switch to overleaf (meh)
  - sample design space (done)

- 11/6 
  - Finish Text Body
    - Citations
      - cite the email with dragonplate material propertiess
      - Cite nyt article? 
      - https://www.arcgis.com/home/item.html?id=ba0631e5417f4c5f8c02a834d121bd81
      - https://apps.nationalmap.gov/viewer/
      - cite nick's thesis for orginal UAV design
    - 1825ft 556.26m
    - State in thesis that I can create a more complex mission profile, uav beams. (not done) 
    - I assume a flat earth
      - Geodesic' 
      -  talk about deflection over span for justiciation. 
      - literature review for justification (other people did it, so i can do it )

    - Sum of rectangles of for I beam mass_momnt of inertia. High priority. 

    - sum moments of inertia with parallel axis theorem to get global moment of inertia of FEA model to use for objective




11/13
- Discuss Stress Failure Criterion: Smith Watson Tepper? Nu-Daniel?
  

- Change the damage calculation to actually takee the maximum bending streess




12/3
- Fix Drone 8 symm model so mass moment works 


- Test more initial guesses


- 1/8

Go over architecture diagrams? 

What to include in results? 
- Optimization Results
- Design Space Diagrams
- Sensitivity Analysis


- Stress state of results? 
- Deflection 
- explore the optimal solution
- Natural freuency anaysis of optiized results?  
- full model analysis 
- Does the result make sense??? From an engineering standpoint. Does it lead to a detailed design. 
- How much battery left? How much payload? Flight time estimate. Ect. ect. ect. 
- started 5%, tightened the constraint, demonstrate that talk about it
- draw.io

- Chapter 3 send to DR West by monday 
Find a book. ASME Journal paper. Metric paper. In the notes I sent you. ME Department formatting.


Show early results Next week
SET UP RESULTS Next week 

- Febrary 20th
feburary 8th

january 28th


To Do:
- Finish BeamElement Formulation 
- Stress calculation 
- Create Achitecture diagrams 
- Picture of segment model
- (Precicely define it (label nodes))
- Results
- Design Space Diagrams
- Optimization Results
- Sensitivity Analysis
- 
- Free body Diagram
- Show 
- Lit Review
- Appendicies
- Introduction