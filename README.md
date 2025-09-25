
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
- Check orientation of the hub
- Test with 8 force model
- Look into optimizer documentation 
- Test with differnt algorithms 
- Move abaqus files into folder

