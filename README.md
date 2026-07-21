<h1>UAV Optimization for Wildfire Management</h1>

<h2>Repository Structure</h2>

<h3><code>Scripts/</code></h3>

<table>
  <tr><th>Module</th><th>Description</th></tr>
  <tr><td><code>drone_geometry</code></td><td>Calculates geometric properties of the UAV.</td></tr>
  <tr><td><code>drone_mass</code></td><td>Calculates UAV mass from a <code>drone_geometry</code> object.</td></tr>
  <tr><td><code>beam_mass</code></td><td>Calculates moment of inertia for a given <code>beam_system</code> object.</td></tr>
  <tr><td><code>mission_profile</code></td><td>Generates the mission profile used to drive loading conditions.</td></tr>
  <tr><td><code>drone_power</code></td><td>Handles calculations related to the UAV's energy consumption during a mission.</td></tr>
  <tr><td><code>optimization</code></td><td>Implements the nonlinear optimization framework.</td></tr>
</table>

<h3><code>Scripts/Finite_Elements/</code></h3>
<table>
  <tr><th>Module</th><th>Description</th></tr>
  <tr><td><code>drone_fea</code></td><td>Initializes the FEA model for a given <code>drone_geometry</code> object.</td></tr>
  <tr><td><code>beam_system</code></td><td>Core FEA framework — beam element formulation and structural solve.</td></tr>
  <tr><td><code>cross_section_properties</code></td><td>Calculates cross-sectional properties for a given beam cross-section.</td></tr>
  <tr><td><code>material_properties</code></td><td>Dictionaries of material properties used throughout the model.</td></tr>
</table>

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
